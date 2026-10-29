// Isolated experimental VAMPIRE 7.0 GNEB program. See GNEB_NOTICES.md.
#include "gneb.hpp"
#include "gneb_core.hpp"
#include "atoms.hpp"
#include "anisotropy.hpp"
#include "exchange.hpp"
#include "material.hpp"
#include "dipole.hpp"
#include "gpu.hpp"
#include "micromagnetic.hpp"
#include "sim.hpp"
#include "vmpi.hpp"
#include "../simulate/internal.hpp"
#include "../exchange/internal.hpp"
#include "../anisotropy/internal.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <tuple>
#include <cstdlib>

namespace gneb {
namespace {
vgneb::Options options;
std::string chain_file,output_prefix="gneb",mode="path";
int threads=1,checkpoint_every=0;
const double mev_j=1.602176634e-22;
double normalization_adjustment=0;

double number(const std::string&s){std::size_t end=0;double x=std::stod(s,&end);if(end!=s.size()||!std::isfinite(x))throw std::runtime_error("Invalid numeric GNEB value: "+s);return x;}
int integer(const std::string&s){double x=number(s);if(x!=std::floor(x)||x<0||x>100000000)throw std::runtime_error("Invalid GNEB integer");return int(x);}
bool boolean(const std::string&s){if(s=="true")return true;if(s=="false")return false;throw std::runtime_error("Expected true or false");}
bool exists(const std::string&s){std::ifstream f(s.c_str());return bool(f);}
typedef std::tuple<long long,long long,long long> Coord;
Coord key(double x,double y,double z){return std::make_tuple(std::llround(x*1e6),std::llround(y*1e6),std::llround(z*1e6));}
std::string trim(std::string s){s=s.substr(0,s.find('#'));std::size_t a=s.find_first_not_of(" \t\r\n\"");if(a==std::string::npos)return "";return s.substr(a,s.find_last_not_of(" \t\r\n\"")-a+1);}

void guard_model(){
#ifdef MPICF
    throw std::runtime_error("GNEB is a single-process program. MPI domain/image parallelism has not been implemented or validated.");
#endif
    if(sim::temperature!=0)throw std::runtime_error("GNEB requires sim:temperature=0 (no stochastic LLG)");
    if(sim::load_checkpoint_flag)throw std::runtime_error("Supply explicit coordinate-labelled GNEB images; checkpoint loading is unsupported");
    // The legacy public four_spin flag defaults true even for zero coupling;
    // enable_fourspin is the actual user-input activation flag.
    if(dipole::activated||gpu::acceleration||sim::lagrange_multiplier||exchange::biquadratic||exchange::internal::enable_fourspin||
       sim::internal::enable_spin_torque_fields||sim::internal::enable_vcma_fields||micromagnetic::discretisation_type!=0)
        throw std::runtime_error("Unsupported GNEB model: dipoles/GPU/constraints/higher exchange/nonconservative torques/micromagnetics");
    if(atoms::num_atoms<=0||atoms::num_non_magnetic_atoms)throw std::runtime_error("GNEB requires a nonempty all-magnetic atomistic system");
    for(int i=0;i<atoms::num_atoms;++i){
        const double mu=mp::material[atoms::type_array[i]].mu_s_SI;
        if(!std::isfinite(mu)||mu<=0)throw std::runtime_error("Invalid magnetic moment");
    }
    if(threads<1||threads>256)throw std::runtime_error("gneb:threads must be between 1 and 256");
#ifndef _OPENMP
    if(threads!=1)throw std::runtime_error("This build has no OpenMP support; rebuild with --openmp or use threads=1");
#else
    omp_set_dynamic(0);
    if(threads>omp_get_num_procs())throw std::runtime_error("Requested threads exceed available CPUs");
#endif
    if(threads>1){
        // The optional parallel native kernel is deliberately limited to the
        // audited range-local bilinear exchange + K2 + static Zeeman model.
        using namespace anisotropy::internal;
        if(enable_uniaxial_fourth_order||enable_biaxial_fourth_order_simple||enable_uniaxial_sixth_order||
           enable_cubic_fourth_order||enable_cubic_sixth_order||enable_cubic_fourth_order_rotation||
           enable_fourth_order_rotational||enable_triaxial_anisotropy_rotated||enable_triaxial_fourth_order_rotated||
           enable_triaxial_anisotropy||enable_triaxial_fourth_order||enable_neel_anisotropy||enable_lattice_anisotropy||enable_random_anisotropy||
           !exchange::internal::four_spin_neighbour_list_array_l.empty())
            throw std::runtime_error("OpenMP GNEB supports only bilinear exchange, second-order uniaxial anisotropy and static Zeeman");
    }
}

// Reuse the actual initialized VAMPIRE neighbour lists, exchange tensors,
// anisotropy energy/field implementations and native magnetic moments.
double evaluate_native(const vgneb::Image&image,vgneb::Image&force){
    if(image.size()!=std::size_t(atoms::num_atoms))throw std::runtime_error("Spin count differs from initialized VAMPIRE geometry");
    for(int i=0;i<atoms::num_atoms;++i){atoms::x_spin_array[i]=image[i].x;atoms::y_spin_array[i]=image[i].y;atoms::z_spin_array[i]=image[i].z;}
#ifdef _OPENMP
    #pragma omp parallel num_threads(threads) if(threads>1)
    {
        int t=omp_get_thread_num(),n=omp_get_num_threads();
        sim::calculate_spin_fields(atoms::num_atoms*t/n,atoms::num_atoms*(t+1)/n);
    }
#else
    sim::calculate_spin_fields(0,atoms::num_atoms);
#endif
    force.resize(image.size());std::vector<double> site_energy(image.size());double energy=0;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads) if(threads>1)
#endif
    for(int i=0;i<atoms::num_atoms;++i){
        int mat=atoms::type_array[i];double mu=mp::material[mat].mu_s_SI;
        vgneb::Vec field(sim::H_applied*sim::H_vec[0],sim::H_applied*sim::H_vec[1],sim::H_applied*sim::H_vec[2]);
        field+=vgneb::Vec(mp::material[mat].applied_field_unit_vector[0],mp::material[mat].applied_field_unit_vector[1],mp::material[mat].applied_field_unit_vector[2])*mp::material[mat].applied_field_strength;
        site_energy[i]=mu*(.5*exchange::single_spin_energy(i,image[i].x,image[i].y,image[i].z)
                    +anisotropy::single_spin_energy(i,mat,image[i].x,image[i].y,image[i].z,0)-vgneb::dot(image[i],field))/mev_j;
        force[i]=(vgneb::Vec(atoms::x_total_spin_field_array[i],atoms::y_total_spin_field_array[i],atoms::z_total_spin_field_array[i])+field)*(mu/mev_j);
    }
    // Identical atom-order summation for every thread count (no reduction tree).
    for(double e:site_energy)energy+=e;
    return energy;
}

vgneb::Chain read_chain(){
    if(chain_file.empty())throw std::runtime_error("Missing gneb:chain-file");
    std::map<Coord,int> indices;
    for(int i=0;i<atoms::num_atoms;++i){if(!indices.emplace(key(atoms::x_coord_array[i],atoms::y_coord_array[i],atoms::z_coord_array[i]),i).second)throw std::runtime_error("Geometry has duplicate coordinates at 1e-6 A resolution");}
    std::ifstream chain(chain_file.c_str());if(!chain)throw std::runtime_error("Cannot read GNEB chain file");
    std::string base;std::size_t slash=chain_file.find_last_of("/\\");if(slash!=std::string::npos)base=chain_file.substr(0,slash+1);
    vgneb::Chain result;std::string line;
    while(std::getline(chain,line)){
        std::string name=trim(line);if(name.empty())continue;
        std::ifstream file((base+name).c_str());if(!file)throw std::runtime_error("Cannot open chain image "+base+name);
        std::string magic;std::size_t n=0;file>>magic>>n;
        if(magic!="GNEB_SPINS_V1"||n!=std::size_t(atoms::num_atoms))throw std::runtime_error("Bad GNEB image header/spin count");
        vgneb::Image spins(n);std::vector<bool> seen(n,false);
        for(std::size_t row=0;row<n;++row){
            double x,y,z,sx,sy,sz;if(!(file>>x>>y>>z>>sx>>sy>>sz))throw std::runtime_error("Truncated GNEB image");
            if(!std::isfinite(x)||!std::isfinite(y)||!std::isfinite(z))throw std::runtime_error("Nonfinite image coordinates");
            auto it=indices.find(key(x,y,z));if(it==indices.end()||seen[it->second])throw std::runtime_error("Image coordinate mismatch or duplicate");
            int i=it->second;seen[i]=true;vgneb::Vec s(sx,sy,sz);double deviation=std::abs(vgneb::norm(s)-1);
            if(!vgneb::finite(s)||deviation>1e-5)throw std::runtime_error("Nonunit/nonfinite input spin (tolerance 1e-5)");
            normalization_adjustment=std::max(normalization_adjustment,deviation);spins[i]=vgneb::unit(s);
        }
        std::string extra;if(file>>extra)throw std::runtime_error("Extra data in GNEB spin image");
        result.push_back(spins);if(result.size()>256)throw std::runtime_error("Too many GNEB images (limit 256)");
    }
    if(result.empty())throw std::runtime_error("Empty GNEB chain");
    return result;
}

void write_image(const std::string&name,const vgneb::Image&s){
    std::ofstream f(name.c_str());if(!f)throw std::runtime_error("Cannot write "+name);
    f<<"GNEB_SPINS_V1\n"<<s.size()<<'\n'<<std::setprecision(17);
    for(std::size_t i=0;i<s.size();++i)f<<atoms::x_coord_array[i]<<' '<<atoms::y_coord_array[i]<<' '<<atoms::z_coord_array[i]<<' '<<s[i].x<<' '<<s[i].y<<' '<<s[i].z<<'\n';
}
void write_force(const std::string&name,const vgneb::Image&force){
    std::ofstream f(name.c_str());if(!f)throw std::runtime_error("Cannot write "+name);
    f<<"# x_A y_A z_A negative_dE_dSx_meV negative_dE_dSy_meV negative_dE_dSz_meV\n"<<std::setprecision(17);
    for(std::size_t i=0;i<force.size();++i)f<<atoms::x_coord_array[i]<<' '<<atoms::y_coord_array[i]<<' '<<atoms::z_coord_array[i]<<' '<<force[i].x<<' '<<force[i].y<<' '<<force[i].z<<'\n';
}

void consistency_gate(const vgneb::Image&image){
    vgneb::Image force;evaluate_native(image,force);double largest=0;
    const std::size_t samples=std::min(std::size_t(12),image.size());
    for(std::size_t j=0;j<samples;++j){
        std::size_t i=j*image.size()/samples;vgneb::Vec axis=std::abs(image[i].z)<.9?vgneb::Vec(0,0,1):vgneb::Vec(1,0,0);
        vgneb::Vec direction=vgneb::unit(vgneb::project(axis,image[i]));const double epsilon=1e-4;
        vgneb::Image plus=image,minus=image,ignored;
        plus[i]=vgneb::exponential(image[i],direction*epsilon);minus[i]=vgneb::exponential(image[i],direction*(-epsilon));
        double fd=(evaluate_native(plus,ignored)-evaluate_native(minus,ignored))/(2*epsilon);
        double expected=-vgneb::dot(force[i],direction);largest=std::max(largest,std::abs(fd-expected));
        double tolerance=1e-6+std::abs(expected)*1e-5+double(image.size())*1e-11;
        if(std::abs(fd-expected)>tolerance)throw std::runtime_error("Native energy/gradient finite-difference check failed");
    }
    evaluate_native(image,force);
    std::cout<<"GNEB native energy-gradient max discrepancy [meV] = "<<std::setprecision(12)<<largest<<std::endl;
}

void execute(){
    guard_model();vgneb::validate_options(options);
    if(mode!="path"&&mode!="evaluate"&&mode!="relax")throw std::runtime_error("gneb:mode must be path, relax or evaluate");
    if(output_prefix.empty()||output_prefix.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")!=std::string::npos)throw std::runtime_error("GNEB output-prefix must be a simple basename");
    if(exists(output_prefix+"_status.json")||exists(output_prefix+"_energies.csv")||exists(output_prefix+"_history.csv"))throw std::runtime_error("GNEB outputs already exist; use a fresh directory/prefix");
    vgneb::Chain chain=read_chain();consistency_gate(chain.front());
    if(mode=="relax"&&chain.size()!=1)throw std::runtime_error("Relax mode requires exactly one image");
    std::ofstream history((output_prefix+"_history.csv").c_str());if(!history)throw std::runtime_error("Cannot write GNEB history");
    history<<"iteration,max_GNEB_force_meV,endpoint_torque_meV,climbing_index\n"<<std::setprecision(17);
    vgneb::State state;
    if(mode=="relax"){
        state=vgneb::relax(chain.front(),evaluate_native,options,[&](const vgneb::Chain&c,const vgneb::State&s){
            history<<s.iteration<<','<<s.max_force<<','<<s.endpoint_max_torque<<",-1\n";history.flush();
            std::cout<<"GNEB relax iteration "<<s.iteration<<" energy "<<std::setprecision(17)<<s.energies[0]<<" max torque "<<s.max_force<<" meV"<<std::endl;
            if(checkpoint_every>0&&(s.iteration%checkpoint_every==0||s.converged||s.iteration==options.max_iterations))
                write_image(output_prefix+"_relax_checkpoint_"+std::to_string(s.iteration)+".txt",c.front());
        });
    }
    if(mode=="path"){
        state=vgneb::optimize(chain,evaluate_native,options,[&](const vgneb::Chain&,const vgneb::State&s){
            history<<s.iteration<<','<<s.max_force<<','<<s.endpoint_max_torque<<','<<s.climbing_index<<'\n';history.flush();
            std::cout<<"GNEB iteration "<<s.iteration<<" max force "<<s.max_force<<" meV CI "<<s.climbing_index<<std::endl;
        });
    }
    std::ofstream energies((output_prefix+"_energies.csv").c_str());if(!energies)throw std::runtime_error("Cannot write energies");
    energies<<"image,rx_rad,energy_meV,delta_energy_meV,physical_torque_meV\n"<<std::setprecision(17);
    std::vector<double> values;double rx=0,e0=0;
    for(std::size_t k=0;k<chain.size();++k){
        vgneb::Image force;double energy=evaluate_native(chain[k],force);values.push_back(energy);if(k==0)e0=energy;if(k)rx+=vgneb::distance(chain[k-1],chain[k]);
        vgneb::Image tangent=force;for(std::size_t i=0;i<force.size();++i)tangent[i]=vgneb::project(force[i],chain[k][i]);
        energies<<k<<','<<rx<<','<<energy<<','<<energy-e0<<','<<vgneb::max_norm(tangent)<<'\n';
        write_image(output_prefix+"_image_"+std::to_string(k)+".txt",chain[k]);
        if(mode=="evaluate")write_force(output_prefix+"_force_"+std::to_string(k)+".tsv",force);
    }
    double top=*std::max_element(values.begin(),values.end());
    std::ofstream status((output_prefix+"_status.json").c_str());if(!status)throw std::runtime_error("Cannot write status");
    status<<std::setprecision(17)<<"{\n  \"implementation\": \"vampire7-gneb-0.2\",\n  \"mode\": \""<<mode<<"\",\n"
          <<"  \"native_kernel_threads\": "<<threads<<",\n"
          <<"  \"images\": "<<chain.size()<<",\n  \"spins_per_image\": "<<chain.front().size()<<",\n"
          <<"  \"numerically_converged\": "<<(mode!="evaluate"&&state.converged?"true":"false")<<",\n"
          <<"  \"iterations\": "<<state.iteration<<",\n  \"max_GNEB_force_meV\": "<<state.max_force<<",\n"
          <<"  \"endpoint_torque_meV\": "<<state.endpoint_max_torque<<",\n  \"climbing_image\": "<<state.climbing_index<<",\n"
          <<"  \"path_barrier_evaluated\": "<<(mode=="path"?"true":"false")<<",\n  \"sampled_forward_barrier_meV\": ";
    if(mode=="path")status<<top-values.front();else status<<"null";
    status<<",\n  \"sampled_reverse_barrier_meV\": ";
    if(mode=="path")status<<top-values.back();else status<<"null";
    status<<",\n  \"maximum_input_norm_adjustment\": "<<normalization_adjustment<<",\n"
          <<"  \"physical_time_integration\": false,\n  \"saddle_index_or_material_stability_certified\": false\n}\n";
    status.close();energies.close();history.close();
    if(mode!="evaluate"&&!state.converged){std::cerr<<"GNEB reached iteration budget without convergence"<<std::endl;std::exit(2);}
}
} // private

bool match_input_parameter(const std::string&key,const std::string&word,const std::string&value,const std::string&unit,int){
    if(key!="gneb")return false;
    try{
        if(!unit.empty())throw std::runtime_error("GNEB parameters use explicit native meV/radian optimizer units, no unit suffix");
        if(word=="chain-file")chain_file=trim(value);
        else if(word=="output-prefix")output_prefix=trim(value);
        else if(word=="mode")mode=value;
        else if(word=="threads")threads=integer(value);
        else if(word=="checkpoint-every")checkpoint_every=integer(value);
        else if(word=="spring-constant")options.spring=number(value);
        else if(word=="force-tolerance")options.force_tolerance=number(value);
        else if(word=="endpoint-tolerance")options.endpoint_tolerance=number(value);
        else if(word=="max-iterations")options.max_iterations=integer(value);
        else if(word=="climbing")options.climbing=boolean(value);
        else if(word=="climb-after")options.climb_after=integer(value);
        else if(word=="report-every")options.report_every=integer(value);
        else if(word=="initial-step")options.initial_step=number(value);
        else if(word=="maximum-step")options.max_step=number(value);
        else if(word=="maximum-rotation")options.max_rotation=number(value);
        else throw std::runtime_error("Unknown gneb parameter: "+word);
        return true;
    }catch(const std::exception&e){std::cerr<<"GNEB input error: "<<e.what()<<std::endl;std::exit(1);}
}
void run(){try{execute();}catch(const std::exception&e){std::cerr<<"GNEB error: "<<e.what()<<std::endl;std::exit(1);}}
}
