// Public benchmark adaptation of Spirit's core/test/test_solvers.cpp and
// core/test/input/solvers.cfg, commit e82250d3b14411c2c2fa292d143f13e3e111ad8c.
// Spirit copyright (c) 2015 Gideon Mueller and collaborators, MIT.
// This independent evaluator/fixture generator is GPL-2.0-or-later.
#include "../extension/src/program/gneb_core.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace vgneb;
const int side=16;
const double muB_spirit=.057883817555;

double square_model(const Image&s,Image&f){
    f.assign(s.size(),Vec(0,0,50*muB_spirit));double e=0;
    for(int x=0;x<side;++x)for(int y=0;y<side;++y){
        int i=x*side+y;e-=50*muB_spirit*s[i].z;
        for(int axis=0;axis<2;++axis){
            int j=axis==0?((x+1)%side)*side+y:x*side+(y+1)%side;
            Vec D=axis==0?Vec(6,0,0):Vec(0,6,0);
            e-=10*dot(s[i],s[j])+dot(D,cross(s[i],s[j]));
            f[i]+=s[j]*10+cross(s[j],D);f[j]+=s[i]*10+cross(D,s[i]);
        }
    }
    return e;
}
Image initial_skyrmion(){
    Image s(side*side,Vec(0,0,1));
    for(int x=0;x<side;++x)for(int y=0;y<side;++y){
        double dx=x-7.5,dy=y-7.5,r=std::hypot(dx,dy);
        if(r<=5){double theta=std::acos(-1.)*(1-r/5),phi=std::atan2(dy,dx)-std::acos(-1.)/2;
            s[x*side+y]=Vec(-std::sin(theta)*std::cos(phi),-std::sin(theta)*std::sin(phi),std::cos(theta));}
    }return s;
}
Image relax(Image s){
    Image v(s.size()),force;double dt=.005,mix=.1;int positive=0;
    for(int iteration=0;iteration<200000;++iteration){
        double e=square_model(s,force);for(std::size_t i=0;i<s.size();++i)force[i]=project(force[i],s[i]);
        double maxf=max_norm(force);
        if(iteration%5000==0)std::cerr<<"minimize "<<iteration<<" E "<<std::setprecision(14)<<e<<" force "<<maxf<<'\n';
        if(maxf<1e-8){std::cerr<<"minimum "<<iteration<<" E "<<std::setprecision(14)<<e<<" force "<<maxf<<'\n';return s;}
        for(std::size_t i=0;i<s.size();++i)v[i]+=force[i]*dt;
        double power=inner(v,force);
        if(power>0){double ratio=field_norm(v)/std::max(field_norm(force),1e-300);for(std::size_t i=0;i<s.size();++i)v[i]=v[i]*(1-mix)+force[i]*(mix*ratio);if(++positive>5){dt=std::min(dt*1.1,.05);mix*=.99;}}
        else{for(Vec&a:v)a=Vec();dt*=.5;positive=0;mix=.1;}
        double maximum=max_norm(v)*dt;double scale=maximum>.05?.05/maximum:1;
        for(std::size_t i=0;i<s.size();++i){Vec old=s[i];s[i]=exponential(old,v[i]*(dt*scale));v[i]=transport(v[i]*scale,old,s[i]);}
        if(dt<1e-13)throw std::runtime_error("Endpoint minimization stalled");
    }throw std::runtime_error("Endpoint minimization iteration limit");
}
void save(const std::string&name,const Image&s){
    if(std::ifstream(name.c_str()))throw std::runtime_error("Output exists; choose a fresh directory");
    std::ofstream out(name.c_str());if(!out)throw std::runtime_error("Cannot write image");
    out<<"GNEB_SPINS_V1\n"<<s.size()<<'\n'<<std::setprecision(17);
    for(int x=0;x<side;++x)for(int y=0;y<side;++y){const Vec&a=s[x*side+y];out<<x+.5<<' '<<y+.5<<" 0.5 "<<a.x<<' '<<a.y<<' '<<a.z<<'\n';}
}
int main(int argc,char**argv){try{
    if(std::ifstream("chain.txt"))throw std::runtime_error("chain.txt exists; choose a fresh directory");
    for(int k=0;k<9;++k)if(std::ifstream(("image_"+std::to_string(k)+".txt").c_str())||std::ifstream(("final_"+std::to_string(k)+".txt").c_str()))throw std::runtime_error("Benchmark output exists");
    Image start=relax(initial_skyrmion()),end(side*side,Vec(0,0,1)),force;
    double minimum=square_model(start,force),uniform=square_model(end,force);
    std::cout<<std::setprecision(15)<<"minimum_energy_meV="<<minimum<<" uniform_energy_meV="<<uniform<<'\n';
    if(std::abs(minimum-(-5849.69140625))>3e-4)throw std::runtime_error("Endpoint differs from public Spirit benchmark");
    int images=9;Chain chain(images,Image(side*side));
    for(int k=0;k<images;++k)for(int i=0;i<side*side;++i)chain[k][i]=slerp(start[i],end[i],double(k)/(images-1));
    for(int k=0;k<images;++k)save("image_"+std::to_string(k)+".txt",chain[k]);
    std::ofstream list("chain.txt");for(int k=0;k<images;++k)list<<"image_"<<k<<".txt\n";list.close();
    if(argc>1&&std::string(argv[1])=="--prepare-only")return 0;
    Options o;o.initial_step=.005;o.max_step=.05;o.force_tolerance=1e-6;o.endpoint_tolerance=1e-6;o.max_iterations=200000;o.climb_after=2000;o.report_every=5000;
    State result=optimize(chain,square_model,o,[](const Chain&,const State&r){std::cerr<<"gneb "<<r.iteration<<" force "<<r.max_force<<" CI "<<r.climbing_index<<'\n';});
    double saddle=*std::max_element(result.energies.begin(),result.energies.end());
    std::cout<<"converged="<<result.converged<<" saddle_energy_meV="<<saddle<<" forward_barrier_meV="<<saddle-minimum<<" expected_saddle_meV=-5811.5244140625\n";
    for(int k=0;k<images;++k)save("final_"+std::to_string(k)+".txt",chain[k]);
    if(!result.converged||std::abs(saddle-(-5811.5244140625))>3e-4)return 2;
    return 0;
}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
