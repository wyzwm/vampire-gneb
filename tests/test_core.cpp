#include "../extension/src/program/gneb_core.hpp"
#include <iostream>
#include <iomanip>
using namespace vgneb;
void require(bool x,const char*m){if(!x)throw std::runtime_error(m);}
template<class F>void rejects(F f){bool thrown=false;try{f();}catch(const std::exception&){thrown=true;}require(thrown,"Expected rejection");}
Chain initial(int images,int spins){
    Chain c(images,Image(spins));
    for(int k=0;k<images;++k){
        double f=std::pow(double(k)/(images-1),1.27),theta=3.141592653589793*(1-f),phi=.35*std::sin(theta);
        for(int i=0;i<spins;++i)c[k][i]=Vec(std::sin(theta)*std::cos(phi),std::sin(theta)*std::sin(phi),std::cos(theta));
    }return c;
}
double model(const Image&s,Image&force){
    force.assign(s.size(),Vec());double e=0;
    for(std::size_t i=0;i<s.size();++i){e+=-s[i].z*s[i].z+.3*s[i].y*s[i].y-.2*s[i].z;force[i]=Vec(0,-.6*s[i].y,2*s[i].z+.2);}
    if(s.size()==2){e-=4*dot(s[0],s[1]);force[0]+=s[1]*4;force[1]+=s[0]*4;}
    return e;
}
int main(){try{
    require(std::abs(distance(Image{Vec(0,0,1)},Image{Vec(1,0,0)})-std::acos(0.))<1e-14,"geodesic distance");
    require(std::abs(distance(Image{Vec(0,0,1)},Image{Vec(0,0,-1)})-std::acos(-1.))<1e-14,"antipodal distance");
    rejects([](){slerp(Vec(0,0,1),Vec(0,0,-1),.5);});
    Vec a=unit(Vec(.2,.4,.7)),v=project(Vec(.3,-.6,.1),a),b=exponential(a,v*.1),w=transport(v,a,b);
    require(std::abs(norm(b)-1)<1e-14&&std::abs(dot(w,b))<1e-14&&std::abs(norm(w)-norm(v))<1e-12,"sphere and transport");
    Options o;o.force_tolerance=1e-8;o.endpoint_tolerance=1e-8;o.max_iterations=30000;o.report_every=1000;
    rejects([&](){Chain z;optimize(z,model,o);});
    rejects([&](){Chain z=initial(7,1);z[3]=z[2];evaluate(z,model,o,false);});
    rejects([&](){Chain z=initial(7,1);z[0][0]=unit(Vec(.2,0,-1));optimize(z,model,o);});
    std::cout<<std::setprecision(15);
    for(int spins: {1,2})for(int images: {9,17}){
        Chain c=initial(images,spins);State result=optimize(c,model,o);
        double barrier=*std::max_element(result.energies.begin(),result.energies.end())-result.energies.front();
        require(result.converged,"Analytic path failed convergence");
        require(std::abs(barrier-.81*spins)<1e-6,"Wrong analytic forward barrier");
        double reverse=*std::max_element(result.energies.begin(),result.energies.end())-result.energies.back();
        require(std::abs(reverse-1.21*spins)<1e-6,"Wrong analytic reverse barrier");
        require(result.climbing_index>0&&result.climbing_index<images-1,"CI index");
        for(const Image&i:c)validate_image(i);
        std::cout<<"PASS spins="<<spins<<" images="<<images<<" iterations="<<result.iteration<<" barrier="<<barrier<<" reverse="<<reverse<<" max_force="<<result.max_force<<'\n';
    }
    Options short_run=o;short_run.max_iterations=1;Chain short_chain=initial(9,1);require(!optimize(short_chain,model,short_run).converged,"Iteration limit incorrectly accepted");
    Chain already_on_path(9,Image(1));
    for(int k=0;k<9;++k){double t=std::acos(-1.)*(1-k/8.);already_on_path[k][0]=Vec(std::sin(t),0,std::cos(t));}
    State ready=optimize(already_on_path,model,o);
    require(ready.converged&&std::abs(*std::max_element(ready.energies.begin(),ready.energies.end())-ready.energies.front()-.81)<1e-6,"CI activation from already-converged normal band");
    rejects([&](){Chain z=initial(7,1);z[3][0]=z[2][0]*(-1);evaluate(z,model,o,false);});
    // Independent Cartesian finite differences of the coupled analytic test evaluator.
    Image s{unit(Vec(.3,.7,.5)),unit(Vec(-.6,.2,.7))},f;model(s,f);
    for(int i=0;i<2;++i)for(int axis=0;axis<3;++axis){Image p=s,m=s,unused;double *pc=axis==0?&p[i].x:axis==1?&p[i].y:&p[i].z;double *mc=axis==0?&m[i].x:axis==1?&m[i].y:&m[i].z;*pc+=1e-6;*mc-=1e-6;double fd=(model(p,unused)-model(m,unused))/2e-6;double fi=axis==0?f[i].x:axis==1?f[i].y:f[i].z;require(std::abs(fd+fi)<1e-8,"analytic gradient mismatch");}
    std::cout<<"PASS geometry, projection, endpoint/duplicate/antipodal guards, iteration-limit and independent gradients\n";
    return 0;
}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
