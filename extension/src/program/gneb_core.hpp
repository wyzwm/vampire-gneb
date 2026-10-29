// GNEB extension for VAMPIRE 7.0. Research/development version, 2026-10-29.
// GPL-2.0-or-later. Algorithm: Bessarab et al., CPC 196 (2015), 335-347.
// Tangent/force construction informed by Spirit (MIT, Gideon Mueller et al.).
// See GNEB_NOTICES.md and the retained Spirit LICENSE.txt for attribution.
#ifndef VAMPIRE_GNEB_CORE_HPP
#define VAMPIRE_GNEB_CORE_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace vgneb {
struct Vec {
    double x,y,z;
    Vec(double a=0,double b=0,double c=0):x(a),y(b),z(c){}
    Vec operator+(Vec b) const{return Vec(x+b.x,y+b.y,z+b.z);}
    Vec operator-(Vec b) const{return Vec(x-b.x,y-b.y,z-b.z);}
    Vec operator*(double a) const{return Vec(x*a,y*a,z*a);}
    Vec operator/(double a) const{return *this*(1/a);}
    Vec& operator+=(Vec b){x+=b.x;y+=b.y;z+=b.z;return *this;}
};
inline double dot(Vec a,Vec b){return a.x*b.x+a.y*b.y+a.z*b.z;}
inline Vec cross(Vec a,Vec b){return Vec(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);}
inline double norm(Vec a){return std::sqrt(dot(a,a));}
inline bool finite(Vec a){return std::isfinite(a.x)&&std::isfinite(a.y)&&std::isfinite(a.z);}
inline Vec unit(Vec a){double n=norm(a);if(!finite(a)||n<1e-15)throw std::runtime_error("Invalid zero/nonfinite spin");return a/n;}
inline Vec project(Vec v,Vec spin){return v-spin*dot(v,spin);}
typedef std::vector<Vec> Image;
typedef std::vector<Image> Chain;
typedef std::function<double(const Image&,Image&)> Evaluator; // E [meV], negative Cartesian gradient [meV]

inline double inner(const Image&a,const Image&b){
    if(a.size()!=b.size())throw std::runtime_error("Image size mismatch");
    double sum=0;for(std::size_t i=0;i<a.size();++i)sum+=dot(a[i],b[i]);return sum;
}
inline double field_norm(const Image&a){return std::sqrt(inner(a,a));}
inline double max_norm(const Image&a){double out=0;for(const Vec&v:a){if(!finite(v))throw std::runtime_error("Nonfinite field");out=std::max(out,norm(v));}return out;}
inline void validate_image(const Image&a){
    if(a.empty())throw std::runtime_error("Empty image");
    for(const Vec&s:a)if(!finite(s)||std::abs(norm(s)-1)>1e-8)throw std::runtime_error("Nonunit/nonfinite image");
}
inline double distance(const Image&a,const Image&b){
    if(a.size()!=b.size())throw std::runtime_error("Image size mismatch");
    double sum=0;for(std::size_t i=0;i<a.size();++i){double t=std::atan2(norm(cross(a[i],b[i])),dot(a[i],b[i]));sum+=t*t;}return std::sqrt(sum);
}
inline Vec exponential(Vec spin,Vec tangent){
    tangent=project(tangent,spin);double theta=norm(tangent);
    if(theta<1e-14)return unit(spin+tangent);
    return unit(spin*std::cos(theta)+tangent*(std::sin(theta)/theta));
}
inline Vec transport(Vec v,Vec oldspin,Vec newspin){
    double denom=1+dot(oldspin,newspin);
    if(denom<1e-10)throw std::runtime_error("Antipodal parallel transport is ambiguous");
    return project(v-(oldspin+newspin)*(dot(v,newspin)/denom),newspin);
}
inline Vec slerp(Vec a,Vec b,double fraction){
    double cosine=std::max(-1.0,std::min(1.0,dot(a,b)));
    if(cosine<-1+1e-10)throw std::runtime_error("Antipodal interpolation requires explicit intermediate image");
    double angle=std::atan2(norm(cross(a,b)),cosine);
    if(angle<1e-10)return unit(a*(1-fraction)+b*fraction);
    return unit((a*std::sin((1-fraction)*angle)+b*std::sin(fraction*angle))/std::sin(angle));
}

// Energy-weighted secants, projected on the product of spin tangent planes and
// normalized over all 3N components (not independently spin by spin).
inline Image tangent(const Image&minus,const Image&mid,const Image&plus,double em,double e,double ep){
    if(minus.size()!=mid.size()||plus.size()!=mid.size())throw std::runtime_error("Image size mismatch");
    Image t(mid.size());double wp=1,wm=1;
    if(ep>e&&e>em){wp=1;wm=0;}
    else if(ep<e&&e<em){wp=0;wm=1;}
    else if((ep<e&&em<e)||(ep>e&&em>e)){
        double hi=std::max(std::abs(ep-e),std::abs(em-e)),lo=std::min(std::abs(ep-e),std::abs(em-e));
        wp=ep>em?hi:lo;wm=ep>em?lo:hi;
    }
    for(std::size_t i=0;i<t.size();++i)t[i]=project((plus[i]-mid[i])*wp+(mid[i]-minus[i])*wm,mid[i]);
    double n=field_norm(t);
    if(n<1e-14)throw std::runtime_error("Degenerate path tangent; supply a nonfolded, resolved chain");
    for(Vec&v:t)v=v/n;
    return t;
}

struct Options {
    double spring=1,force_tolerance=1e-5,endpoint_tolerance=1e-5;
    double initial_step=.05,max_step=.25,max_rotation=.1;
    int max_iterations=20000,climb_after=100,report_every=50;
    bool climbing=true;
};
struct State {
    std::vector<double> energies,rx,physical_torque;
    Chain forces;
    double max_force=0,endpoint_max_torque=0;
    int climbing_index=-1,iteration=0;
    bool converged=false;
};
inline State evaluate(const Chain&images,const Evaluator&eval,const Options&o,bool climb){
    if(images.size()<3)throw std::runtime_error("GNEB needs two endpoints and at least one internal image");
    State r;std::size_t n=images.size();r.energies.resize(n);r.rx.assign(n,0);r.forces.resize(n);r.physical_torque.resize(n);
    Chain trueforce(n);
    for(std::size_t k=0;k<n;++k){
        validate_image(images[k]);if(images[k].size()!=images[0].size())throw std::runtime_error("Unequal image sizes");
        r.energies[k]=eval(images[k],trueforce[k]);
        if(!std::isfinite(r.energies[k])||trueforce[k].size()!=images[k].size())throw std::runtime_error("Invalid evaluator result");
        for(std::size_t i=0;i<images[k].size();++i)trueforce[k][i]=project(trueforce[k][i],images[k][i]);
        r.physical_torque[k]=max_norm(trueforce[k]);r.forces[k].assign(images[k].size(),Vec());
        if(k){
            for(std::size_t i=0;i<images[k].size();++i)if(dot(images[k-1][i],images[k][i])<-1+1e-10)throw std::runtime_error("Antipodal adjacent spins need an additional path image");
            double d=distance(images[k-1],images[k]);if(d<1e-10)throw std::runtime_error("Duplicate adjacent images");r.rx[k]=r.rx[k-1]+d;
        }
    }
    r.endpoint_max_torque=std::max(r.physical_torque.front(),r.physical_torque.back());
    if(climb)r.climbing_index=int(std::max_element(r.energies.begin()+1,r.energies.end()-1)-r.energies.begin());
    for(std::size_t k=1;k+1<n;++k){
        Image t=tangent(images[k-1],images[k],images[k+1],r.energies[k-1],r.energies[k],r.energies[k+1]);
        double component=inner(trueforce[k],t);
        double spring=o.spring*((r.rx[k+1]-r.rx[k])-(r.rx[k]-r.rx[k-1]));
        for(std::size_t i=0;i<t.size();++i){
            r.forces[k][i]=(int(k)==r.climbing_index)?trueforce[k][i]-t[i]*(2*component):trueforce[k][i]+t[i]*(spring-component);
            if(std::abs(dot(r.forces[k][i],images[k][i]))>1e-8*(1+norm(r.forces[k][i])))throw std::runtime_error("GNEB force left tangent space");
        }
        r.max_force=std::max(r.max_force,max_norm(r.forces[k]));
    }
    return r;
}

typedef std::function<void(const Chain&,const State&)> Reporter;
inline void validate_options(const Options&o){
    if(!std::isfinite(o.spring)||o.spring<=0||!std::isfinite(o.force_tolerance)||o.force_tolerance<=0||!std::isfinite(o.endpoint_tolerance)||o.endpoint_tolerance<=0||
       !std::isfinite(o.initial_step)||o.initial_step<=0||!std::isfinite(o.max_step)||o.max_step<o.initial_step||
       !std::isfinite(o.max_rotation)||o.max_rotation<=0||o.max_rotation>.5||o.max_iterations<1||o.report_every<1||o.climb_after<0)
        throw std::runtime_error("Invalid GNEB options");
}

// Free single-image FIRE minimization. No frozen spins, fields or topology
// constraints are introduced. Stationarity alone does not certify a minimum.
inline State relax(Image&image,const Evaluator&eval,const Options&o,const Reporter&report=Reporter()){
    validate_options(o);validate_image(image);
    Image velocity(image.size()),force;State r;
    double step=o.initial_step,mixing=.1;int positive=0;
    for(int iteration=0;iteration<=o.max_iterations;++iteration){
        double energy=eval(image,force);
        if(!std::isfinite(energy)||force.size()!=image.size())throw std::runtime_error("Invalid relaxation evaluator result");
        for(std::size_t i=0;i<image.size();++i)force[i]=project(force[i],image[i]);
        r.iteration=iteration;r.max_force=max_norm(force);r.endpoint_max_torque=r.max_force;
        r.energies.assign(1,energy);r.physical_torque.assign(1,r.max_force);
        r.converged=r.max_force<o.force_tolerance;
        if(report&&(iteration%o.report_every==0||r.converged||iteration==o.max_iterations))report(Chain(1,image),r);
        if(r.converged||iteration==o.max_iterations)break;
        for(std::size_t i=0;i<image.size();++i)velocity[i]+=force[i]*step;
        double power=inner(velocity,force);
        if(power>0){
            double ratio=field_norm(velocity)/std::max(field_norm(force),1e-300);
            for(std::size_t i=0;i<image.size();++i)velocity[i]=velocity[i]*(1-mixing)+force[i]*(mixing*ratio);
            if(++positive>5){step=std::min(step*1.1,o.max_step);mixing*=.99;}
        }else{
            for(Vec&v:velocity)v=Vec();
            step*=.5;positive=0;mixing=.1;
            if(step<1e-12)throw std::runtime_error("Endpoint FIRE step collapsed before convergence");
        }
        double largest=max_norm(velocity)*step,scale=largest>o.max_rotation?o.max_rotation/largest:1.;
        for(std::size_t i=0;i<image.size();++i){
            Vec old=image[i];image[i]=exponential(old,velocity[i]*(step*scale));velocity[i]=transport(velocity[i]*scale,old,image[i]);
        }
    }
    validate_image(image);return r;
}

// FIRE minimization on (S^2)^N. Its step is an optimizer parameter, NOT LLG time.
inline State optimize(Chain&images,const Evaluator&eval,const Options&o,const Reporter&report=Reporter()){
    validate_options(o);
    if(images.size()<3)throw std::runtime_error("GNEB needs at least three images");
    Chain endpoints{images.front(),images.back()},velocity=images;
    for(Image&im:velocity)for(Vec&v:im)v=Vec();
    double step=o.initial_step,mixing=.1;int positive=0,previous_ci=-2;
    State r;bool climbing_started=false;
    for(int iteration=0;iteration<=o.max_iterations;++iteration){
        if(o.climbing&&iteration>=o.climb_after)climbing_started=true;
        bool climbing=climbing_started;
        r=evaluate(images,eval,o,climbing);r.iteration=iteration;
        // A normal band can already be stationary before the warm-up ends.
        // Activate CI then; do not repeatedly halve FIRE's step on zero force.
        if(o.climbing&&!climbing&&r.max_force<o.force_tolerance){
            climbing_started=climbing=true;r=evaluate(images,eval,o,true);r.iteration=iteration;
        }
        if(r.endpoint_max_torque>o.endpoint_tolerance)throw std::runtime_error("Endpoints are not stationary; relax them separately first");
        r.converged=r.max_force<o.force_tolerance&&(!o.climbing||climbing);
        if(report&&(iteration%o.report_every==0||r.converged||iteration==o.max_iterations))report(images,r);
        if(r.converged||iteration==o.max_iterations)break;
        if(previous_ci!=r.climbing_index){for(Image&im:velocity)for(Vec&v:im)v=Vec();positive=0;step=o.initial_step;mixing=.1;previous_ci=r.climbing_index;}
        double power=0,vnorm2=0,fnorm2=0;
        for(std::size_t k=1;k+1<images.size();++k)for(std::size_t i=0;i<images[k].size();++i){
            velocity[k][i]+=r.forces[k][i]*step;
            power+=dot(velocity[k][i],r.forces[k][i]);vnorm2+=dot(velocity[k][i],velocity[k][i]);fnorm2+=dot(r.forces[k][i],r.forces[k][i]);
        }
        if(power>0){
            double ratio=std::sqrt(vnorm2/std::max(fnorm2,1e-300));
            for(std::size_t k=1;k+1<images.size();++k)for(std::size_t i=0;i<images[k].size();++i)velocity[k][i]=velocity[k][i]*(1-mixing)+r.forces[k][i]*(mixing*ratio);
            if(++positive>5){step=std::min(step*1.1,o.max_step);mixing*=.99;}
        }else{
            for(Image&im:velocity)for(Vec&v:im)v=Vec();
            step*=.5;positive=0;mixing=.1;
            if(step<1e-12)throw std::runtime_error("FIRE step collapsed before convergence");
        }
        double largest=0;for(std::size_t k=1;k+1<images.size();++k)largest=std::max(largest,max_norm(velocity[k])*step);
        double scale=largest>o.max_rotation?o.max_rotation/largest:1.;
        for(std::size_t k=1;k+1<images.size();++k)for(std::size_t i=0;i<images[k].size();++i){
            Vec before=images[k][i];images[k][i]=exponential(before,velocity[k][i]*(step*scale));velocity[k][i]=transport(velocity[k][i]*scale,before,images[k][i]);
        }
    }
    if(distance(endpoints[0],images.front())>1e-12||distance(endpoints[1],images.back())>1e-12)throw std::runtime_error("Endpoint changed");
    return r;
}
} // namespace vgneb
#endif
