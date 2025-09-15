#ifndef BEAVER_TRILOG_HPP
#define BEAVER_TRILOG_HPP
#include "log.hpp"
#include "config.hpp"
#include <cmath>

namespace beaver {
    namespace internals::trilog{
    //universal constants
  const double ln2=0.69314718055994530942;
  const double zeta2=1.6449340668482264365;
  const double zeta3=1.2020569031595942854;
  //MiniMax numerator coefficients of log(1-x)/x
  const double P1[]={-0.9999999999999998672,2.7601168443679114620,-2.7923621883637091036,1.2497503713725827027,-0.23173044836039025061,0.012188323559602391476};
    //MiniMax denominator coefficients of log(1-x)/x
  const double Q1[]={1.0000000000000000000,-3.2601168443679653008,4.0890872772181186328,-2.4575883952960021967,0.71252476615980667079,-0.086169703890931582878,0.0026519629643638051268};
    //MiniMax numerator coefficients of Li3(x)/x on 10^(-4)<x<0.5
  const double P3a[]={0.99999999999999998889,-2.5225087681785685178,2.3205667855092737010,-0.93986039506026480652,0.15730234392961045107,-0.0075494230157797519206};        
  //MiniMax denominator coefficients of Li3(x)/x on 10^(-4)<x<0.5
  const double Q3a[]={1.0000000000000000000,-2.6475087681785732399,2.6144683444949039205,-1.1842380578292087625,0.24186726490232288650,-0.018222790964808075074,0.00024931477353312848423}; 
    //MiniMax numerator coefficients of Li3(x)/x on -1<x<10^(-4)
  const double P3b[]={0.99999999999999997872,-2.0281455911060640193,1.4363491768494536537,-0.42238055117778581583,0.047292383318533615224,-0.0013451800040195675605};        
  //MiniMax denominator coefficients of Li3(x)/x on -1<x<10^(-4)
  const double Q3b[]={1.0000000000000000000,-2.1531455911060558732,1.6684553387011928972,-0.56681633549854450482,0.081992683111646187121,-0.0040751460198142638061,0.000034311405527458791883}; 
    //MiniMax numerator coefficients of (-Li3(x/(x-1))-Li3(x))/x on 10^(-4)<x<0.5
  const double P3c[]={0.0000000000000000000,0.74999999999999999816,-2.1777582585501254342,2.2971231247321237024,-1.0455892723397490220,0.18551293485607298693,-0.0087864196233895762452};        
  //MiniMax denominator coefficients of (-Li3(x/(x-1))-Li3(x))/x on 10^(-4)<x<0.5
  const double Q3c[]={1.0000000000000000000,-3.9036776780668347591,6.0266936962282502098,-4.6317082415371928343,1.8223918699839630482,-0.33647047688023663835,0.023067775463027464700,-0.00028585989330444185814}; 
    }
  //Li3 kernel around x=0
  static double Li3_Taylor(double x){
        double x2=x*x;
        double x4=x2*x2;
        return x+(1.0/8+1.0/27*x)*x2+1.0/64*x4;
    }
    BEAVER_NODISCARD inline double  trilog(double x) noexcept {
    namespace  LOC=internals::trilog;
    double taylorswitch=1e-4;//switch up to which an expansion about x=0 or x=1 is used. Guarantees precission for small y.
    double invtaylorswitch=1/taylorswitch;
  //Catch non-finite input
  if (!std::isfinite(x)) return std::numeric_limits<double>::quiet_NaN();
  if(x<-invtaylorswitch){
    double l=beaver::log(-x);
    double l3=l*l*l;
    double c=-l*LOC::zeta2-1.0/6*l3;
    double xinv=1.0/x;
    double li3inv=Li3_Taylor(xinv);
    return li3inv+c;
  }else if(x<-1){
    double l=beaver::log(-x);
    double l3=l*l*l;
    double c=-l*LOC::zeta2-1.0/6*l3;
    double x2=x*x;
    double x4=x2*x2;
    double x6=x4*x2;
    double p=LOC::P3b[5]+x*LOC::P3b[4]+x2*(LOC::P3b[3]+x*LOC::P3b[2])+x4*(LOC::P3b[1]+x*LOC::P3b[0]);
    double q=LOC::Q3b[6]+x*LOC::Q3b[5]+x2*(LOC::Q3b[4]+x*LOC::Q3b[3])+x4*(LOC::Q3b[2]+x*LOC::Q3b[1])+x6*LOC::Q3b[0];
    return p/q+c;
  }else if(x==-1){
    return -0.75*LOC::zeta3;
  }else if(x<-taylorswitch){
    double x2=x*x;
    double x4=x2*x2;
    double x6=x4*x2;
    double p=x*LOC::P3b[0]+x2*(LOC::P3b[1]+x*LOC::P3b[2])+x4*(LOC::P3b[3]+x*LOC::P3b[4])+x6*LOC::P3b[5];
    double q=LOC::Q3b[0]+x*LOC::Q3b[1]+x2*(LOC::Q3b[2]+x*LOC::Q3b[3])+x4*(LOC::Q3b[4]+x*LOC::Q3b[5])+x6*LOC::Q3b[6];
    return p/q;
  }else if(x<taylorswitch){
    return Li3_Taylor(x);
  }else if(x<0.5){
    double x2=x*x;
    double x4=x2*x2;
    double x6=x4*x2;
    double p=x*LOC::P3a[0]+x2*(LOC::P3a[1]+x*LOC::P3a[2])+x4*(LOC::P3a[3]+x*LOC::P3a[4])+x6*LOC::P3a[5];
    double q=LOC::Q3a[0]+x*LOC::Q3a[1]+x2*(LOC::Q3a[2]+x*LOC::Q3a[3])+x4*(LOC::Q3a[4]+x*LOC::Q3a[5])+x6*LOC::Q3a[6];
    return p/q;
  }else if(x==0.5)
  {
    return 21.0/24*LOC::zeta3+1.0/6*LOC::ln2*LOC::ln2*LOC::ln2-0.5*LOC::zeta2*LOC::ln2;
  }else if(x<1-taylorswitch){
    double y=1-x;
    double y2=y*y;
    double y4=y2*y2;
    double y6=y2*y4;
    double ly=beaver::log(y);
    //calculate l=log(x) with MiniMax (argument is guaranteed to be between 0.5 and 1)
    //double l=beaver::log(x); //is slower
    double pl=y*LOC::P1[0]+y2*std::fma(LOC::P1[2],y,LOC::P1[1])+y4*std::fma(LOC::P1[4],y,LOC::P1[3])+y6*LOC::P1[5];
    double ql=LOC::Q1[0]+y*LOC::Q1[1]+y2*std::fma(LOC::Q1[3],y,LOC::Q1[2])+y4*std::fma(LOC::Q1[5],y,LOC::Q1[4])+y6*LOC::Q1[6];
    double l=pl/ql;
    double l2=l*l;
    double l3=l2*l;
    double c=LOC::zeta3+l*LOC::zeta2-1.0/2*l2*ly+1.0/6*l3;
    double p=y2*(LOC::P3c[1]+y*LOC::P3c[2])+y4*(LOC::P3c[3]+y*LOC::P3c[4])+y6*(LOC::P3c[5]+y*LOC::P3c[6]);
    double q=LOC::Q3c[0]+y*LOC::Q3c[1]+y2*(LOC::Q3c[2]+y*LOC::Q3c[3])+y4*(LOC::Q3c[4]+y*LOC::Q3c[5])+y6*(LOC::Q3c[6]+y*LOC::Q3c[7]);
    return p/q+c;
  } else if(x<1){
    double y=1-x;
    double y2=y*y;
    double y4=y2*y2;
    double ly=beaver::log(y);
    double taylor=LOC::zeta3-y*LOC::zeta2+y2*(0.75-LOC::zeta2/2+y*(7.0/12-1.0/3*LOC::zeta2))+y4*(131.0/288-1.0/4*LOC::zeta2);
    double logtaylor=y2*(-0.5-0.5*y)-11.0/24*y4;
    return taylor+ly*logtaylor;
  }else if(x==1){
    return LOC::zeta3;
  }else{//Argument on branch-cut
    return std::numeric_limits<double>::quiet_NaN(); 
  }
}

}
#endif