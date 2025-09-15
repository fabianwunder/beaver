#ifndef BEAVER_DILOG_HPP
#define BEAVER_DILOG_HPP
#include "log.hpp"
#include "config.hpp"
#include <cmath>

namespace beaver {
  namespace internals::dilog{
    //universal constants
  const double ln2=0.69314718055994530942;
  const double zeta2=1.6449340668482264365;
  //MiniMax numerator coefficients of log(1-x)/x
  const double P1[]={-0.9999999999999998672,2.7601168443679114620,-2.7923621883637091036,1.2497503713725827027,-0.23173044836039025061,0.012188323559602391476};
    //MiniMax denominator coefficients of log(1-x)/x
  const double Q1[]={1.0000000000000000000,-3.2601168443679653008,4.0890872772181186328,-2.4575883952960021967,0.71252476615980667079,-0.086169703890931582878,0.0026519629643638051268};
    //MiniMax numerator coefficients of Li3(x)/x on 10^(-4)<x<0.5
  const double P2a[]={0.9999999999999999482,-2.6884303049105920491,2.6478043064679810600,-1.1539162530148495100,0.20887751439651881693,-0.010861077434674471020};        
  //MiniMax denominator coefficients of Li3(x)/x on 10^(-4)<x<0.5
  const double Q2a[]={1.0000000000000000000,-2.9384303049106135224,3.2713007715860567301,-1.7077491898534406365,0.41598884303385831214,-0.039805247256814134926,0.00082755014107320443509}; 
    //MiniMax numerator coefficients of Li2(x)/x on -1<x<10^(-4)
  const double P2b[]={0.9999999999999999277,-2.1758446306576192119,1.6555942695751114745,-0.52288981465311876594,0.062597549060802191130,-0.0018723824826943140006};        
  //MiniMax denominator coefficients of Li2(x)/x on -1<x<10^(-4)
  const double Q2b[]={1.0000000000000000000,-2.4258446306575919627,2.1509443161301692038,-0.85358760134515718081,0.14861592651355131752,-0.0093613047200351782876,0.00011533394212474942253}; 
     //Li2 kernel around x=0
  }
  static double Li2_Taylor(double x){
        double x2=x*x;
        double x4=x2*x2;
        return x+(1.0/4+1.0/9*x)*x2+1.0/16*x4;
    }
    BEAVER_NODISCARD inline double  dilog(double x) noexcept {
    namespace  LOC=internals::dilog;
    double taylorswitch=1e-4;//switch up to which an expansion about x=0 or x=1 is used. Guarantees precission for small y.
    double invtaylorswitch=1/taylorswitch;
  //Catch non-finite input
  if (!std::isfinite(x)) return std::numeric_limits<double>::quiet_NaN();
  if(x<-invtaylorswitch){
    double l=beaver::log(-x);
    double c=-LOC::zeta2-1.0/2*l*l;
    double xinv=1.0/x;
    double li2inv=Li2_Taylor(xinv);
    return -li2inv+c;
  }else if(x<-1){
    double l=beaver::log(-x);
    double c=-LOC::zeta2-1.0/2*l*l;
    double x2=x*x;
    double x4=x2*x2;
    double x6=x4*x2;
    double p=LOC::P2b[5]+x*LOC::P2b[4]+x2*(LOC::P2b[3]+x*LOC::P2b[2])+x4*(LOC::P2b[1]+x*LOC::P2b[0]);
    double q=LOC::Q2b[6]+x*LOC::Q2b[5]+x2*(LOC::Q2b[4]+x*LOC::Q2b[3])+x4*(LOC::Q2b[2]+x*LOC::Q2b[1])+x6*LOC::Q2b[0];
    return -p/q+c;
  }else if(x==-1){
    return -0.5*LOC::zeta2;
  }else if(x<-taylorswitch){
    double x2=x*x;
    double x4=x2*x2;
    double x6=x4*x2;
    double p=x*LOC::P2b[0]+x2*(LOC::P2b[1]+x*LOC::P2b[2])+x4*(LOC::P2b[3]+x*LOC::P2b[4])+x6*LOC::P2b[5];
    double q=LOC::Q2b[0]+x*LOC::Q2b[1]+x2*(LOC::Q2b[2]+x*LOC::Q2b[3])+x4*(LOC::Q2b[4]+x*LOC::Q2b[5])+x6*LOC::Q2b[6];
    return p/q;
  }else if(x<taylorswitch){
    return Li2_Taylor(x);
  }else if(x<0.5){
    double x2=x*x;
    double x4=x2*x2;
    double x6=x4*x2;
    double p=x*LOC::P2a[0]+x2*(LOC::P2a[1]+x*LOC::P2a[2])+x4*(LOC::P2a[3]+x*LOC::P2a[4])+x6*LOC::P2a[5];
    double q=LOC::Q2a[0]+x*LOC::Q2a[1]+x2*(LOC::Q2a[2]+x*LOC::Q2a[3])+x4*(LOC::Q2a[4]+x*LOC::Q2a[5])+x6*LOC::Q2a[6];
    return p/q;
  }else if(x==0.5)
  {
    return 0.5*LOC::zeta2-0.5*LOC::ln2*LOC::ln2;
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
    double c=LOC::zeta2-l*ly;
    double p=y*LOC::P2a[0]+y2*(LOC::P2a[1]+y*LOC::P2a[2])+y4*(LOC::P2a[3]+y*LOC::P2a[4])+y6*LOC::P2a[5];
    double q=LOC::Q2a[0]+y*LOC::Q2a[1]+y2*(LOC::Q2a[2]+y*LOC::Q2a[3])+y4*(LOC::Q2a[4]+y*LOC::Q2a[5])+y6*LOC::Q2a[6];
    return -p/q+c;
  } else if(x<1){
    double y=1-x;
    double y2=y*y;
    double y4=y2*y2;
    double ly=beaver::log(y);
    double taylor=LOC::zeta2-y-y2*(0.25+1.0/9*y)-1.0/16*y4;
    double logtaylor=y+y2*(1.0/2+1.0/3*y)+1.0/4*y4;
    return taylor+ly*logtaylor;
  }else if(x==1){
    return LOC::zeta2;
  }else{//Argument on branch-cut
    return std::numeric_limits<double>::quiet_NaN(); 
  }
}

}
#endif