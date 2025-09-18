#ifndef BEAVER_ARCTAN_HPP
#define BEAVER_ARCTAN_HPP
//#include "log.hpp"
#include "config.hpp"
#include <cmath>


namespace beaver {
  namespace internals::arctan{
    //universal constants
  const double pi=3.1415926535897932385;
    //MiniMax numerator coefficients of arctan(x)/x on 10^(-4)<x<1
  const double P1[]={1.0000000000000000031,1.8342348086649654905,3.3257594029617757468,3.5990739720281070481,3.2250546833891888322,2.0790695069014332191,1.0011354591290797518,0.31803640663671448932,0.052547928299761760621}; 
    //MiniMax denominator coefficients of arctan(x)/x on 10^(-4)<x<1
  const double Q1[]={1.0000000000000000000,1.8342348086649679506,3.6590927362947704617,4.2104855749348947173,4.2447522616248727101,3.1157177460107374180,1.8270913698071677130,0.77654625981176704113,0.22423950752923047543,0.033419821882336195679};
    //arctan kernel around x=0
  }
  static double arctan_Taylor(double x){
        double x2=x*x;
        return x-1.0/3*x2*x;
    }
  

  BEAVER_NODISCARD inline double  arctan(double x) noexcept {
    namespace  LOC=internals::arctan;
    double taylorswitch=1e-3;//switch up to which an expansion about x=0 or x=1 is used. Guarantees precission for small y.
    double invtaylorswitch=1/taylorswitch;
    //Catch non-finite input
    if (!std::isfinite(x)) return std::numeric_limits<double>::quiet_NaN();
    if(x==0){
      return 0.0;
    }
      // Branchless magnitude and sign
    double   sign = std::copysign(1.0,x);
    double   y=std::fabs(x);
    //
    if(y<taylorswitch){
      return sign*arctan_Taylor(y);
    }else if(y<1){
      double y2=y*y;
      double y4=y2*y2;
      double y6=y4*y2;
      double y8=y4*y4;
      double p=y*LOC::P1[0]+y2*(LOC::P1[1]+y*LOC::P1[2])+y4*(LOC::P1[3]+y*LOC::P1[4])+y6*(LOC::P1[5]+y*LOC::P1[6])+y8*(LOC::P1[7]+y*LOC::P1[8]);
      double q=LOC::Q1[0]+y*LOC::Q1[1]+y2*(LOC::Q1[2]+y*LOC::Q1[3])+y4*(LOC::Q1[4]+y*LOC::Q1[5])+y6*(LOC::Q1[6]+y*LOC::Q1[7])+y8*(LOC::Q1[8]+y*LOC::Q1[9]);
      return sign*p/q;
    }else if(y==1.0){
      return sign*0.25*LOC::pi;
    }else if(y<invtaylorswitch){
      double c=sign*0.5*LOC::pi;
      double y2=y*y;
      double y4=y2*y2;
      double y6=y4*y2;
      double y8=y4*y4;
      double p=LOC::P1[8]+y*LOC::P1[7]+y2*(LOC::P1[6]+y*LOC::P1[5])+y4*(LOC::P1[4]+y*LOC::P1[3])+y6*(LOC::P1[2]+y*LOC::P1[1])+y8*LOC::P1[0];
      double q=LOC::Q1[9]+y*LOC::Q1[8]+y2*(LOC::Q1[7]+y*LOC::Q1[6])+y4*(LOC::Q1[5]+y*LOC::Q1[4])+y6*(LOC::Q1[3]+y*LOC::Q1[2])+y8*(LOC::Q1[1]+y*LOC::Q1[0]);
      return -sign*p/q+c;
    }else{
      double c=sign*0.5*LOC::pi;
      double yinv=1.0/y;
      double arctaninv=arctan_Taylor(yinv);
      return -sign*arctaninv+c;
    }
  }
}
#endif

