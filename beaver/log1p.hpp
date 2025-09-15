#ifndef BEAVER_LOG1P_HPP
#define BEAVER_LOG1P_HPP
#include "log.hpp"
#include "config.hpp"
#include <cmath>


namespace beaver {
  namespace internals::log1p{
    //MiniMax numerator coefficients of arctan(x)/x on -0.3<x<0.3
  const double P1[]={1.0000000000000000088,2.5354992169744272295,2.3128650957787127388,0.90911801705539161145,0.14229189572592744963,0.0059080279683984783396};
      //MiniMax denominator coefficients of arctan(x)/x on -0.3<x<0.3
  const double Q1[]={1.0000000000000000000,3.0354992169744259496,3.4972813709325962709,1.8959256301972141936,0.48336905809034197355,0.049504512928543790754,0.0012137637778754891357};
  }
   //log(1+x) kernel around x=0
  static double log1p_Taylor(double x){
        double x2=x*x;
        double x4=x2*x2;
        return x+x2*(-1.0/2+1.0/3*x)+x4*(-1.0/4+1.0/5*x);
    }
  BEAVER_NODISCARD inline double  log1p(double x) noexcept {
    namespace  LOC=internals::log1p;
    double taylorswitch=1e-3;//switch up to which an expansion about x=0 . Guarantees precission for small x.
    double logswitch=0.3;
    //Catch non-finite input
    if (!std::isfinite(x)) return std::numeric_limits<double>::quiet_NaN();
    if(x==0){
      return 0.0;
    }
    if(x==-1){
      return -INFINITY;
    }
    //on branch-cut
    if(x<-1){
      return NAN;
    }
      // Branchless magnitude and sign
    double   y=std::fabs(x);
    //
    if(y<taylorswitch){
      return log1p_Taylor(x);
    }else if(y<logswitch){
        double x2=x*x;
      double x4=x2*x2;
      double x6=x4*x2;
      double p=x*LOC::P1[0]+x2*(LOC::P1[1]+x*LOC::P1[2])+x4*(LOC::P1[3]+x*LOC::P1[4])+x6*LOC::P1[5];
      double q=LOC::Q1[0]+x*LOC::Q1[1]+x2*(LOC::Q1[2]+x*LOC::Q1[3])+x4*(LOC::Q1[4]+x*LOC::Q1[5])+x6*LOC::Q1[6];
      return p/q;
    }else{
      return beaver::log(1+x);
    }
  }
}
#endif

