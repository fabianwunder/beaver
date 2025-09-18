#ifndef BEAVER_ATANINT_HPP
#define BEAVER_ATANINT_HPP
#include "log.hpp"
#include "config.hpp"
#include <cmath>


namespace beaver {
    //universal constants
  const double pi=3.1415926535897932385;
  const double catalan=0.91596559417721901505;
    //MiniMax numerator coefficients of Ti2(x)/x on 10^(-4)<x<1
  const double P2[]={0.99999999999999999886,1.0100093542230733512,2.1455220673317650277,1.5305076352812549690,1.4078551331250402304,0.64227878391676890953,0.28789298099160559927,0.063249445901043777286,0.0090534760702843347330};
    //MiniMax denominator coefficients of Ti2(x)/x on 10^(-4)<x<1
  const double Q2[]={1.0000000000000000000,1.0100093542230725212,2.2566331784429801196,1.6427308968564321342,1.6185921530867935222,0.78440406282611283524,0.39787941177428306790,0.10530849411572513926,0.022227486952458724944,0.0013776205701865119100};
         //Ti2 kernel around x=0
  static double Ti2_Taylor(double x){
        double x2=x*x;
        return x-1.0/9*x2*x;
    }
  BEAVER_NODISCARD inline double  atanint(double x) noexcept {
    double taylorswitch=1e-4;//switch up to which an expansion about x=0 or x=1 is used. Guarantees precission for small y.
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
      return sign*Ti2_Taylor(y);
    }else if(y<1){
      double y2=y*y;
      double y4=y2*y2;
      double y6=y4*y2;
      double y8=y4*y4;
      double p=y*P2[0]+y2*(P2[1]+y*P2[2])+y4*(P2[3]+y*P2[4])+y6*(P2[5]+y*P2[6])+y8*(P2[7]+y*P2[8]);
      double q=Q2[0]+y*Q2[1]+y2*(Q2[2]+y*Q2[3])+y4*(Q2[4]+y*Q2[5])+y6*(Q2[6]+y*Q2[7])+y8*(Q2[8]+y*Q2[9]);
      return sign*p/q;
    }else if(y==1.0){
      return sign*catalan;
    }else if(y<invtaylorswitch){
      double l=beaver::log(y);
      double c=sign*0.5*pi*l;
      double y2=y*y;
      double y4=y2*y2;
      double y6=y4*y2;
      double y8=y4*y4;
      double p=P2[8]+y*P2[7]+y2*(P2[6]+y*P2[5])+y4*(P2[4]+y*P2[3])+y6*(P2[2]+y*P2[1])+y8*P2[0];
      double q=Q2[9]+y*Q2[8]+y2*(Q2[7]+y*Q2[6])+y4*(Q2[5]+y*Q2[4])+y6*(Q2[3]+y*Q2[2])+y8*(Q2[1]+y*Q2[0]);
      return sign*p/q+c;
    }else{
      double l=beaver::log(y);
      double c=sign*0.5*pi*l;
      double yinv=1.0/y;
      double ti2inv=Ti2_Taylor(yinv);
      return sign*ti2inv+c;
    }
  }
}
#endif

