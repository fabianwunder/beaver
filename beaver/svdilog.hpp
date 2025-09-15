#ifndef BEAVER_SVDILOG_HPP
#define BEAVER_SVDILOG_HPP
#include "log.hpp"
#include "config.hpp"
#include <cmath>

namespace beaver {
    BEAVER_NODISCARD inline double  svdilog(double x) noexcept {

  //Catch non-finite input
  if (!std::isfinite(x)) return std::numeric_limits<double>::quiet_NaN();
  //universal constants
  const double zeta2=1.6449340668482264365;
  const double ln2=0.69314718055994530942;
  //MiniMax numerator coefficients of log(1-x)/x
  const double P1[]={-0.9999999999999998672,
                            2.7601168443679114661,
                          -2.7923621883637091130,
                          1.2497503713725827098,
                          -0.23173044836039025264,
                          0.012188323559602391638};
  //MiniMax denominator coefficients of log(1-x)/x
  const double Q1[]={1.0000000000000000000,
                          -3.2601168443679653050,
                           4.0890872772181186442,
                          -2.4575883952960022081,
                           0.71252476615980667576,
                          -0.086169703890931583743,
                           0.0026519629643638051651};
  //MiniMax numerator coefficients of Li2(x)/x
  const double P2[]={0.9999999999999999482,
                    -2.6884303049105920449,
                    2.6478043064679810510,
                    -1.1539162530148495033,
                    0.20887751439651881509,
                    -0.010861077434674470876};        
  //MiniMax denominator coefficients of Li2(x)/x
  const double Q2[]={1.0000000000000000000,
                     -2.9384303049106135183,
                      3.2713007715860567200,
                     -1.7077491898534406278,
                     0.41598884303385830899,
                     -0.039805247256814134496,
                      0.00082755014107320442259}; 
  //Define arguments of SVPs and constants in the mapping formula. Argument always land in (0,0.5).
  double y=0; 
  double c=0; 
  double s=1; 
    if(x < -1.0){
        y=1/(1-x); //SVP argument in mapping
        c=-zeta2; //Constant in mapping
        s=1; //Signum in mapping
    } else if(x < 0.0){
        y=x/(x-1); //SVP argument in mapping
        c=0; //Constant in mapping
        s=-1; //Signum in mapping
    } else if(x==0.0){
        return 0.0;
    } else if(x <= 0.5){
        y=x; //SVP argument in mapping
        c=0; //Constant in mapping
        s=1; //Signum in mapping
    } else if(x < 1.0){
        y=1-x; //SVP argument in mapping
        c=zeta2; //Constant in mapping
        s=-1; //Signum in mapping
    } else if(x==1.0){
        return zeta2;  
    } else if(x < 2.0){
        y=1-1/x; //SVP argument in mapping
        c=zeta2; //Constant in mapping
        s=1; //Signum in mapping
    } else {
        y=1/x; //SVP argument in mapping
        c=2.0*zeta2; //Constant in mapping
        s=-1; //Signum in mapping
    }
  const double y2=y*y;
  const double y4=y2*y2;
  const double tinyswitch=1e-4;//switch up to which an expansion about y=0 is used. Guarantees precission for small y.
    if(y<tinyswitch) {
        const double l=beaver::log(y);//Fast implementation of log(y) that is stable near 0
        //C= y+1/4*y^2+1/9*y^3+1/16*y^4 (non-log part from Li2)
        const double C=y+std::fma(1.0/9.0,y,0.25)*y2+1.0/16.0*y4;//polynomial part of small y expansion
        // D = -1/2*y - 1/4*y^2 - 1/6*y^3 - 1/8*y^4 (log part from log(1-x))
        const double D = -0.5*y + std::fma(-1.0/6.0, y, -0.25)*y2 - 0.125*y4;
        return c + s*std::fma(l,D,C); // small y expression including mapping constants and signs
    }
    else {
        //log(1-x)/x MiniMax numerator polynomial in Estrin scheme
        const double p1=std::fma(P1[1],y,P1[0])+y2*std::fma(P1[3],y,P1[2])+y4*std::fma(P1[5],y,P1[4]);

        //log(1-x)/x MiniMax denominator polynomial in Estrin scheme
        const double q1_t0=std::fma(Q1[1],y,Q1[0]);
        const double q1_t2=std::fma(Q1[3],y,Q1[2]);
        const double q1_t4=std::fma(Q1[5],y,Q1[4]);
        const double q1=q1_t0 + y2*q1_t2 + y4*std::fma(Q1[6], y2, q1_t4);

        //Li2(x)/x MiniMax numerator polynomial in Estrin scheme
        const double p2=std::fma(P2[1],y,P2[0])+y2*std::fma(P2[3],y,P2[2])+y4*std::fma(P2[5],y,P2[4]);

        //Li2(x)/x MiniMax denominator polynomial in Estrin scheme
        const double q2_t0=std::fma(Q2[1],y,Q2[0]);
        const double q2_t2=std::fma(Q2[3],y,Q2[2]);
        const double q2_t4=std::fma(Q2[5],y,Q2[4]);
        const double q2=q2_t0 + y2*q2_t2 + y4*std::fma(Q2[6], y2, q2_t4);

        //Denominator for log(1-x)/x
        const double denomq1half=1.0/(2.0*q1);
        //Denominator for Li2(x)/x
        const double denomq2=1.0/q2;
        //const double denomql=1.0/ql;
        //Final result for svp(2,x)=Li2(y)+1/2*log(y)*log(1-y) (argument is apped to fall in 0<y<=0.5)
        return std::fma(s*y*p1*denomq1half,beaver::log(y),std::fma(s*denomq2,y*p2,c));
    }
    }

}
#endif