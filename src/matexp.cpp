#include <iostream>
//#include <armadillo>
#include <cmath>
#include "RcppArmadillo.h"
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <chrono>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

using namespace std;

/*
 void printdvec(double *pv, int n) {
 int i;
 for (i=0;i<n;i++) {
 cout << pv[i]<<", ";
 }
 cout <<"\n";
 }
 void printivec(int *pv, int n) {
 int i;
 for (i=0;i<n;i++) {
 cout << pv[i]<<", ";
 }
 cout <<"\n";
 }
 */

inline double h(double x) {
  return 1.0-x+x*log(x);
}

inline double hifunc(double rho, double B,double loge) {
  return rho +(B-loge)/3.0*(1.0+sqrt(1.0+18.0*rho/(B-loge)));
}
inline double lofunc(double rho, double A, double loge) {
  const double logroot2pi=0.5*log(2*3.14159265);
  return rho+sqrt(2*rho)*sqrt(-logroot2pi-loge-1.5*log(A)+log(A-1));
}
int get_mlo(int mhi, double rho) {
  int mlo;
  mlo=2*(int)(rho-0.5) - mhi;
  if (mlo<0) {
    mlo=0;
  }
  return mlo;
}
// Take the L1 norm of the matrix and return the mmax for the req precision
// [[Rcpp::export]]
int get_m(double rho, double prec, int mlo=-1) {
  const double logprec=log(prec), pi=3.14159265;
  double dmhi, dmlo;
  int mhi;
  
  dmhi= hifunc(rho,0.0,logprec)-1;
  dmlo=lofunc(rho,2*rho*h(dmhi/rho),logprec);
  if ((int)dmlo > mlo) {
    mlo=(int)dmlo;
  }
  
  if (log(boost::math::gamma_p((double)(mlo+1),rho))<logprec) {
    return mlo;
  }
  else {
    const double B=-0.5*log(4*pi*rho*h(dmlo/rho));
    if (B>logprec) {
      dmhi=hifunc(rho,B,logprec);
    }
    mhi=(int)(dmhi+1);
    
    //    cout<<mhi<<"-"<<mlo<<"=width: "<< mhi-mlo<<"\n";
    
    while (mhi-mlo>1) {
      int mmid=(mhi+mlo)/2; // rounds down
      double dm=(double)mmid;
      double loginv;
      
      loginv=log(boost::math::gamma_p(dm+1,rho));
      //    cout <<mlo<<", "<<mmid <<", "<<mhi<<", "<<loginv<<"\n";
      
      if (loginv<logprec) {
        mhi=mmid;
      }
      else {
        mlo=mmid;
      }
    }
  }
  return mhi;
}

// [[Rcpp::export]]
arma::mat SPS_v_exp_spQ(const arma::mat v, const arma::sp_mat Q, double prec, bool renorm=true, bool t2=true) {
  const double rho=-Q.min(), BIG=1e100;
  const int d=v.n_cols, nc=v.n_rows, mhi=get_m(rho,prec/(1.0+(double)t2)), mlo=t2? get_mlo(mhi,rho): 0;
  const double dnc=(double)nc;
  int j;
  const arma::sp_mat P=Q+rho*speye(size(Q));
  arma::mat vsum(nc,d),vpro=v;
  double fac=1.0, final_mult, szvpro, rncumu=0;
  
  //cout << "rho="<<rho<<", mlo="<< mlo<< " ,mhi="<<mhi << ", d="<<d<<"\n";
  
  szvpro=arma::accu(v)/dnc;
  final_mult=szvpro;
  if (szvpro>BIG) {
    vpro/=szvpro;
    szvpro=1; rncumu=log(szvpro);
  }
  if (mlo==0) {
    vsum=vpro;
  }
  else {
    vsum.zeros();
  }
  
  j=1;
  while (j<=mhi) { // do vsum <- v^T e^P
    vpro=vpro*P;
    vpro/=fac;
    if (j>=mlo) {
      vsum+=vpro; 
    }
    szvpro*=rho/fac++; 
    if (szvpro>BIG) {
      vpro/=szvpro; 
      if (j>=mlo) { 
        vsum/=szvpro; 
      }
      if (!renorm) {
        rncumu+=log(szvpro);
      }
      szvpro=1;
    }
    j++;
  }
  
  // if (pm!=NULL) {
  //   *pm=mhi;
  // }
  if (renorm) {
    final_mult/=(arma::accu(vsum)/dnc);
  }
  else {
    final_mult=exp(rncumu-rho);
  }
  return vsum*final_mult;
}

// Matrix action on a single vector but multiple time points supplied
// in a vector with t_j>=t_i for j>i.
// #if RN == 0 // MUSPS2 i.e. two-tailed but no renormalisation
// [[Rcpp::export]]
arma::mat MUSPS2_v_exp_spQ(arma::rowvec v, arma::sp_mat Q, double prec, arma::vec ts) {
  const double rho=-Q.min(), BIG=1e100, SMALL=1e-100;
  const int d=v.n_cols, nt=ts.n_elem;
  double tmax=ts(nt-1);
  const arma::sp_mat P=Q+rho*speye(size(Q));
  int j,t,tlo,thi;
  double fac=1.0, mmax;
  arma::mat vsum(d,nt); // output
  arma::rowvec vpro;
  int *pmshi=new int[nt], *pmslo=new int[nt+1];
  double szvpro;
  double *prncu=new double[nt], *ptmultfac=new double[nt], *pts=ts.memptr();
  //  cout << "rho="<<rho<<", m="<< m<<"\n";
  
  // Since the ts are ordered, ms(j-1) is a lower bound for ms(j)
  pmshi[0]=get_m(rho*pts[0],prec/2.0);
  pmslo[0]=get_mlo(pmshi[0],rho*pts[0]);
  for (t=1;t<nt;t++) {
    pmshi[t]=get_m(rho*pts[t],prec/2.0,pmshi[t-1]);
    pmslo[t]=get_mlo(pmshi[t],rho*pts[t]);
  }
  // NB: sequence of pmslo need not be monotonic, but this doesn't matter
  mmax=pmshi[nt-1];
  pmslo[nt]=mmax+1; // makes code logic easier later
  //  printivec(pmshi,nt);
  //  printivec(pmslo,nt+1);
  
  szvpro=arma::sum(v);
  if (szvpro>BIG) {
    const double lsvp=log(szvpro);
    vpro=v/szvpro;
    fill(prncu,prncu+nt,lsvp);
    szvpro=1.0;
  }
  else {
    vpro=v;
    fill(prncu,prncu+nt,0.0);
  }
  fill(ptmultfac,ptmultfac+nt,1.0);
  
  vsum.zeros();
  t=0;
  while ((t<nt) && (pmslo[t]==0)) {
    vsum.col(t)=vpro.t();
    t++;
  }
  thi=t-1; tlo=0;
  
  j=1;
  while (j<=mmax) { // do (e^P) vsum
    while (pmslo[thi+1]<=j) {  // add any new times that have come into range
      thi++;          // need <= rather than == because pmslo not nec. monotone
    }
    while (pmshi[tlo]<j) { // remove any that have gone out of range
      tlo++;
    }
    //    cout << "j, tlo and thi="<<j<<", "<<tlo << ", "<<thi<<"\n";
    vpro*=P;
    vpro*=(tmax/fac); // tmax gives simpler (but less efficient) implementation
    szvpro*=(rho*tmax/fac);  // as it leads to more under/overflows
    
    if ((szvpro>BIG) || (ptmultfac[tlo]<SMALL)) {
      const double lsvp=log(szvpro);
      for (t=tlo;t<=thi;t++) {
        vsum.col(t)/=(ptmultfac[t]*szvpro);
      }
      for (t=tlo;t<nt;t++) {
        prncu[t]=prncu[t]+lsvp+log(ptmultfac[t]);
      }
      fill(ptmultfac+tlo,ptmultfac+nt,1.0);
      vpro/=szvpro;
      szvpro=1.0;
    }
    
    for (t=tlo;t<nt;t++) { // keep track for all that could be relevant
      ptmultfac[t]=ptmultfac[t]*pts[t]/tmax;
    }
    for (t=tlo;t<=thi;t++) {
      vsum.col(t)+=ptmultfac[t]*vpro.t();
    }
    fac++; j++;
  }
  
  for (t=0;t<nt;t++) {
    vsum.col(t)*=exp(prncu[t]-rho*pts[t]);
  }
  
  delete[] pmslo; delete[] pmshi; delete[] prncu; delete[] ptmultfac;
  
  
  return vsum.t();
  
}
// #else // MUSPS2r
// [[Rcpp::export]]
arma::mat MUSPS2r_v_exp_spQ(arma::rowvec v, arma::sp_mat Q, double prec, arma::vec ts) {
  const double rho=-Q.min(), BIG=1e100, SMALL=1e-100;
  const int d=v.n_cols, nt=ts.n_elem;
  const arma::sp_mat P=Q+rho*speye(size(Q));
  int j,t,tlo,thi;
  double fac=1.0, mmax;
  arma::mat vsum(d,nt); // output
  arma::rowvec vpro, vcolsums;
  int *pmshi=new int[nt], *pmslo=new int[nt+1];
  double szvpro, final_mult;
  double *ptmultfac=new double[nt], *pts=ts.memptr();
  //  cout << "rho="<<rho<<", m="<< m<<"\n";
  
  // Since the ts are ordered, ms(j-1) is a lower bound for ms(j)
  pmshi[0]=get_m(rho*pts[0],prec/2.0);
  pmslo[0]=get_mlo(pmshi[0],rho*pts[0]);
  for (t=1;t<nt;t++) {
    pmshi[t]=get_m(rho*pts[t],prec/2.0,pmshi[t-1]);
    pmslo[t]=get_mlo(pmshi[t],rho*pts[t]);
  }
  //printivec(pmshi,nt);
  //printivec(pmslo,nt);
  mmax=pmshi[nt-1];
  pmslo[nt]=mmax+1; // makes code logic easier later
  
  szvpro=arma::sum(v);
  final_mult=szvpro;
  if (szvpro>BIG) {
    vpro=v/szvpro;
    szvpro=1.0;
  }
  else {
    vpro=v;
  }
  fill(ptmultfac,ptmultfac+nt,1.0);
  
  vsum.zeros();
  t=0;
  while (pmslo[t]==0) {
    vsum.col(t)=vpro.t();
    t++;
  }
  thi=t-1; tlo=0;
  
  j=1;
  while (j<=mmax) { // do (e^P) vsum
    while (pmslo[thi+1]<=j) {  // add any new times that have come into range
      thi++;
    }
    while (pmshi[tlo]<j) { // remove any that have gone out of range
      tlo++;
    }
    //    cout << "j, tlo and thi="<<j<<", "<<tlo << ", "<<thi<<"\n";
    
    vpro*=P;
    if (thi>-1) { // slightly more efficient (but slightly more complicated)
      const double mul=pts[thi]/fac;
      vpro*=mul;
      szvpro*=(rho*mul);
    }
    else {
      szvpro*=rho;
    }
    
    if ((szvpro>BIG) || (ptmultfac[tlo]<SMALL)) {
      for (t=tlo;t<=thi;t++) {
        vsum.col(t)/=(ptmultfac[t]*szvpro);
      }
      fill(ptmultfac+tlo,ptmultfac+thi+1,1.0);
      vpro/=szvpro;
      szvpro=1.0;
    }
    
    //    for (t=tlo;t<=thi;t++) {
    // }
    for (t=tlo;t<=thi;t++) {
      ptmultfac[t]=ptmultfac[t]*pts[t]/pts[thi];
      vsum.col(t)+=ptmultfac[t]*vpro.t();
    }
    fac++; j++;
  }
  
  vcolsums=arma::sum(vsum);
  for (t=0;t<nt;t++) {
    const double this_mult=final_mult/vcolsums(t);
    vsum.col(t)*=this_mult;
  }
  
  delete[] pmslo; delete[] pmshi; delete[] ptmultfac;
  
  
  return vsum.t();
  
}
// #endif

// arma::mat SPS_exp_spQ_v(arma::sp_mat Q, arma::mat v, double prec, bool ren=true, bool tt=true, int *pm=NULL) {
//   //  arma::sp_mat Qt=Q.t();
//   //arma::mat vt=v.t();
//   return SPS_v_exp_spQ(v.t(), Q.t(), prec, ren, tt, pm).t();
// }
// arma::mat MUSPS_exp_spQ_v(arma::sp_mat Q, arma::mat v, double prec, arma::vec ts) {
//   return MUSPS_v_exp_spQ(v.t(), Q.t(), prec, ts).t();
// }

/*
 int main(int argc, const char** pargv) {
 arma::vec xs;
 arma::vec LVthetas={.3,.4,.01};
 arma::sp_mat Q;
 arma::mat v;
 int i;
 int tst=0;
 
 if (argc>1) {
 tst=atoi(pargv[1]);
 }
 if (tst==-2) {
 double rho=262.0, eps=1e-5;
 int m1,m2;
 m1=get_m(rho,eps);
 m2=get_m(rho,eps,86);
 cout << m1 <<", "<<m2<<"\n";
 }
 if (tst==-1) { // understand creating matrices from pointers in armadillo
 double *pmy=new double[4];
 pmy[0]=0; pmy[1]=1; pmy[2]=2; pmy[3]=3;
 arma::mat my(pmy,2,2);
 cout <<my<<"\n";
 cout << arma::sum(my)<<"\n";
 delete[] pmy;
 }
 
 if (tst==0) { // check get_m against the system's qpois - values and timing
 double rho=0.1;
 double eps=1e-16;
 auto t1 = std::chrono::high_resolution_clock::now();
 double q1;
 int q3;
 int j;
 
 for (j=0;j<8;j++) {
 cout <<"rho="<<rho<<"\n";
 
 for (i=0;i<10000;i++) {
 boost::math::poisson dist(rho);
 q1=quantile(complement(dist,eps));
 }
 auto t2 = std::chrono::high_resolution_clock::now();
 cout<< "Timing for boost: "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() <<"\n";
 
 for (i=0;i<10000;i++) {
 q3=get_m(rho,eps);
 }
 t1 = std::chrono::high_resolution_clock::now();
 cout<< "Timing for get_m: "<<std::chrono::duration_cast<std::chrono::milliseconds>(t1-t2).count() <<"\n";
 cout << q1 << ", "<< q3<<"\n";
 rho=rho*10.0;
 }
 }
 
 if (tst==1) { // test matrix exponential is producing the right answer
 arma::mat LoHi={{1,3},{3,4}};
 arma::mat v(7,1), vout;
 int i;
 double rho;
 for (i=0;i<7;i++) {
 v(i,0)=i+1;
 }
 cout << v.t()<< "\n";
 Q=LoHi_create_Q(LoHi,LVthetas,LV_S,&LV_Rates);
 cout << arma::mat(Q) << "\n";
 rho=-Q.min();
 cout << "rho="<< rho << ", d=" << get_m(rho,1e-16) << "\n";
 cout << SPS_exp_spQ_v(Q.t(),v,1e-16).t() << "\n";
 cout << (arma::expmat(arma::mat(Q.t())) * v).t() << "\n";
 cout << SPS_v_exp_spQ(v.t(),Q,1e-16) << "\n";
 }
 
 if (tst==2) { // check code for multiple time points
 arma::mat LoHi={{1,3},{3,4}};
 arma::mat v(7,1), vout;
 int i;
 double eps=1e-9;
 for (i=0;i<7;i++) {
 v(i,0)=i+1;
 }
 v=v/28.0;
 Q=LoHi_create_Q(LoHi,LVthetas,LV_S,&LV_Rates);
 cout << SPS_v_exp_spQ(v.t(),Q/10.0,eps) << "\n";
 cout << SPS_v_exp_spQ(v.t(),Q/5.0,eps) << "\n";
 cout << SPS_v_exp_spQ(v.t(),Q/2,eps) << "\n";
 cout << SPS_v_exp_spQ(v.t(),Q*2.0,eps) << "\n";
 cout << SPS_v_exp_spQ(v.t(),Q*2.1,eps) << "\n";
 cout << SPS_v_exp_spQ(v.t(),Q*2.2,eps) << "\n";
 cout << SPS_v_exp_spQ(v.t(),Q*20,eps) << "\n";
 cout << SPS_v_exp_spQ(v.t(),Q*100.0,eps,true,true) << "\n";
 arma::vec ts={0.1,0.2,0.5,2.0,2.1,2.2,20,100.0};
 cout << "Times\n"<< ts.t() << "\n";
 cout << "MUSPS\n"<<MUSPS_v_exp_spQ(v.t(),Q,eps,ts) << "\n";
 }
 
 if (tst==3) { // timing with a big Q
 int npop=60, nts=100;
 arma::vec SEIRthetas={.01,.75,.25}; // R0=4
 arma::sp_mat Q=SEIR_create_Qfull(npop,SEIRthetas,SEIR_S,&SEIR_Rates)/100.0;
 int d=Q.n_cols;
 arma::mat v=arma::zeros(d,1);
 arma::mat v1,v2;
 arma::mat storage1=arma::mat(d,nts),storage2=arma::mat(d,nts);
 arma::mat storage3=arma::mat(d,nts);
 arma::vec ts=arma::vec(nts);
 v(d-2)=1.0;
 double eps=1e-14;
 
 cout << "d="<<Q.n_rows <<", non-zero="<<Q.n_nonzero <<"\n";
 
 auto t1 = std::chrono::high_resolution_clock::now();
 for (i=0;i<nts;i++) {
 //      cout << i<<"\n";
 storage1.col(i)=SPS_v_exp_spQ(v.t(),Q*(double)(i+1),eps/100).t();
 }
 auto t2 = std::chrono::high_resolution_clock::now();
 cout<< "Timing for naive: "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() <<"\n";
 
 for (i=0;i<nts;i++) {
 ts(i)=(double)(i+1);
 }
 storage2=MUSPS_v_exp_spQ(v.t(),Q,eps,ts).t();
 t1 = std::chrono::high_resolution_clock::now();
 cout<< "Timing for SPS: "<<std::chrono::duration_cast<std::chrono::milliseconds>(t1-t2).count() <<"\n";
 
 cout<<"Discrepancy: "<<arma::sum(arma::sum(abs(storage2-storage1)))<< "\n";
 
 t1 = std::chrono::high_resolution_clock::now();
 storage3.col(0)=SPS_v_exp_spQ(v.t(),Q,eps).t();
 for (i=1;i<nts;i++) {
 storage3.col(i)=SPS_v_exp_spQ((storage3.col(i-1)).t(),Q,eps).t();
 }
 t2 = std::chrono::high_resolution_clock::now();
 cout<< "Timing for iterative: "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() <<"\n";
 cout<<"Discrepancy: "<<arma::sum(arma::sum(abs(storage3-storage1)))<< "\n";
 
 }
 
 if (tst==4) { // more thorough test of get_m
 int nrhos=7, neps=3;
 double prhos[]={.1,1,10,100,1000,10000,100000};
 double peps[]={1e-5,1e-9,1e-16};
 int nreps=1000,i,j,k;
 arma::mat pprob=arma::mat(neps,nrhos);
 arma::mat pms_CS=arma::mat(neps,nrhos);
 arma::mat pms_boost=arma::mat(neps,nrhos);
 
 if (argc>2) {
 nreps=atoi(pargv[2]);
 cout << "Nreps= "<< nreps << "\n";
 }
 
 for (i=0;i<neps;i++) {
 double eps=peps[i];
 for (j=0;j<nrhos;j++) {
 double rho=prhos[j];
 int m;
 auto t1 = std::chrono::high_resolution_clock::now();
 for (k=0;k<nreps;k++) {
 m=get_m(rho,eps);
 }
 auto t2 = std::chrono::high_resolution_clock::now();
 auto tdiffgetm = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
 pms_CS(i,j)=m-rho;
 
 t1 = std::chrono::high_resolution_clock::now();
 for (k=0;k<nreps;k++) {
 boost::math::poisson dist(rho);
 m=quantile(complement(dist,eps));
 }
 t2 = std::chrono::high_resolution_clock::now();
 auto tdiffboost=std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
 cout<< "Time with (rho,eps)="<<rho<<","<<eps<<"= [getm, boost] " << tdiffgetm<<", "<<tdiffboost<<"\n";
 pms_boost(i,j)=m-rho;
 }
 }
 cout << "ms_CS:\n"<< pms_CS <<"\n";
 //cout << "ms_boost:\n"<< pms_boost <<"\n";
 cout << "diff in ms:\n"<<pms_CS-pms_boost<<"\n";
 // so accuracy is same - NB R's qpois() gets the wrong answer for high rho and low epsilon
 }
 
 return 0;
 }
 */
