#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

//#f1 function 
// [[Rcpp::export]]
double f1_EM_Cpp_(double x)
{
  double ret=0;
  if(x<log(20)){
    ret=10*exp(x);
  }
  if(x>=log(20)){
    ret=400/(1+400*exp(-2*x));
  }
  return ret;
};

//#f2 function 
// [[Rcpp::export]]
double f2_EM_Cpp_(double x)
{
  double ret=0;
  if(x<log(20)){
    ret=exp(x);
  }
  if(x>=log(20)){
    ret=40/(1+400*exp(-2*x));
  }
  return ret;
};

//Drift component
// [[Rcpp::export]]
NumericVector drift_EM_Cpp_(NumericVector X, double c1, double c2, double ny1, double ny2, double eta1, double eta2)
{
  double kappa=eta1+eta2+2;
  NumericVector ret(kappa);
  
  for(int i=0;i<(eta1);i++)
  {
    ret(i)=-ny1*X(i)+X(i+1);
  }
  ret(eta1)=-ny1*X(eta1)+c1*f2_EM_Cpp_(X(eta1+1));
  for(int i=(eta1+1);i<(kappa-1);i++)
  {
    ret(i)=-ny2*X(i)+X(i+1);
  }
  ret(kappa-1)=-ny2*X(kappa-1)+c2*f1_EM_Cpp_(X(0));
  return ret;
};

//Diffusion component
// [[Rcpp::export]]
NumericVector diffusion_EM_Cpp_(NumericVector X, NumericVector randvec, double c1, double c2, double p1, double p2, double eta1, double eta2)
{
  double kappa=eta1+eta2+2;
  NumericVector ret(kappa);
  ret(eta1)=(c1/sqrt(p2))*sqrt(f2_EM_Cpp_(X(eta1+1)))*randvec(0); 
  ret(kappa-1)=(c2/sqrt(p1))*sqrt(f1_EM_Cpp_(X(0)))*randvec(1);
  return ret;
};

//Euler-Maruyama 
// [[Rcpp::export]]
NumericMatrix Hawkes_EM_Cpp_(NumericVector grid_i, double h_i, NumericVector start_i, double c1_i, double c2_i, double ny1_i, double ny2_i, double p1_i, double p2_i, double N_i, double eta1_i, double eta2_i)
{
  double h=h_i;
  NumericVector start=start_i;
  NumericVector grid=grid_i;
  int iter=grid.size();
  
  NumericMatrix randarr(2,iter);
  randarr(0,_)=sqrt(h)*rnorm(iter);
  randarr(1,_)=sqrt(h)*rnorm(iter);

  //parameter values
  double c1=c1_i;
  double c2=c2_i;
  double ny1=ny1_i;
  double ny2=ny2_i;
  double p1=p1_i;
  double p2=p2_i;
  double N=N_i;
  double eta1=eta1_i;
  double eta2=eta2_i;
  
  double kappa=eta1+eta2+2;

  NumericMatrix sol(kappa,iter);
  sol(_, 0)=start;
  NumericVector newv=start;
  NumericVector randvec;
  
  for(int i=1;i<iter;i++)
  {
    randvec=randarr(_,i);
    newv=newv+h*drift_EM_Cpp_(newv,c1,c2,ny1,ny2,eta1,eta2)+(1/sqrt(N))*diffusion_EM_Cpp_(newv,randvec,c1,c2,p1,p2,eta1,eta2);
    sol(_,i)=newv;
  }
  
  NumericMatrix ret=sol;
  
  return sol;
};


