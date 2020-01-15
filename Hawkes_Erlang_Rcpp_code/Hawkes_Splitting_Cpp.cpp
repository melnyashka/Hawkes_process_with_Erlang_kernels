#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

//#Matrix-Vector Multiplication mv_mult
// [[Rcpp::export]]
NumericVector mv_mult_Hawkes_(NumericMatrix mat, NumericVector vec)
{
  NumericVector ret(mat.nrow());
  double temp=0;
  
  for(int i = 0; i < mat.nrow(); i++)
  {
    for(int j = 0; j < vec.size(); j++)
    {
      temp = temp + mat(i,j) * vec[j];
    }
    ret[i]=temp;
    temp=0;
  }
  return ret;
};

// #f1 function 
// [[Rcpp::export]]
double f1_SP_Cpp_(double x)
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
double f2_SP_Cpp_(double x)
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

//ODE subsystem
// [[Rcpp::export]]
NumericVector ODE_SP_Cpp_(NumericVector X, NumericMatrix expmatA, double ny1, double ny2)
{
  NumericVector ret=mv_mult_Hawkes_(expmatA,X);
  return ret;
};

//SDE subsystem drift
// [[Rcpp::export]]
NumericVector SDE_drift_SP_Cpp_(NumericVector X, NumericVector randvec, double c1, double c2, double eta1, double eta2)
{
  double kappa = eta1+eta2+2;
  NumericVector ret(kappa);
  ret(eta1)=c1*f2_SP_Cpp_(X(eta1+1));
  ret(kappa-1)=c2*f1_SP_Cpp_(X(0));
  return ret;
};

//SDE subsystem diffusion
// [[Rcpp::export]]
NumericVector SDE_diffusion_SP_Cpp_(NumericVector X, NumericVector randvec, double c1, double c2, double p1, double p2, double eta1, double eta2)
{
  double kappa=eta1+eta2+2;
  NumericVector ret(kappa);
  ret(eta1)=(c1/sqrt(p2))*sqrt(f2_SP_Cpp_(X(eta1+1)))*randvec(0);
  ret(kappa-1)=(c2/sqrt(p1))*sqrt(f1_SP_Cpp_(X(0)))*randvec(1);
  return ret;
};

//Splitting: Strang
// [[Rcpp::export]]
NumericMatrix Hawkes_SP_ST_Cpp_(NumericVector grid_i, double h_i, NumericVector start_i, double c1_i, double c2_i, double ny1_i, double ny2_i, double p1_i, double p2_i, double N_i, NumericMatrix expmatA_i, double eta1_i, double eta2_i)
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
  
  //Matrix Exponential
  NumericMatrix expmatA=expmatA_i;

  NumericMatrix sol(kappa,iter);
  sol(_, 0)=start;
  NumericVector newv=start;
  NumericVector randvec;
  
  for(int i=1;i<iter;i++)
  {
    randvec=randarr(_,i);
    newv=ODE_SP_Cpp_(newv,expmatA,ny1,ny2);
    newv=newv+h*SDE_drift_SP_Cpp_(newv,randvec,c1,c2,eta1,eta2)+(1/sqrt(N))*SDE_diffusion_SP_Cpp_(newv,randvec,c1,c2,p1,p2,eta1,eta2);
    newv=ODE_SP_Cpp_(newv,expmatA,ny1,ny2);
    sol(_,i)=newv;
  }
  
  NumericMatrix ret=sol;
  
  return sol;
};

//Splitting: Lie-Trotter
// [[Rcpp::export]]
NumericMatrix Hawkes_SP_LT_Cpp_(NumericVector grid_i, double h_i, NumericVector start_i, double c1_i, double c2_i, double ny1_i, double ny2_i, double p1_i, double p2_i, double N_i, NumericMatrix expmatA_i, double eta1_i, double eta2_i)
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
  
  //Matrix Exponential
  NumericMatrix expmatA=expmatA_i;
  
  NumericMatrix sol(kappa,iter);
  sol(_, 0)=start;
  NumericVector newv=start;
  NumericVector randvec;
  
  for(int i=1;i<iter;i++)
  {
    randvec=randarr(_,i);
    newv=newv+h*SDE_drift_SP_Cpp_(newv,randvec,c1,c2,eta1,eta2)+(1/sqrt(N))*SDE_diffusion_SP_Cpp_(newv,randvec,c1,c2,p1,p2,eta1,eta2);
    newv=ODE_SP_Cpp_(newv,expmatA,ny1,ny2);
    sol(_,i)=newv;
  }
  
  NumericMatrix ret=sol;
  
  return sol;
};

