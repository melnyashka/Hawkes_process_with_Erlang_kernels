
#Exponential matrix: Splitting 
expmatA_SP<-function(t,ny1,ny2,eta1,eta2)
{
  eAny1<-matrix(0,nrow=eta1+1,ncol=eta1+1)
  eat<-exp(-ny1*t)
  for (i in 1:(eta1+1)){
    for (j in 1:(eta1+1-i+1)){
      k<-j-1
      l<-j+i-1
      eAny1[i,l]<-eat*(t^k/factorial(k))
    }
  }
  Null1<-matrix(0,nrow=eta1+1,ncol=eta2+1)
  Mat1<-cbind(eAny1,Null1)
  
  eAny2<-matrix(0,nrow=eta2+1,ncol=eta2+1)
  ebt<-exp(-ny2*t)
  for (i in 1:(eta2+1)){
    for (j in 1:(eta2+1-i+1)){
      k<-j-1
      l<-j+i-1
      eAny2[i,l]<-ebt*(t^k/factorial(k))
    }
  }
  Null2<-matrix(0,nrow=eta2+1,ncol=eta1+1)
  Mat2<-cbind(Null2,eAny2)
  ret<-rbind(Mat1,Mat2)
  
  return (ret)
}


