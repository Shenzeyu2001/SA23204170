#' @title A banding parameter estimator for banded autocovariance matrix 
#' @description A banding parameter estimator for banded autocovariance matrix 
#' @param data  observation vector
#' @param K  a parameter that generate K by K  submatrix
#' @param s block size parameter
#' @return  estimated band
#' @export
BACM<-function(data,K,s){
  stopifnot(is.vector(data))
  stopifnot(K>0)
  stopifnot(s>0)
  stopifnot(s>K)
  n <- length(data)
  ###### computer sample autocovariance function######
  gamma_hat<-numeric(K+1)
  gamma_hat[1]<-sum(data^2)
  for(k in 1:K){
    gamma_hat[k+1] <- sum( data[-(1:k)] * data[1:(n-k)] )
  }
  gamma_hat <- round(gamma_hat/n,4)
  Sigma_nhat <- matrix(0,nrow=K+1,ncol=K+1)
  for(rno in 1:(K+1)){
    for(cno in rno :(K+1)){
      Sigma_nhat[rno,cno]<-gamma_hat[cno-rno+1]
    }
  }
  Sigma_nhat<- Sigma_nhat+t( Sigma_nhat)-diag(rep(gamma_hat[1],K+1))
  risk<-numeric(K)
  for(l in 1:K){
    matrix_l1norm<-numeric(n-s+1)
    for(t in 1:(n-s+1)){
      residual_truncated<-data[t:(t+s-1)]
      truncated_gamma_hat<-numeric(l+1)
      truncated_gamma_hat[1]<-sum(residual_truncated^2)
      for(ind in 1 : l ){
        truncated_gamma_hat[ind+1]<-
          sum( residual_truncated[-(1:ind)] * residual_truncated[1:(s-ind)] )
      }
      truncated_gamma_hat <- round(truncated_gamma_hat/s,4)
      Sigma_slt_hat<-matrix(0,nrow=K+1,ncol=K+1)
      for(rno in 1:(K+1)){
        for(cno in rno :(K+1)){
          if(abs(rno-cno)<=l){
            for(band in 1:l){
              if(abs(rno-cno)==0)
                Sigma_slt_hat[rno,cno]=truncated_gamma_hat[1]
              if(abs(rno-cno)==band)
                Sigma_slt_hat[rno,cno]=truncated_gamma_hat[band+1]
            }
          }
        }
      }
      
      Sigma_slt_hat<-Sigma_slt_hat+t(Sigma_slt_hat)-diag(rep(truncated_gamma_hat[1],K+1))
      maxcol<-apply(abs(Sigma_slt_hat-Sigma_nhat),1,sum)
      matrix_l1norm[t]<-maxcol[which.max(maxcol)]
    }
    risk[l]<-round(sum(matrix_l1norm)/(n-s+1),4)
  }
  
  return(band=which.min(risk))
}

#' @title Benchmark R and Rcpp functions.
#' @name benchmarks
#' @description Use R package \code{microbenchmark} to compare the performance of R function (\code{BACM_R}) and Cpp functions (\code{BACM_C}).
#' @import microbenchmark
#' @import Rcpp
#' @import stats 
NULL

#' @title homework
#' @name homework
#' @description All the homework in this semester
#' @import boot
#' @import coda
NULL

#' @title A dataset used for illustration.
#' @name data
#' @description This dataset includes  observations about weight and height, used in homework.
NULL