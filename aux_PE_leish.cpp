#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

// this function calculates the invasion pressure variable
// [[Rcpp::export]]
NumericMatrix CalcInvasionPressure(IntegerMatrix z, NumericMatrix OneOverDist, 
                                   int nanos, int nlocs,
                                   NumericVector SomaOneOverDist) {
  
  NumericMatrix IP(nlocs,nanos);
  double res;
  for(int i=0; i<nlocs;i++){
    for (int j=0; j<(nanos-1); j++){
      res=0;
      for (int k=0; k<nlocs; k++){
        res=res+z(k,j)*OneOverDist(i,k);
      }
      IP(i,j)=res/SomaOneOverDist[i];
    }
  }
  return IP;
}

// this function identifies the leading 1's
// notice that this function only works properly if z is such that, once invaded, always invaded
// [[Rcpp::export]]
IntegerVector IdentifLead1(IntegerMatrix z, int nanos, int nlocs,
                           IntegerMatrix Identifiers) {

 IntegerVector id(nlocs);
 for(int i=0; i<nlocs;i++){
   for (int j=1; j<nanos; j++){
     if ((z(i,j)==1) & (z(i,j-1)==0)){
       id[i]=Identifiers(i,j);
     }
   }
 }
 return id;
}

