/*****************************************************************
copyright GPL2. Please send comments, bug 
          report etc. to mai@ms.uky.edu  or
 
            Mai Zhou
            Department of Statistics
            University of Kentucky
            Lexington, KY 40506-0027
*******************************************************************/


#include "utils.h"

void wur(int * d, double * wt, double * sur, double * jum, int n){  
   int i;
   double a=0.0;

   for(i=0; i<n; i++)
       a+=wt[i];        
/* this computes the total weight a, which may not=1 */
/* Next compute the first sur and jum */
   if(d[0]==0){
       sur[0]=1.0;
       jum[0]=0.0;
   }
   else {
       sur[0]=1.0-wt[0]/a;
       jum[0]=1.0-sur[0];
   }
/* Now compute the rest of sur and jum in a loop */ 
   for(i=1; i<n; i++){
       a-=wt[i-1];
       if(d[i]==0){
            sur[i]=sur[i-1];
            jum[i]=0.0;
       }
       else {
            sur[i]=sur[i-1]*(1.0-wt[i]/a);
            jum[i]=sur[i-1]-sur[i];
       }
   }
}

