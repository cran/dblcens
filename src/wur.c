/*****************************************************************
copyright GPL2. Please send comments, bug 
          report etc. to mai@ms.uky.edu  or
 
            Mai Zhou
            Department of Statistics
            University of Kentucky
            Lexington, KY 40506-0027
*******************************************************************/
/* The following function wur is to compute weighted Kaplan-Meier.      */
/* Although z did not appear, the data, d, wt, etc. are                 */ 
/* sorted according to z value.                                         */ 
/*                                                                      */ 
/* It computes the weighted Kaplan-Meier estimator, did not             */ 
/* remove the censored observation (i.e. leave the observations         */ 
/* as is. which is z, d, and wt. where wt do not have to sum to one.    */ 
/* The censoring indicator, d[i], do not have to be 0/1 value,          */ 
/* in fact, the only thing matters is wheather d[i]=0 or not.           */ 
/*                                                                      */ 
/* The results are in sur, jum. (which is the survival prob and jump    */
/* at zi). n is included for convinience of C programing.               */

void wur(d, wt, sur, jum, n)    
double wt[], sur[], jum[];   /*  *wt, *sur, *jum[]? */
int d[], n;            /* should that be *n ? int *n ?  d[] ? */ 
                       /* since this is not directly interface with S */
                       /* may be we do not need them to be pointers   */
{  int i;
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

