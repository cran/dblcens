/* The following function ma is to calculate the maximum 
absolute value among the first n entries of the array s-t (both double),
where s and t are arrays of length n. */

#include <math.h>

double ma(t, s, n)
double t[], s[];     /*  *t, *s ? */ 
long n;              /* should that be *n ?  or long *n ? */ 
{
    long i;
    double a=0.0, b;

    for(i=0; i<n; i++) 
         if((b=fabs(t[i]-s[i])) >a)
              a=b;
/*     return a;   */
    return(a);
}

