/*************************************************************************
  Program Name : urnew010.c
  Programmer   : Li Lee and Mai Zhou (mai@ms.uky.edu) Gnu type copyright.
  Description  : On HP700 computers, compile by cc -c -O urnew010.c  
                 and then inside Splus use 
                 > dyn.load2("urnew010.o") to load it. 
                 Ported to R, change long to int to work on AMD64. 2/2005
  Date 1st ver. : June 9, 1993
  Last revision : June 12, 1995. 1999. 2002. called by R function d011()
  Depend on     : wur.c  and  ma.c
**************************************************************************/

/* #include <stdio.h>  */ 
/* #include <math.h>   */ 
#include <stdlib.h>

void urnew010(z, d, dup, sur, jum, max, err, r, s, rs, zext, dext, wext)
int *max, *r, *s, *rs;
int d[], dext[];
double *err;
double z[], sur[], jum[], zext[], wext[];
char *dup[];        /* I used char *dup and char dup[], they do not work! */ 

/* z, d, max and err are z012, d012, maxiter and error in Splus.
   dup, sur, jum  are as in Splus.
   Everything should be ordered according to (z,-d) before calling this
   function.
   r is the length of dup which >= length of z or d.
   s is the length of d. (and z)
   rs is length of d[d=2] 
   zext, dext, wext are for output extended data, they have length s+rs  */

{    int i, j, h, mm, nn, en= *r, n= *s, el= *rs, num, *k, *dadd, m= *max;
     int *d01;
     double u, *o, *w, *wadd, *zadd, *w01, *z01, *w2, a= *err;
     double ma(); 
     void wur();   

      k=(void *)malloc((el+1)*sizeof(int));
      o=(void *)malloc((n+1)*sizeof(double));
      w=(void *)malloc((n+1)*sizeof(double));
      wadd=(void *)malloc((n+el+1)*sizeof(double));
      zadd=(void *)malloc((n+el+1)*sizeof(double));
      dadd=(void *)malloc((n+el+1)*sizeof(int));

      for(i=0;i<n; i++) { 
           w[i]=1.0;
           o[i]=0.0;         /* do I need to end the o[n]='\0'; ?*/
      } 
      o[n]= '\0';

                             /* this block is to compute w. After  */
        j=0;                 /* that, length(w) should = length(z) */
        for(i=1; i<en; i++)
            if(*dup[i] == 'T')  /* am I compare pointer with int? dup[i]... */
               w[j]+=1;         /* SB-95                                    */
            else
               ++j; 
            
        w[j+1]= '\0'; 
 

/* Now add a phony point of d=-1 w=0 when the d= 0, 2 pattern happen */
/* the length of this amended vector is mm                          */
 
       j=0;
       for(i=0; i<n-1; i++) {
           zadd[j]=z[i];
           wadd[j]=w[i];
           dadd[j]=d[i];
           if(d[i]==0 && d[i+1] ==2) {
               ++j;
               zadd[j]=(z[i]+z[i+1])/2;
               wadd[j]=0.0;
               dadd[j]=(-1);   /* use dadd[i]=-1 to indicate this is foney */
           }
        ++j;
       }
      zadd[j]=z[n-1];
      wadd[j]=w[n-1]; 
      dadd[j]=d[n-1];
      zadd[j+1]= '\0';
      wadd[j+1]= '\0';
      dadd[j+1]= '\0'; 
      mm=j+1;


                         /* this block is to compute k .                */
        i=0;             /* k[i] is the # of d=0,1,-1 obs that are less */
        j=0;             /* the ith d=2 obs.                            */
        for(h=0; h<mm; h++)
            if(dadd[h]<2)
                ++j;
            else {
                k[i]=j;     /*  k[i++]=j; ? */ 
                ++i;
            }
        k[i]= '\0';
 
      
      w01=(void *)malloc((mm-el+1)*sizeof(double));
      z01=(void *)malloc((mm-el+1)*sizeof(double));
      d01=(void *)malloc((mm-el+1)*sizeof(int));
      w2=(void *)malloc((el+1)*sizeof(double));


                            /* This block is to separate data according */
         i=0;               /* to d value: d=0,1,-1 and d=2.  */
         j=0;               /* the length of z01 is nn */
         for(h=0; h<mm; h++)
             if(dadd[h]!=2) {
                z01[i]=zadd[h];
                d01[i]=dadd[h];
                w01[i]=wadd[h];
                ++i;
             } 
             else {
                w2[j]=wadd[h];    /* w2[j++]=wadd[h] ?? */ 
                ++j;
             }
          z01[i]= '\0';
          d01[i]= '\0';
          w01[i]= '\0';
          w2[j]= '\0';
          nn=i; 
          

/* Now compute an initial F(), we chose F to be the KM estimator */
/* with the phony point's weight elevate to 0.5, others gets 0.5 more too */

      for(j=0; j<nn; j++)
          w[j]=w01[j]+0.5;

      w[nn]= '\0';
      sur[nn]= '\0';
      jum[nn]= '\0';
      wur(d01, w, sur, jum, nn);
/* This is just one way of getting the initial F()            */ 

/* Now this is the main E-M iteration computations, surprisingly short */
/* compared to the other preparetory computations.                     */
      num=1;
      while((num<= m)  && (a < ma(sur, o, nn))) {
           for(j=0; j<nn; ++j){
                w[j]= w01[j];
                o[j]= sur[j];
           }
           for(i=0; i<el; ++i) {
                u= 1.0-sur[k[i]-1];
                for(j=0; j<k[i]; ++j)
                       w[j]+= w2[i]*(jum[j])/u;
           }
           wur(d01, w, sur, jum, nn);
          ++num;
      }

/* getting ready to exit: make sure the variables get their correct */
/* values, and free up spaces for the intermediate variables        */
      for(i=0;i<nn;i++) { 
          z[i]=z01[i];
          d[i]=d01[i];
      }
      *err= ma(sur, o, nn); 
      *max=num;
      *s=nn;
      for(i=0;i<mm;i++) {
          zext[i]=zadd[i];
          dext[i]=dadd[i];
          wext[i]=wadd[i];
      } 
      *r=mm;
      free(k);
      free(wadd);
      free(zadd);
      free(dadd);
      free(w01);
      free(z01);
      free(d01);
      free(w2);
      free(o);
      free(w); 
}
           

/*****************************************************************
copyright The software can be freely used and freely distributed
          for non-commercial purposes. Please send comments, bug 
          report etc. to mai@ms.uky.edu  or
 
            Mai Zhou
            Department of Statistics
            University of Kentucky
            Lexington, KY 40506-0027
*******************************************************************/
