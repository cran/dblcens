/*************************************************************************
  Program Name : File name: d011ch.c Function name: d011ch [9/23/99email]
  Programmer   : Li Lee, Mai Zhou (mai@ms.uky.edu) and 
                 Kun Chen(chen@ms.uky.edu) 
  Copyright    : Gnu GPL2 type copyright.
  Description  : This program is used in R function d011ch, replacing 
                 (old) d009.
  First ver.   : June 9, 1993
  Revisions    : 6/27/1999, 2/2002. This C program should be used with
                 the Splus program d011ch.s. 
                 2/2005 changed type long to int to work better in AMD64.
**************************************************************************/

#include <math.h>
#include <stdlib.h>

void d011ch(z, d, dup, sur, jum, max, err, r, s, rs, zext, dext, wext,
             pK, pkonst, konstdist, konstjump, llratio)
int *max, *r, *s, *rs;
int d[], dext[];
double *err, *pK, *pkonst, *llratio;
double z[], sur[], jum[], zext[], wext[];
double konstdist[], konstjump[];
char *dup[];        /* I used char *dup and char dup[], they do not work! */ 


/* z, d, max and err are z012, d012, maxiter and error in Splus.
   dup, sur, jum  are as in Splus.
   Everything should be ordered according to (z,-d) before calling this
   function.
   r is the length of dup which >= length of z or d.
   s is the length of d. (and z)
   rs is length of d[d=2] 
   zext, dext, wext are for output extended data, they have length s+rs  */

{
  int i, j, h, mm, nn, en= *r, n= *s, el= *rs, num, *k, *dadd, m= *max;
  int *d01;
  double u, *o, *w, *wadd, *zadd, *w01, *z01, *w2, a= *err;

  double ma(); 
  void wur(); 

  void selfbeforeT(), selfafterT();
  double loglik1(), loglik2(); 

  double *jadd, *fadd, *sadd;
  int count=0, size1=0, cn1=0,cn2=0;
  double *extdist1, *extdist2;
  double K, konst, logliksc, loglikct, loglikratio;

  K = *pK;
  konst = *pkonst;

  k=(void *)malloc((el+1)*sizeof(int));
  o=(void *)malloc((n+1)*sizeof(double));
  w=(void *)malloc((n+1)*sizeof(double));
  wadd=(void *)malloc((n+el+1)*sizeof(double));
  zadd=(void *)malloc((n+el+1)*sizeof(double));
  dadd=(void *)malloc((n+el+1)*sizeof(int));

  /* This block is to compute w. After that, length(w) should =
     length(z).
  */
  for(i=0;i<n; i++)
  { 
    w[i]=1.0;
    o[i]=0.0;
  } 
  o[n]= '\0';

  j=0;
  for(i=1; i<en; i++)
  {
    if(*dup[i] == 'T')
      w[j]+=1;
    else
      ++j; 
  }
           
  w[j+1]= '\0'; 
 
  /* Now add a phony point of d=-1 w=0 when the d=0, 2 pattern 
     happen. The length of this amended vector is mm.
  */
  j=0;
  for(i=0; i<n-1; i++) 
  {
    zadd[j]=z[i];
    wadd[j]=w[i];
    dadd[j]=d[i];
    if(d[i]==0 && d[i+1] ==2)
    {
      ++j;
      zadd[j]=(z[i]+z[i+1])/2;
      wadd[j]=0.0;
      dadd[j]=(-1);
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

  /* This block is to compute k. k[i] is the # of d=0,1 -1 obs
     that are less than ith d=2 obs.
  */
  i=0;
  j=0;
  for(h=0; h<mm; h++)
  {
    if(dadd[h]<2)
      ++j;
    else {
      k[i]=j;
      ++i;
    }
  }
  k[i]= '\0';
 
      
  w01=(void *)malloc((mm-el+1)*sizeof(double));
  z01=(void *)malloc((mm-el+1)*sizeof(double));
  d01=(void *)malloc((mm-el+1)*sizeof(int));
  w2=(void *)malloc((el+1)*sizeof(double));

  /* This block is to seperate data according to d vlaue:
     d=0, 1, -1 and d=2. The length of z01 is nn 
  */
  i=0;
  j=0;
  for(h=0; h<mm; h++)
    if(dadd[h]!=2)
    {
      z01[i]=zadd[h];
      d01[i]=dadd[h];
      w01[i]=wadd[h];
      ++i;
    } else {
      w2[j]=wadd[h];
      ++j;
    }
    z01[i]= '\0';
    d01[i]= '\0';
    w01[i]= '\0';
    w2[j]= '\0';
    nn=i; 
   
    /* Now compte an initial F(), we chose F to be the KM estimator
       with the phony point's weight elevate to 0.5, others gets 0.5 
       more too.
    */      
    for(j=0; j<nn; j++)
      w[j]=w01[j]+0.5;

    w[nn]= '\0';
    sur[nn]= '\0';
    jum[nn]= '\0';
    wur(d01, w, sur, jum, nn);
   /* this is just one way of getting the initial F()   */


   /* Now this is the main E-M iteration computations, surprisingly
      short compared to the other preparetory computations.
   */
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

    for(i=0;i<nn;i++)
    { 
      z[i]=z01[i];
      d[i]=d01[i];
    }
      
    *err= ma(sur, o, nn); 
    *max=num;
    *s=nn;
    
    jadd=(void *)malloc((mm+1)*sizeof(double));
    fadd=(void *)malloc((mm+1)*sizeof(double));
    sadd=(void *)malloc((mm+1)*sizeof(double));

    /* This block gets extended time, status, weights and jumps */ 
    j=0; 
    for(i=0;i<mm;i++) 
    {
      zext[i]=zadd[i];
      dext[i]=dadd[i];
      wext[i]=wadd[i];
      if(dadd[i]!=2)
      {
        jadd[i]=jum[j];
        j++;
      }
      else jadd[i]=0;
    } 
    jadd[mm]='\0';
 
    /* After getting jumps, we can obtain extended distribution
       and then survival funtion. 
    */ 
    fadd[0]=jadd[0]; 
    sadd[0]=1-fadd[0]; 
    for(i=1; i<mm;i++)
    { 
      fadd[i]=fadd[i-1]+jadd[i];
      sadd[i]=1-fadd[i];
    }
    fadd[mm]='\0';
    sadd[mm]='\0';
 
    /* count how many observations are less than or equal to 
       the constrained time. */
    for(i=0;i<mm;i++)
    {
      if(zadd[i] <= K)
        count++;
    }

    for(i=0;i<count;i++)
    {
      if(abs(dadd[i])==1)
        size1++;
    }                                

    extdist1 = (void *)malloc((count+1)*sizeof(double));
    extdist2 = (void *)malloc((mm-count+1)*sizeof(double));
   
    /* Compute the modify self-consistency distribution before 
       and after constrained time T.  */
    for(i=0;i<mm;i++)
    {   
      if(zadd[i] <= K)
        cn1++;
      if(zadd[i] > K)
        cn2++;
    }

    if(cn2==0)
    {
      selfbeforeT(zadd,dadd,wadd,fadd, extdist1, K,konst,mm,count,m,size1,a);
      for(i=0;i<mm;i++)
        konstdist[i]=extdist1[i];
      konstjump[0]=konstdist[0]; 
    } 
    else if(cn1==0)
    {
      selfafterT(zadd,dadd,wadd,sadd, extdist1, extdist2, K,konst,mm,
                 count,m,size1,a); 
      for(i=0;i<mm;i++)
      {
        konstdist[i]=extdist2[i];
        konstjump[0]=konstdist[0]-konst;
      }
    } 
    else
    { 
      selfbeforeT(zadd,dadd,wadd,fadd, extdist1, K,konst,mm,count,m,size1,a);
      selfafterT(zadd,dadd,wadd,sadd, extdist1, extdist2, K,konst,mm,
                 count,m,size1,a); 
 
      /* getting constrained jumps from constrained distribution. */ 
      for(i=0;i<count;i++)
      {
        konstdist[i]=extdist1[i];

      }
      for(i=count;i<mm;i++)
      {
        konstdist[i]=extdist2[i-count];
      }
      konstjump[0]=konstdist[0]; 
    }
    for(i=1;i<mm;i++)
    {  
      konstjump[i]=konstdist[i]-konstdist[i-1];
/****************************************************************
      printf("konstdist[%d] = %lf\n", i, konstdist[i]);   
      printf("dadd[%d] = %ld\n", i, dadd[i]);   
      printf("knostjump[%d] = %lf\n", i, konstjump[i]);   
*****************************************************************/
    }    
    konstdist[mm]='\0';
    konstjump[mm]='\0';
   
    /* Compute the loglikelihood for general self-consisteny
       distribution and constrained self-consistency dist.
    */  
    logliksc=loglik1(dadd,sadd,jadd,mm);
    loglikct=loglik2(konstdist,dadd,konstjump,mm); 
/*    
    printf("logliksc=%lf\n", logliksc);  
    printf("loglikct=%lf\n", loglikct);
*/
    /* Then compute the -2loglikilihood ration. */
    loglikratio=2*(logliksc-loglikct);
/*    printf("loglikratio=%lf\n", loglikratio); */ 
 
    /* Ready to exit. Make sure the variable get their correct
       values, and free up spaces for the intermediate variables.*/
    
    *r=mm;
    *llratio = loglikratio;
 
    free(extdist1);
    free(extdist2);
    free(jadd);
    free(sadd);
    free(fadd); 
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
copyright:The software is Gnu, means it can be freely used and freely 
          distributed for non-commercial purposes. Please send comments,
          bug report etc. to mai@ms.uky.edu  or
 
            Mai Zhou
            Department of Statistics
            University of Kentucky
            Lexington, KY 40506-0027
*******************************************************************/
/* This function computes the loglikihood for general  */
/* self-consistency distribution.                      */
      
double loglik1(d,sur,jum,n)
int d[],n;
double sur[], jum[];
{
  double *sur2, *sur0, *sur1;
  double sum1, sum2, sum0;
  double logliksc;
  int k, l, i, j, m;

  k=l=i=j=m=0; 
  sum1 = sum2 = sum0 = 0.0; 
  for(i=0;i<n;i++)
  {
    if(d[i]==2)
      k++; 
    if(d[i]==0)
      l++;
    if(d[i]==1)
      m++;
  }

  sur2=(void *)malloc((k+1)*sizeof(double));
  sur0=(void *)malloc((l+1)*sizeof(double));
  sur1=(void *)malloc((m+1)*sizeof(double));

  j=0; 
  for(i=0;i<n;i++)
  {
    if(d[i]==2)
    {
      sur2[j]=1-sur[i]; 
/*      printf("sur2[%d]=%lf\n",j,sur2[j]);   */
      sum2+=log(sur2[j]); 
      j++;
    }
  }
  sur2[k]='\0';

  j=0;
  for(i=0;i<n;i++)
  {
    if(d[i]==0)
    {
      sur0[j]=sur[i];
/*      printf("sur0[%d]=%lf\n",j,sur0[j]); */
      sum0+=log(sur0[j]);
      j++;
    }
  }
  sur0[l]='\0';

  j=0;
  for(i=0;i<n;i++)
  {
    if(d[i]==1)
    {
      sur1[j]=jum[i];
      sum1+=log(sur1[j]);
      j++;
    }
  }
  sur1[m]='\0';

  logliksc=sum0+sum1+sum2; 

  free(sur1);
  free(sur2);
  free(sur0);

  return(logliksc);
}

/* This function computes loglikelihood for constrained    */
/* self-consistency distribution.                          */

double loglik2(dist, d,jump,n)
double dist[], jump[];
int d[],n;
{
  double *sur1, *sur0, *sur2;
  double sum1, sum0, sum2; 
  double loglikct;  
  int i,j,k,l,m;

  sum1 = sum2 = sum0 =0.0;
  i = j = k = l = m =0;

  /* count how many right-censored obs. how many left-censored */
  /* obs. how many uncensored obs.                             */
  for(i=0;i<n;i++)
  {
    if(d[i]==2)
      k++; 
    if(d[i]==0)
      l++;
    if(d[i]==1)
      m++;
  }
  sur2=(void *)malloc((k+1)*sizeof(double));
  sur0=(void *)malloc((l+1)*sizeof(double));
  sur1=(void *)malloc((m+1)*sizeof(double));

  j=0; 
  for(i=0;i<n;i++)
  {
    if(d[i]==2)
    {
      sur2[j]=dist[i];   
/*      printf("sur2[%d]=%lf\n",j,sur2[j]);   */
      sum2+=log(sur2[j]);
      j++;
    }
  }
  sur2[k]='\0';

  j=0;
  for(i=0;i<n;i++)
  {
    if(d[i]==0)
    {
      sur0[j]=1-dist[i];
/*      printf("sur0[%d]=%lf\n",j,sur0[j]); */
      sum0+=log(sur0[j]);
/*    printf("sur0[%d]=%lf\n",j,sur0[j]); */
      j++;
    }
  }
  sur0[l]='\0';

  j=0;
  for(i=0;i<n;i++)
  {
    if(d[i]==1)
    {
      sur1[j]=jump[i];
/*      printf("sur1[%d]=%lf\n",j,sur1[j]);  */
      sum1+=log(sur1[j]);
      j++;
    }
  }
  sur1[m]='\0';

  loglikct=sum0+sum1+sum2;  
  
  free(sur1);
  free(sur2);
  free(sur0);
  
  return(loglikct);
}


/*****************************************************************
copyright:The software is Gnu, means it can be freely used and freely 
          distributed for non-commercial purposes. Please send comments,
          bug report etc. to mai@ms.uky.edu  or
 
            Mai Zhou
            Department of Statistics
            University of Kentucky
            Lexington, KY 40506-0027
*******************************************************************/

/****************************************/
/* Compute the NPMLE after time T */
/* written by Kun Chen            */
/* depend on ma()                 */

void selfafterT(z,d,wt,sur, extdist1, extdist2, K,konst,n,count,m,size1, a)
double z[], wt[], sur[], extdist1[],extdist2[], konst, K, a;
int d[], n, count,m,size1;
{
  double value;
  int count2, size2, w2, *newd;
  double *newtime, *theta, *Ftheta;
  double *newdist, *newFtheta, *neww;
  int num, q;
  double sum1, sum2, sum3, sum4;
  int i, j, h, l, k;
  double ma();

  count2 = size2 = w2 = 0;
  num = 1;
  q = 0;
  h = 1;
  value=sur[0];
 
  for(i=0;i<n;i++)
  {
    if(z[i]>K)
      count2++;
  } 
  newtime=(void *)malloc((count2+1)*sizeof(double));
  newd=(void *)malloc((count2+1)*sizeof(int)); 
  newdist=(void *)malloc((count2+1)*sizeof(double));
  neww=(void *)malloc((count2+1)*sizeof(double));

  for(i=0;i<count2;i++){
    newtime[i]=z[count+i];
    newd[i]=d[count+i];
    neww[i]=wt[count+i];
    if(value != 0)
      newdist[i]=1-sur[count+i]*konst/value;
    else
      newdist[i]=1-sur[count+i];
/*    printf("newdist[%d]=%lf\n",i,newdist[i]);  */
    w2+=wt[count+i];
  }
  newtime[count2]='\0';
  newd[count2]='\0';
  newdist[count2]='\0';
  neww[count2]='\0';
 
  for(i=0;i<count2;i++){
    if(abs(newd[i])==1)
      size2++;
  } 

  theta=(void *)malloc((size2+2)*sizeof(double));
  Ftheta=(void *)malloc((size2+2)*sizeof(double));
  newFtheta=(void *)malloc((size2+2)*sizeof(double)); 

  if(size2 == 0)
  { 
    for(i=0;i<count2;i++)
    {
      extdist2[i] = extdist1[count-1];  
    } 
    extdist2[count2]='\0';
  }
  else
  { 
 
    theta[0]=K;
    Ftheta[0]=konst;
   
    if(size1 == 0)
      newFtheta[0]=0;
    else
      newFtheta[0]=konst;

    h=1;
    for(i=0;i<count2;i++)
    {
      if(abs(newd[i])==1)
      {
        theta[h]=newtime[i];
        Ftheta[h]=newdist[i]; 
        h++; 
      }
    } 
    theta[size2+1]='\0';
    Ftheta[size2+1]='\0';
  
    while(num<= m)
    {
      for(j=1;j<size2+1;j++)
      {
        sum1=0.0;
        sum2=0.0;
        sum3=0.0;
        sum4=0.0;
        for(i=0;i<count2;i++)
        {
          q=0;
          if((newtime[i]<=theta[j])&&((newd[i])==1))
            sum1 += neww[i];
          for(l=0;l<size2+1;l++)
          {
            if(theta[l]<=newtime[i])
              q++;
          } 
          if((newd[i]==0)&&(newtime[i]<=theta[j]))
            sum2+= neww[i]*(Ftheta[j]-Ftheta[q-1])/(1-Ftheta[q-1]);
          if(newd[i]==2)
          {
            sum3+=neww[i]*konst*(Ftheta[j]-konst)/(Ftheta[q-1]*(1-konst)); 
            if((newtime[i]<=theta[j]))
              sum4+= neww[i]*(1-konst/Ftheta[q-1]);
            else
              sum4+= neww[i]*(Ftheta[j]-konst)/Ftheta[q-1];
          }
        }
        newFtheta[j] = konst + ((1-konst)/w2)*(sum1+sum2+sum3+sum4);
      }
      if(a < ma(newFtheta, Ftheta, size2))
      {
        for(k=1;k<size2+1;k++)
        {
          Ftheta[k]=newFtheta[k];
        } 
        num++;
      }
      else break;
    }
    
    j=1;
    for(i=0;i<count2;i++)
    {
      if(abs(newd[i])!=1)
        extdist2[i]=newFtheta[j-1];
      else
      {
        extdist2[i]=newFtheta[j];
/*        printf("extdist2[%d]=%lf\n",i,extdist2[i]); */ 
        j++;
      }
    }
    extdist2[count2]='\0';
  } 

  free(newtime);
  free(newd);
  free(neww);
  free(newdist);
  free(theta);
  free(Ftheta);
  free(newFtheta);
}

/*****************************************************************
copyright:The software is Gnu, means it can be freely used and freely 
          distributed for non-commercial purposes. Please send comments,
          bug report etc. to mai@ms.uky.edu  or
 
            Mai Zhou
            Department of Statistics
            University of Kentucky
            Lexington, KY 40506-0027
*******************************************************************/
/* Compute the NPMLE before time T */
/* written by Kun Chen             */
/* depend on ma()                  */

void selfbeforeT(z,d,wt,dist, extdist1, K,konst, n, count, m, size1,a)
double z[], wt[], dist[], extdist1[], konst, K, a;
int d[],n, count,m, size1;
{
  double value, ans;
  int  *newd;
  double *newtime, *theta, *Ftheta;
  double *newdist ,*newFtheta, *extFtheta;
  int num=1, q=0, count1=0;
  double sum1, sum2, sum3, sum4;
  double numerator, w1=0.0;
  int h, i, j, l;
  double ma();

  value=dist[count-1];
 
  newtime=(void *)malloc((count+1)*sizeof(double));
  newd=(void *)malloc((count+1)*sizeof(int)); 
  newdist=(void *)malloc((count+1)*sizeof(double)); 

  for(i=0;i<count;i++){
    newtime[i]=z[i];
    newd[i]=d[i];
    w1+=wt[i];
    if(value != 0)
      newdist[i]=dist[i]*konst/value;
    else
      newdist[i]=dist[i];
  }
  newtime[count]='\0';
  newd[count]='\0';
  newdist[count]='\0';
  
  theta=(void *)malloc((size1+1)*sizeof(double));
  Ftheta=(void *)malloc((size1+1)*sizeof(double));
  newFtheta=(void *)malloc((size1+1)*sizeof(double));
  extFtheta=(void *)malloc((size1+2)*sizeof(double));
  if(size1 == 0)
  {  
    for(i=0;i<count;i++)
    {
      extdist1[i]=0;
    }
  }
  else
  {
    if((count <= 1) && (newd[0]==1))  
    {
      for(i=0;i<count;i++)
      { 
        extdist1[i]=konst;
      } 
    }
    else if(count >=2)
    { 
      for(i=1;i<count;i++)
      {
        if(abs(newd[i])==1)
          count1++;
      }
      if((newd[0]==1) && (count1==0))
      {
        for(i=0;i<count;i++)
        { 
          extdist1[i]=konst;
        } 
      }
      else
      {
        h=0; 
        for(i=0;i<count;i++)
        {
          if(abs(newd[i])==1)
          {
            theta[h]=newtime[i];
            Ftheta[h]=newdist[i]; 
/*
            printf("theta[%d]=%lf\n",h,theta[h]);
            printf("Ftheta[%d]=%lf\n",h,Ftheta[h]);
*/
            h++; 
          }
        } 
        theta[size1]='\0';
        Ftheta[size1]='\0';

        while(num<= m) 
        {
          for(j=0;j<size1;j++)
          {
            sum1=0.0;
            sum2=0.0;
            sum3=0.0;
            sum4=0.0;
            numerator=((1-konst)/konst)*Ftheta[j];
            for(i=0;i<count;i++)
            {
              q=0; 
              if((newtime[i]<=theta[j])&&(newd[i]==1))
              {
                sum1 += wt[i];
              }
              for(l=0;l<size1;l++)
              {
                if(theta[l]<newtime[i])
                q++; 
              }
              if(q!=0)
                ans = Ftheta[q-1];
              else
                ans = 0;  
              if(newd[i]==0)
              {
                sum2+= wt[i]*numerator/(1-ans);
              }
              if((newtime[i]<=theta[j])&&(newd[i]==0))
              {
                sum3+= wt[i]*(Ftheta[j]-ans)/(1-ans);
        
              } 
              if(newd[i]==2)
              {
                if((newtime[i]<=theta[j]))
                  sum4+= wt[i];
                else
                  sum4+= wt[i]*Ftheta[j]/ans;
              }
            }
            newFtheta[j] = (konst/w1)*(sum1+sum2+sum3+sum4);
          }
          newFtheta[size1]='\0';
    
          if(a < ma(newFtheta, Ftheta, size1))
          {
            for(i=0;i<size1;i++)
              Ftheta[i]=newFtheta[i];
            num++;
          }
          else break;
        }
  
        if(newd[0] == 1)
        {
          extdist1[0]=newFtheta[0];
          for(i=0;i<size1;i++)
            extFtheta[i]=newFtheta[i];
          extFtheta[size1+1]='\0';
        }
        else
        {
          extdist1[0]=0;
          extFtheta[0]=0;
          for(i=1;i<(size1+1);i++)
            extFtheta[i]=newFtheta[i-1];
          extFtheta[size1+1]='\0';
        } 
        j=1;
        if(count >= 2)
        {
          for(i=1;i<count;i++)
          {
            if(abs(newd[i])!=1)
              extdist1[i]=extFtheta[j-1];
            else 
            {
              extdist1[i]=extFtheta[j];
              j++;
            }
          } 
        } 
  
        extdist1[count]='\0';
      }
    }
  } 
  free(newtime);
  free(newd);
  free(newdist);
  free(theta);
  free(Ftheta);
  free(newFtheta);
  free(extFtheta);

}

/*****************************************************************
copyright:The software is Gnu, means it can be freely used and freely 
          distributed for non-commercial purposes. Please send comments,
          bug report etc. to mai@ms.uky.edu  or
 
            Mai Zhou
            Department of Statistics
            University of Kentucky
            Lexington, KY 40506-0027
*******************************************************************/
