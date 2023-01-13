#ifndef __DBLCENS_H_
/*****************************************************************
copyright GPL2. Please send comments, bug 
          report etc. to mai@ms.uky.edu  or
 
            Mai Zhou
            Department of Statistics
            University of Kentucky
            Lexington, KY 40506-0027
*******************************************************************/

#include <math.h>

#define INT_MALLOC_(L) (int * ) malloc(sizeof(int)*(L))
#define DOUBLE_MALLOC_(L) (double * ) malloc(sizeof(double)*(L))


/** This function computes the loglikihood for general self-consistency distribution.
 * @param dist double array
 * @param d censoring indicator
 * @param jump jump size
 * @param n sample size
 * copyright:The software is Gnu, means it can be freely used and freely 
 * distributed for non-commercial purposes. Please send comments,
 * bug report etc. to mai@ms.uky.edu  or
 * Mai Zhou
 *  Department of Statistics
 *  University of Kentucky
 *  Lexington, KY 40506-0027
*/
double loglik1(int * d, double * sur, double * jum, int n);

/** This function computes loglikelihood for constrained self-consistency distribution.
 * @param dist double array
 * @param d censoring indicator
 * @param jump jump size
 * @param n sample size
*/
double loglik2(double * dist, int * d, double * jump, int n);


/** The following function wur is to compute weighted Kaplan-Meier.      
 * Although z did not appear, the data, d, wt, etc. are                
 * sorted according to z value.                             
 * It computes the weighted Kaplan-Meier estimator, did not            
 * remove the censored observation (i.e. leave the observations        
 * as is. which is z, d, and wt. where wt do not have to sum to one.   
 * The censoring indicator, d[i], do not have to be 0/1 value,         
 * in fact, the only thing matters is wheather d[i]=0 or not.   
 * 
 * The results are in sur, jum. (which is the survival prob and jump  
 * at zi). n is included for convinience of C programing.               
 * 
 * @param d Censoring indicator
 * @param wt weight
 * @param sur Survival Probabilities
 * @param jum Jump at z_i
 * @param n sample size
 * */
void wur(int * d, double * wt, double * sur, double * jum, int n);


/** The following function ma is to calculate the maximum 
 * absolute value among the first n entries of the array s-t (both double),
 * where s and t are arrays of length n. 
 * @param t double array of length n
 * @param s double array of length n
 * @param n sample size
 * */
double ma(double * t, double * s, int n);


/** selfafterT
 *  Compute the NPMLE after time T         
 * 
 * copyright:The software is Gnu, means it can be freely used and freely 
 * distributed for non-commercial purposes. Please send comments,
 * bug report etc. to mai@ms.uky.edu  or
 * Mai Zhou
 *  Department of Statistics
 *  University of Kentucky
 *  Lexington, KY 40506-0027
 * 
 * @author: Kun Chen 
*/
void selfafterT(
  double * z, int * d, double *  wt, double * sur, double * extdist1,
  double * extdist2, double K, double konst, int n,
  int count, int m, int size1, double a);


/** selfbeforeT
 *  Compute the NPMLE before time T         
 * 
 * copyright:The software is Gnu, means it can be freely used and freely 
 * distributed for non-commercial purposes. Please send comments,
 * bug report etc. to mai@ms.uky.edu  or
 * Mai Zhou
 *  Department of Statistics
 *  University of Kentucky
 *  Lexington, KY 40506-0027
 * 
 * @author: Kun Chen 
*/
void selfbeforeT(
  double * z, int *d,double * wt, double * dist, double * extdist1, 
  double K, double konst, int n, int count, int m, int size1, double a);


/** z, d, max and err are z012, d012, maxiter and error in Splus.
 * dup, sur, jum  are as in Splus. Everything should be ordered 
 * according to (z,-d) before calling this function.
 * 
 * - r is the length of dup which >= length of z or d.
 * - s is the length of d. (and z)
 * - rs is length of d[d=2] 
 * - zext, dext, wext are for output extended data, they have length s+rs 
 * */
void urnew010(
    double * z, int * d, int *dup, double * sur, double * jum, int * max, 
    double * err, int * r, int * s, int * rs, 
    double * zext, int * dext, double * wext);



/** z, d, max and err are z012, d012, maxiter and error in Splus. dup, sur, jum are as in Splus.
 * Everything should be ordered according to (z,-d) before calling this function.
 * 
 * - r is the length of dup which >= length of z or d.
 * - s is the length of d. (and z)
 * - rs is length of d[d=2] 
 * - zext, dext, wext are for output extended data, they have length s+rs  */
void d011ch(
  double * z, 
  int * d, 
  int * dup, 
  double *sur, 
  double *jum, 
  int * max, 
  double * err, 
  int * r, 
  int * s, 
  int * rs, 
  double *zext, 
  int * dext, 
  double *wext,
  double * pK, 
  double * pkonst, 
  double *konstdist, 
  double *konstjump, 
  double * llratio);


#define __DBLCENS_H_
#endif
