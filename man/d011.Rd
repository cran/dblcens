\name{d011}
\alias{d011}
\title{Compute NPMLE of CDF from doubly censored data}
\description{
    \code{d011} computes the NPMLE of CDF from doubly censored data
         via EM algorithm starting from an initial estimator that
         have jumps at (1) uncensored points; (2) (mid-point of) 
         consecutive survival times with censoring indicator 
         pattern of (0,2), (see below for definition).

  When there are ties, the left (right)
 censored points are treated as happened slightly before (after), to break tie.
 Also when the last observation happens to be right censored and/or 
 when the first observation happens to be left censored, they are changed 
 to uncensored. This is to ensure we obtain
 a proper distribution as the CDF estimator. (though this can be modified
 easily as they are written in R language).

It also computes the NPMLE of the two censoring distributions.
There is an option that you may also try to compute the 
three influence functions (but could slow and memory hungry).
}
\usage{
d011(z, d, identical = rep(0, length(z)),
         maxiter = 49, error = 0.00001, influence.fun = FALSE)
}
\arguments{
    \item{z}{a vector of length n denoting observed times, (ties permitted)}
    \item{d}{a vector of length n that contains censoring indicator:
            d= 2 or 1 or 0, (according to z being left, not, right censored)}
    \item{identical}{optional. A vector of length n that has values 
               either 0 or 1.
             identical[i]=1 means: even if (z[i],d[i]) is identical to
              (z[j],d[j]), for some \eqn{j \not= i}, they still stay as 2 
          observations,
 (not 1 obs. with weight 2, which only happen if identical[i]=0
 and identical[j] =0). One reason for this is because they may have
 different covariates not shown here. This adds more flexibility
 for regression applications. 
 Default value is  identical = 0, (i.e. collapse
 if identical observations). }
    \item{maxiter}{optional integer value. default to 49}
    \item{error}{optional. Default to 0.00001}
    \item{influence.fun}{optional. Default to FALSE. If TRUE, the code will
 try to compute the influence functions (3 of them) at the censored times.
 This computation can be very slow and memory intensive (for data with
 >500 censored times). }
}
\value{ 
    a list contain the NPMLE of CDF and other information.
\item{time}{Times of input z, with time corresponding to status=2 removed.}
\item{status}{Censoring status of the above times. Status = -1 means this is
an added time because of the censoring pattern (0,2).}
\item{surv}{Survival probability at the above times.}
\item{jump}{Jumps of the NPMLE at the above times.}
\item{exttime}{Similar to times but those with status =2 not removed.}
\item{extstatus}{status of exttime}
\item{extjump}{jump pf NPMLE at exttime.}
\item{extsurv.Sx}{Estimated lifetime distribution.}
\item{surv0.Sy}{One of the censoring distributions.}
\item{jump0}{Jump of surv0.Sy}
\item{surv2.Sz}{Another censoring distribution.}
\item{jump2}{Jump of surv2.Sz}
\item{conv}{ A vector of length 2: the actual number of iterations, and the 
  actual error of successive iteration. If the iteration number equal to 
               the maxiter you set, then the iteration has not converged.}
\item{Nodes}{Points where the influence function is computed.}
\item{IC1tu}{Influence function value at the nodes. See Chang (1990)
for details.}
\item{IC1tu2}{Influence function value at other points. See Chang (1990)
for details.}
\item{IC2tu}{ditto IC1tu}
\item{IC3tu}{ditto IC1tu}
\item{VarFt}{Estimated variances of F(t) at the Nodes.}
}
\author{Mai Zhou \email{mai@ms.uky.edu}, Li Lee.}
\references{
Chang, M. N. and Yang, G. L. (1987). Strong consistency
of a nonparametric estimator of the survival function with doubly
censored data. Ann. Statist. 15, 1536-1547.

Turnbull (1976) The empirical distribution function with
arbitrarily grouped, censored and truncated data. JRSS B 290-295.

Chang, M. N. (1990). Weak convergence in doubly censored data.
Ann. Statist. 18, 390-405.

Chen, K. and Zhou, M. (2003). 
Nonparametric Hypothesis Testing and
Confidence Intervals with Doubly Censored Data.
Lifetime Data Analysis,  9, 71-91. 
}
\examples{
d011(z=c(1,2,3,4,5), d=c(1,0,2,2,1))
#
# you should get something like below (and more)
#
#       $time:
#       [1] 1.0 2.0 2.5 5.0    (notice the times, (3,4), corresponding
#                                   to d=2 are removed, and time 2.5 added
#       $status:               since there is a (0,2) pattern at
#       [1]  1  0 -1  1        times 2, 3. The status indicator of -1
#                                   show that it is an added time )
#       $surv
#       [1] 0.5000351 0.5000351 0.3333177 0.0000000
#
#       $jump
#       [1] 0.4999649 0.0000000 0.1667174 0.3333177
#
#       $exttime
#       [1] 1.0 2.0 2.5 3.0 4.0 5.0
#
#       $extstatus
#       [1]  1  0 -1  2  2  1
#
#       ...... 
#
#       $conv
#       [1] 3.300000e+01  8.788214e-06  ### did 33 iterations
#
# BTW, the true NPMLE of surv is (1/2, 1/2, 1/3, 0) at times (1,2,2.5,5).
###### Example 2. 
d011(c(1,2,3,4,5), c(1,2,1,0,1),influence.fun=TRUE)
#     we get
# ......
#$conv:
#[1] 3 0
#
#$Nodes:
#[1] 2 4
#
#$IC1tu:
#     [,1] [,2]
#[1,]   -1    0
#[2,]   -1   -2
#
#$IC2tu:
#           [,1] [,2]
#[1,]  0.0000000    0
#[2,] -0.3333333    0
#
#$IC3tu:
#     [,1]       [,2]
#[1,]   -1 -0.6666667
#[2,]   -1 -1.0000000
#
#$VarFt:
#[1] 0.24 0.24           ## est var of F(t) at t=nodes
#######################################################
}
\keyword{survival}
\keyword{nonparametric}
