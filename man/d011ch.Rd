\name{d011ch}
\alias{d011ch}
\title{Compute NPMLE of CDF from doubly censored data, with a constraint}
\description{
    \code{d011ch} computes the NPMLE of CDF, with and without a 
              constraint, from doubly censored data.
         It uses EM algorithm starting from an initial estimator that
         have jumps at uncensored points as well as the mid point of those
 censoring times that has pattern of (0,2), (see below for definition 
 and example.)

When there are ties, the left (right)
 censored points are treated as happened before (after), to break tie.
 Also the last right censored observation and first left censored
 observations are changed to uncensored, in order to obtain
 a proper distribution as estimator. (though this can be modified
 easily as they are written in R language).

It also computes the -2 log likelihood ratio based on the empirical
likelihood 
for testing the constraint via empirical likelihood theorem, i.e.
under Ho it should be chi-square with df=1.
}
\usage{
 d011ch(z, d, K, konst, identical=rep(0,length(z), maxiter=49, error=0.00001)
}
\arguments{
    \item{z}{a vector of length n denoting observed times, (ties permitted)}
    \item{d}{a vector of length n that contains censoring indicator:
            d= 2 or 1 or 0, (according to z being left, not, right censored)}
    \item{K}{the constraint time.}
    \item{konst}{the constraint value, i.e. F(K)=konst.}
    \item{identical}{optional. a vector of length n that has values 
               either 0 or 1.
             identical[i]=1 means even if $(z[i],d[i])$ is identical with
              $(z[j],d[j])$, for some $j \not= i$, they still stay as 2 
           observations,
 not 1 obs. with weight 2, which only happen if identical[i]=0
 and identical[j] =0. One reason for this is because they may have
 different covariates not shown here. This flexibility may be useful
 for regression applications. Default value is  identical = 0. }
    \item{maxiter}{optional integer value. default to 49}
    \item{error}{optional. Default to 0.00001}
}
\value{ 
    a list contain the NPMLE of CDF with and without the constraint,
 -2loglik ratio and other informations.
\item{time}{ survival times. Those corresponding to d=2 are removed.
Those corresponding to (0,2) censoring pattern are added, at mid-point.}
\item{status}{Censoring status of the above times. Since left censored 
times are removed, there is no status =2. There may be -1, indicating 
that this is an added time for (0,2) censoring pattern.
}
\item{surv}{The survival function at the above times.}
\item{jump}{Jumps of NPMLE at the above times.}
\item{exttime}{Similar to time but now include the left censored times. }
\item{extstatus}{Censoring status of exttime. -1 has same meaning as status
before.}
\item{extjump}{Jumps of the unconstrained NPMLE on extended times.}
\item{extsurv.Sx}{Survival probability at exttime.}
\item{konstdist}{The constrained NPMLE of distribution.}
\item{konstjump}{Jumps of the constrained NPMLE of CDF.}
\item{konsttime}{Location of the constraint, same as K in the input.}
\item{theta}{is the same value \code{konst} in the input.}
\item{"-2loglikR"}{the Wilks statistics. 
                   Distributed approximately chi-square df=1 under Ho}
\item{maxiter}{the actual number of iterations for the unconstrained NPNLE.
         The constrained NPMLE usually took less iterations to converge.}
}
\author{Kun Chen, 
Mai Zhou \email{mai@ms.uky.edu} }
\references{
Chang, M. N. and Yang, G. L. (1987). Strong consistency
of a nonparametric estimator of the survival function with doubly
censored data. Ann. Statist. 15, 1536-1547.

Murphy, S. and Van der Varrt. (1997).
Semiparametric Likelihood Ratio Inference.
Ann. Statist. 25, 1471-1509.

Chen, K. and Zhou, M. (2000). 
Nonparametric Hypothesis Testing and
Confidence Intervals with Doubly Censored Data. Tech Report, 
Univ. of Kentucky.
}
\examples{
d011ch(z=c(1,2,3,4,5), d=c(1,0,2,2,1), K=3.5, konst=0.6)
#
# you should get something like
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
#       [1] 1.0 2.0 2.5 3.0 4.0 5.0       (exttime include all the times,
#                                         censor or not, plus the added time)
#       $extstatus
#       [1]  1  0 -1  2  2  1
#
#       $extjump
#       [1] 0.4999649 0.0000000 0.1667174 0.0000000 0.0000000 0.3333177
#
#       $extsurv.Sx
#       [1] 0.5000351 0.5000351 0.3333177 0.3333177 0.3333177 0.0000000
#
#       $konstdist
#       [1] 0.4999365 0.4999365 0.6000000 0.6000000 0.6000000 1.0000000
#
#       $konstjump
#       [1] 0.4999365 0.0000000 0.1000635 0.0000000 0.0000000 0.4000000
#
#       $konsttime
#       [1] 3.5
#
#       $theta
#       [1] 0.6
#
#       $"-2loglikR"                  (the Wilks statistics to test Ho:
#       [1] 0.05679897                  F(K)=konst)
#
#       $maxiter
#       [1] 33
}
\keyword{survival}
\keyword{nonparametric}
