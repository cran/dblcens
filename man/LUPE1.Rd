\name{LUPE1}
\alias{LUPE1}
\title{Compute self-consistant estimator of CDF from doubly censored data}
\description{
    \code{LUPE1} computes the self-consistent estimator of CDF from doubly censored data
         via an algorithm specified in LUPE 2008.
         It only have jumps on the observed failure times (where d=1).

  When there are ties in the observed times, the left (right)
 censored points are treated as happened before (after) the failure time, 
 to break tie.
}
\usage{
LUPE1(survtime, status, maxiter=30, eps=1e-8)
}
\arguments{
    \item{survtime}{a vector of length N denoting observed times, 
                   (ties permitted)}
    \item{status}{a vector of length N that contains censoring indicator:
            d= 2 or 1 or 0, (according to z being left, not, right censored)}
    \item{maxiter}{optional integer value. default to 30}
    \item{eps}{optional. control iteration. Default to 1e-8}
}
\value{ 
    a list contain the self-consistant estimator of CDF and other informations.
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
}
\author{Mai Zhou \email{mai@ms.uky.edu} }
\references{
Chang, M. N. and Yang, G. L. (1987). Strong consistency
of a nonparametric estimator of the survival function with doubly
censored data. Ann. Statist. 15, 1536-1547.

LUPE (2009)
}
\examples{
LUPE1(survtime=c(1,2,3,4,5), status=c(1,0,1,2,1))
}
\keyword{survival}
\keyword{nonparametric}
