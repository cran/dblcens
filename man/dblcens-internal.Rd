\name{dblcens-internal}
\alias{lupe}
\alias{LUPEiteration}
\alias{Wdataclean2}
\title{Internal dblcens functions}
\usage{
lupe(d, mu, lam, jump, N)
LUPEiteration(d, mu, lam, jump, N) 
Wdataclean2(z,d,wt=rep(1,length(z)))
}
\description{ Internal dblcens functions }
\details{
These are not intended to be called by the user. 

\code{lupe} is one iteration of the computation.

\code{LUPEiteration} is doing two iterations at a time and take the average
as output, to combat oscillations.

\code{Wdataclean2} will sort the data and collaps those that are true ties,
and the number of tied value is in the weights. 
Same code as in the package emplik.

\code{Wdataclean2} is used by LUPE1.
}
\keyword{internal}
