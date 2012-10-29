\name{dblcens-internal}
\alias{Wdataclean2}
\title{Internal dblcens functions}
\usage{
Wdataclean2(z,d,wt=rep(1,length(z)))
}
\description{ Internal dblcens functions }
\details{
These are not intended to be called by the user. 

\code{Wdataclean2} will sort the data and collaps those that are true ties,
and the number of tied value is in the weights. 
Same code as in the package emplik.

}
\keyword{internal}
