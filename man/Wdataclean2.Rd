\name{Wdataclean2}
\alias{Wdataclean2}
\title{Internal dblcens functions}
\usage{
Wdataclean2(z, d, wt=rep(1,length(z)))
}
\description{ Internal dblcens functions }
\details{
These are not intended to be called by the user. 

\code{Wdataclean2} will sort the data and collaps those that are true ties,
and the number of tied value is in the weights. 
Same code as in the package emplik.

}

\arguments{
    \item{z}{a vector of length n denoting observed times, (ties permitted)}
    \item{d}{a vector of length n that contains censoring indicator:
            d = 2 or 1 or 0, (according to z being left, not, right censored)}
    \item{wt}{a vector of length n that is used to derive the number of ties. 
    By default it is a ``1'' vector of length n.}
}
\value{ 
    \item{value}{ Cleaned survival times.}
    \item{d}{ Cleaned Censoring status of the above times.}
    \item{weight}{a vector that indicate the number of ties.}
}

\examples{
z <- c(0.312 ,0.808 ,0.793 ,2.119 ,0.152 ,0.104 ,1.002 ,0.82 ,0.356 ,0.618)
d <- c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0)
Wdataclean2(z,d)
}
\keyword{Wdataclean2}
