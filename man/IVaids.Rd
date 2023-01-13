\name{IVaids}
\alias{IVaids}
\title{Data: AIDS patient among IV drug user}
\description{
Time to AIDS among 232 patients infected with HIV. 
136 left AIDS-free (right censored). 14 died with AIDS without prior diagnoses (left censored). 
82 had AIDS while in the program (non-censored).
}
\usage{ data(IVaids) }
\details{
A date set with 232 rows and 7 variables.

"A doubly censoring scheme occurs when the lifetimes T being measured, from a well-known time origin, are exactly observed only within a window [L, R] of observational time and are otherwise censored either from above by R (right-censored observations) or below by L (left-censored observations). 

Sample data consists on the pairs (U, delta) where U = min[R, max(T, L)] and delta indicates whether T is exactly observed (delta = 0), right-censored (delta = 1) or left-censored (delta = -1). We are interested in the estimation of the marginal behaviour of the three random variables T, L and R based on the observed pairs (U, delta)." ------ quote from the below reference paper.

The definition of the censoring indicator, (delta), here is different from the one we use, (d), in the functions d011( ) and d011ch( ).

  - (delta=0) corresponds to (d=1);

  - (delta=1) corresponds to (d=0);

  - (delta= -1) corresponds to (d=2).

Therefore we need to make the change before call the function d011( ) or d011ch( ). 

"The data set is from a cohort of drug users recruited in a detoxification program in Badalona (Spain). For these data we may estimate the survival function for the elapsed time from starting IV-drugs to AIDS diagnosis, as well as the potential follow-up time." ------quote from the below reference paper. 


The entry "AIDSDate-FIRST_IV" is the observed AIDS-free times, U, in the above definition. The unit of measurement is "days".

According to the paper the estimated median of U is 15.44 years. And at time of 10 years, the probability of AIDS-free is about 0.7. (from my reading of the plot from the paper, this probability is more closer to 0.75).
}
\references{ 
Julia, Olga and Gomez, Guadalupe (2011)
Simultaneous marginal survival estimators when doubly censored data is present.
Lifetime Data Analysis,
July 2011, Volume 17, Issue 3, pp 347-372.
}
\keyword{datasets}
