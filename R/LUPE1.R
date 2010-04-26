########################################
#  survtime = observed times, 
#  status = 0/1/2 for right/none/left censor
#  Now only work for d1=dN=1 case.
########################################
LUPE1 <- function(survtime, status, maxiter=30, eps=1e-8) {
#######################
## checking inputs   ##
#######################
survtime <- as.vector(survtime)
N <- length(survtime)
if(N <=1) stop("only one observation?")
if(length(status)!=N) stop ("length of status and survtime must agree")
if(any((status !=0)&(status !=1)&(status !=2)))
   stop("status must be 0(right-censored) or 1(uncensored) or 2(left-censored)")
if(!is.numeric(survtime)) stop("survtime must be numeric")
#####################################
## sort data according to survtime ##
#####################################
sortdata <- Wdataclean2(z=survtime, d=status)
exttime <- tt <- sortdata$value
extstatus <- d <- sortdata$dd
extweight <- we <- sortdata$weight

tt1 <- tt[d == 1]
dd1 <- we[d == 1]
n <- length(tt1) 
if( n < 2 ) stop("only 1 distinct uncensored observation")
########################
## compute mu and lam ##
########################
mu <- lam <- rep(0, n)
ind <- as.integer( (1:length(d))[ d == 1 ] )
for( k in 2:n ) {
wei <- we[ ind[k-1]:ind[k] ]
di  <-  d[ ind[k-1]:ind[k] ]
lam[k] <- sum( wei[ di == 0 ])
mu[k] <- sum( wei[ di ==2 ])
}
################### what about lam[1] and mu[1]???
## Initialize  ## 
###################
temp1 <- rep(1, n)/n
iternum <- 0
error <- N

#print( tt1 )
#print( dd1)
#print(lam)
#print(mu)
######################################################################
##  N is the sample size, n is the number of distinct uncensored times
######################################################################
##  Finally the iterations
#############################
while( iternum <= maxiter & error > eps) {
     temp2 <- LUPEiteration(d=dd1, mu=mu, lam=lam, jump=temp1, N=N)
     iternum <- iternum +1
     error <- sum(abs(temp2 - temp1))
     temp1 <- temp2
}

extjump <- rep(0, length(extstatus) )
extjump[extstatus == 1] <- temp2
extsurv <- 1 - cumsum(extjump)

dzero <- as.numeric(extstatus == 0)
jumpzero <- (extweight * dzero)/(N * extsurv)
jumpzero[is.na(jumpzero)] <- 0      
survzero <- 1 - cumsum(jumpzero)
  

dtwo <- as.numeric(extstatus == 2)
jumptwo <- (extweight * dtwo)/(N*(1 - extsurv))  
jumptwo[is.na(jumptwo)] <- 0           
survtwo <- rev(cumsum(rev(jumptwo)))


list(time=tt1, prob=temp2, iterations=iternum, L1error=error, 
     exttime=exttime, extjump1=extjump, extS1=extsurv,
     extjump0=jumpzero, extS0=survzero,
     extjump2=jumptwo,  extS2=survtwo)
}

######################################################
## since it oscilates, we do 2 iterations at a time ##
######################################################
LUPEiteration <- function(d, mu, lam, jump, N) {
jump <- jump/sum(jump)
jump1 <- jump2 <- jump
jump1 <- lupe(d=d, mu=mu, lam=lam, jump=jump, N=N )
jump2 <- lupe(d=d, mu=mu, lam=lam, jump=jump1, N=N )
return( (jump1+jump2)/2 )
}

lupe <- function(d, mu, lam, jump, N) {
 # d[1] must be > 0. All other d[i] must >0
 # mu[i] = sum mu in   ( tt(i-1), tt(i) ]  left censor
 # lam[i] = sum lam in [ tt(i-1), tt(i) )  right censor
 # therefore mu[1] and lam[1] never got used, just put zero there.
Sur <- rev( cumsum(rev(jump)) )
CDF <- 1-Sur

A <- lam/Sur - mu/CDF

f <- jump
f[1] <- d[1]/( N - sum( mu[-1]/CDF[-1] ) )

for(j in 2:length(jump)) {
f[j] <- (d[j]*f[j-1])/(d[j-1] - A[j]*f[j-1])
}
return( f/sum(f))
}

#################################################
## this is the same function in emplik package ##
#################################################
Wdataclean2 <- function (z, d, wt = rep(1, length(z))) 
{
    niceorder <- order(z, -d)
    sortedz <- z[niceorder]
    sortedd <- d[niceorder]
    sortedw <- wt[niceorder]
    n <- length(sortedd)
    y1 <- sortedz[-1] != sortedz[-n]
    y2 <- sortedd[-1] != sortedd[-n]
    y <- y1 | y2
    ind <- c(which(y | is.na(y)), n)
    csumw <- cumsum(sortedw)
    list(value = sortedz[ind], dd = sortedd[ind], weight = diff(c(0, 
        csumw[ind])))
}

