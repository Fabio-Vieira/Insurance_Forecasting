data <- matrix(c(3323,8332,9572,10172,7631,3855,3252,4433,2188,
              	 3785,10342,8330,7849,2839,3577,1404,1721,NA,
              	 4677,9989,8746,10228,8572,5787,3855,NA,NA,
              	 5288,8089,12839,11829,7560,6383,NA,NA,NA,
              	 2294,9869,10242,13808,8775,NA,NA,NA,NA,
              	 3600,7514,8247,9327,NA,NA,NA,NA,NA,
              	 3642,7394,9838,NA,NA,NA,NA,NA,NA,
              	 2463,5033,NA,NA,NA,NA,NA,NA,NA,
              	 2267,NA,NA,NA,NA,NA,NA,NA,NA), byrow = T, ncol = 9)

test <- matrix(c(NA,NA,NA,NA,NA,NA,NA,NA,NA,
		             NA,NA,NA,NA,NA,NA,NA,NA,1065,
		             NA,NA,NA,NA,NA,NA,NA,1445,1612,
		             NA,NA,NA,NA,NA,NA,4118,3016,1575,
		             NA,NA,NA,NA,NA,5419,2424,1597,4149,
		             NA,NA,NA,NA,8584,4245,4096,3216,2014,
		             NA,NA,NA,9733,6377,4884,11920,4188,4492,
		             NA,NA,6980,7722,6702,7834,5579,3622,1300,
		             NA,5959,6175,7051,8102,6339,6978,4396,3107), byrow = T,
               ncol = 9)

N <- 9
Y <- NULL
for(i in 1:N){Y <- c(Y, data[1:(N-i+1),i])}

Z <- log(Y)
hist(Z)

################################################################################
############################Matching Functions##################################
################################################################################

VecAlpha <- function(Alpha,N){
  AlphaVec <- NULL
  for(i in 1:N){
    AlphaVec <- c(AlphaVec, Alpha[1:(N-i+1)])
  }
  return(AlphaVec)
}

AlphaMatrix <- diag(1, nrow = N, ncol = N)
for(i in 2:N){
  Aux <- cbind(diag(1, nrow = (N-i+1), ncol = (N-i+1)), matrix(0, nrow = (N-i+1), ncol = i-1))
  AlphaMatrix <- rbind(AlphaMatrix, Aux)
}


RowBeta <- function(N){
  BetaRow <- NULL
  for(i in 1:N){
    BetaRow <- c(BetaRow,seq(1,(N-i+1)))
  }
  return(BetaRow)
}

################################################################################
############################Conditional Posteriors##############################
################################################################################

updateMu <- function(m,Sig2Mu,Y,Alpha,Beta,Sig2,Lambda,N){
  V.Alpha <- VecAlpha(Alpha,N)
  
  sig2.post <- 1/(sum(Lambda/Sig2) + 1/Sig2Mu)
  m.post <- sig2.post * (m/Sig2Mu + sum((Lambda/Sig2)*(Y-V.Alpha-Beta)))
  return(rnorm(1,m.post,sqrt(sig2.post)))
}

updateSig2 <- function(a,b,Y,Alpha,Beta,Mu,Lambda,N){
	V.Alpha <- VecAlpha(Alpha,N)

	a.post <- (N*(N+1)/4) + a
	sumsq <- sum(Lambda*((Y-Mu-V.Alpha-Beta)^2))/2
	b.post <- sumsq + b
	return(1/rgamma(1,a.post,b.post))
}

updateAlpha <- function(Alpha,Sig2Al,Y,Beta,Mu,Sig2,Lambda,AlphaMatrix,N){
	alpha <- array(0,dim=N)
	sumdata <- (((Lambda/Sig2)*(Y - Mu - Beta))%*%AlphaMatrix)
	auxdata <- (Lambda/Sig2)%*%AlphaMatrix

	for(i in 1:N){
		if(i==1){
			sig2.post <- 1/(auxdata[i] + 2/Sig2Al)
			m.post <- (Alpha[i+1]/Sig2Al + sumdata[i]) * sig2.post
			alpha[i] <- rnorm(1,m.post,sqrt(sig2.post))
		} else if(i==N){
			sig2.post <- 1/(auxdata[i] + 1/Sig2Al)
			m.post <- (Alpha[i-1]/Sig2Al + sumdata[i]) * sig2.post
			alpha[i] <- rnorm(1,m.post,sqrt(sig2.post))
		} else {
			sig2.post <- 1/(auxdata[i] + 2/Sig2Al)
			m.post <- ((Alpha[i-1] + Alpha[i+1])/Sig2Al + sumdata[i]) * sig2.post
			alpha[i] <- rnorm(1,m.post,sqrt(sig2.post))
		}
	}
	
	alpha <- c(0,alpha[-1])
	return(alpha)
}


updateBeta <- function(Beta,Sig2Be,Y,Alpha,Mu,Sig2,Lambda,N){
	V.Alpha <- VecAlpha(Alpha,N)
	Index <- RowBeta(N)
	beta <- array(0,dim=length(N*(N+1)/2))
	auxdata <- (Lambda/Sig2)*(Y-Mu-V.Alpha)
	aux <- Lambda/Sig2

	for(i in 1:(N*(N+1)/2)){
		if(Index[i]==1){
			sig2.post <- 1/(aux[i] + 2/Sig2Be)
			if(is.na(Beta[i+1])==F){
				m.post <- (Beta[i+1]/Sig2Be + auxdata[i]) * sig2.post
			} else {
				m.post <- auxdata[i]*sig2.post
			}
			beta[i] <- rnorm(1,m.post,sqrt(sig2.post))
		} else if(Index[i+1]==1){
			sig2.post <- 1/(aux[i] + 1/Sig2Be)
			m.post <- (Beta[i-1]/Sig2Be + auxdata[i]) * sig2.post
			beta[i] <- rnorm(1,m.post,sqrt(sig2.post))
		} else {
			sig2.post <- 1/(aux[i] + 2/Sig2Be)
			m.post <- ((Beta[i-1] + Beta[i+1])/Sig2Be + auxdata[i]) * sig2.post
			beta[i] <- rnorm(1,m.post,sqrt(sig2.post))
		}
	}
	beta <- c(rep(0,N),beta[-c(1:N)])
	return(beta)
}


updateSig2Be <- function(Beta,N){
	Index <- RowBeta(N)
	a.post <- (N*(N+1)/4) + 1
	sumsq <- array(0,dim=(N*(N+1)/2))
	
	for(i in 1:(N*(N+1)/2)){
		if(Index[i]==1){
			sumsq[i] <- Beta[i]^2
		} else {
			sumsq[i] <- (Beta[i] - Beta[i-1])^2
		}
	}	
	b.post <- sum(sumsq)/2
	sig2be <- 1/rgamma(1,a.post,b.post)
	return(sig2be)
}


updateLambda <- function(Nu,Y,Alpha,Beta,Mu,Sig2,N){
	V.Alpha <- VecAlpha(Alpha,N)

	c.post <- rep((Nu + 1)/2, length(Y))
	d.post <- ((Y-Mu-V.Alpha-Beta)^2)/(2*Sig2) + Nu/2
	return(rgamma(length(Y),c.post,d.post))
}

postNu <- function(Nu,Lambda,N){
  n <- N*(N+1)/2
  prior <- 0.5*(Nu/(Nu+3)) + 0.5*log(trigamma(Nu/2) - trigamma((Nu+1)/2)
						- 2*(Nu+3)/(Nu*(Nu+1)^2))
  like <- (Nu/2)*sum(log(Lambda)) - Nu/2 * sum(Lambda) + n*Nu/2*log(Nu/2) -
    n*log(gamma(Nu/2))
  return(prior+like)
}

updateNu <- function(Nu,Lam,Sm,Sigma,N,t,NuMean){#Nu is the starting value
  accept <- NULL
  count <- 0
  c <- .8
  ggamma <- 1/t^c
  
  I <- 2
  while(I>0){
    proposal <- Nu + sqrt(Sm*Sigma)*rnorm(1,0,1)
    if(proposal>0.1 & proposal<40){
      I <- 0
    }
  }
  
  probab <- min(1,exp(postNu(proposal,Lam,N)-postNu(Nu,Lam,N)))
  
  if(runif(1)<probab){
    accept <- proposal
    count <- 1
  } else{
    accept <- Nu
  }
  
  lSm <- log(Sm)+ggamma*(probab-.234)
  Sm <- exp(lSm)	
  
  Sigma <- Sigma+ggamma*(((Nu-NuMean)^2)-Sigma)
  NuMean <- NuMean+ggamma*(Nu-NuMean)
  
  return(list(accept,Sm,Sigma,NuMean,count))
}


################################################################################
################################Gibbs Algorithm#################################
################################################################################

Niter <- 200000
N <- 9

Mu.out <- array(NA,dim=Niter)
Sig2.out <- array(NA,dim=Niter)
Alpha.out <- array(NA,dim=c(Niter,N))
Beta.out <- array(NA, dim = c(Niter, (N*(N+1)/2)))
Sig2Be.out <- array(NA, dim = Niter)
Lambda.out <- array(NA, dim = c(Niter, (N*(N+1)/2)))
Nu.out <- array(NA, dim = Niter)
Sm <- array(NA, dim = Niter)
Sigma <- array(NA, dim = Niter)
NuMean <- array(NA, dim = Niter)
Count <- array(0,dim=Niter)

#Initial Values
Mu.out[1] <- 8
Sig2.out[1] <- 0.3
Alpha.out[1,] <- rep(1,N)
Beta.out[1,] <- length(c(rep(0,N),rep(1,(N*(N+1)/2)-N)))
Sig2Be.out[1] <- 1
Lambda.out[1,] <- rep(1,(N*(N+1)/2))
Nu.out[1] <- 2
Sm[1] <- 2.4^2
Sigma[1] <- 1
NuMean[1] <- 2


t1 <- Sys.time()
for(i in 2:Niter){
	NuUpdate <- updateNu(Nu.out[i-1],Lambda.out[i-1,],Sm[i-1],Sigma[i-1],N,i,NuMean[i-1])
	Nu.out[i] <- NuUpdate[[1]]
	Sm[i] <- NuUpdate[[2]]
	Sigma[i] <- NuUpdate[[3]]
	NuMean[i] <- NuUpdate[[4]]
	Count[i] <- NuUpdate[[5]]

	Mu.out[i] <- updateMu(0,1000,Z,Alpha.out[i-1,],Beta.out[i-1,],Sig2.out[i-1],Lambda.out[i-1,],N)
	Sig2.out[i] <- updateSig2(0.01,0.01,Z,Alpha.out[i-1,],Beta.out[i-1,],Mu.out[i],Lambda.out[i-1,],N)
	Alpha.out[i,] <- updateAlpha(Alpha.out[i-1,],100,Z,Beta.out[i-1,],Mu.out[i],Sig2.out[i],Lambda.out[i-1,],AlphaMatrix,N)
	Beta.out[i,] <- updateBeta(Beta.out[i-1,],Sig2Be.out[i-1],Z,Alpha.out[i,],Mu.out[i],Sig2.out[i],Lambda.out[i-1,],N)
	Sig2Be.out[i] <- updateSig2Be(Beta.out[i,],N)
	Lambda.out[i,] <- updateLambda(Nu.out[i],Z,Alpha.out[i,],Beta.out[i,],Mu.out[i],Sig2.out[i],N)
	print(i)
}
(Sys.time() - t1)

Nburn <- 100000
Thinning <- 10

#Inspecting convergence of some of the parameters
hist(Mu.out[seq((Nburn+1),Niter,by=Thinning)])
mean(Mu.out[seq((Nburn+1),Niter,by=Thinning)])
plot(Mu.out[seq((Nburn+1),Niter,by=Thinning)],type='l')

hist(Sig2.out[seq((Nburn+1),Niter,by=Thinning)])
mean(Sig2.out[seq((Nburn+1),Niter,by=Thinning)])
plot(Sig2.out[seq((Nburn+1),Niter,by=Thinning)],type='l')

hist(Nu.out[seq((Nburn+1),Niter,by=Thinning)])
mean(Nu.out[seq((Nburn+1),Niter,by=Thinning)])
sum(Count)/Niter
plot(Nu.out[seq((Nburn+1),Niter,by=Thinning)],type='l')

hist(Alpha.out[seq((Nburn+1),Niter,by=Thinning),4])
mean(Alpha.out[seq((Nburn+1),Niter,by=Thinning),4])
plot(Alpha.out[seq((Nburn+1),Niter,by=Thinning),4],type='l')

hist(Beta.out[seq((Nburn+1),Niter,by=Thinning),12])
mean(Beta.out[seq((Nburn+1),Niter,by=Thinning),12])
plot(Beta.out[seq((Nburn+1),Niter,by=Thinning),126],type='l')


hist(Sig2Be.out[seq((Nburn+1),Niter,by=Thinning)])
sigbe <- sqrt(Sig2Be.out)
mean(sigbe)

################################################################################
#################################Burn in Phase##################################
################################################################################

AlphaNames <- NULL
BetaNames <- NULL
LambdaNames <- NULL
for(i in 1:N){
	AlphaNames <- c(AlphaNames, paste("alpha", i, sep = ""))
}
for(i in 1:(N*(N+1)/2)){
	BetaNames <- c(BetaNames, paste("beta", i, sep = ""))
	LambdaNames <- c(LambdaNames, paste("lambda", i, sep = ""))
}

colnames(Alpha.out) <- AlphaNames
colnames(Beta.out) <- BetaNames
colnames(Lambda.out) <- LambdaNames

MuSim <- Mu.out[(Nburn+1):Niter]
Sig2Sim <- Sig2.out[(Nburn+1):Niter]
AlphaSim <- Alpha.out[(Nburn+1):Niter,]
BetaSim <- Beta.out[(Nburn+1):Niter,]
Sig2BeSim <- Sig2Be.out[(Nburn+1):Niter]
NuSim <- Nu.out[(Nburn+1):Niter]
LambdaSim <- Lambda.out[(Nburn+1):Niter,]

#In case it is necessary to exclude the MCMC samples after burn-in
#rm(Mu.out);rm(Sig2.out);rm(Alpha.out);rm(Beta.out);rm(Sig2Be.out);rm(Lambda.out);rm(Nu.out)

################################################################################
###############################Forecasting Claims###############################
################################################################################
#Matrix where the forecasts will be stored
Claimfcast <- array(NA, dim = c(dim(MuSim),(N*N-(N*(N+1)/2))))

row1 <- NULL
col1 <- NULL
k <- 1
N <- 9
h <- 9

for(i in 1:(N-1)){
	row1 <- c(row1,seq(N,(i+1)))
}

col1[1] <- 17
for(i in 2:(N-1)){
	col1[i] <- col1[i-1]+N-i
}

#This function samples betas outside of the upper triangle, once they are
#needed to do the forecasts
Betafcast <- function(Beta,Sig2Be,H){
  BetaF <- array(NA,dim=c(dim(MuSim),H))
  
  for(i in 1:dim(MuSim)){
	if(is.matrix(Beta)==FALSE){
		BetaF[i] <- rnorm(1,Beta[i],sqrt(Sig2Be[i]))
	}
    	for(j in 1:H){
		if(is.matrix(Beta)==FALSE){break}
      	BetaF[i,j] <- rnorm(1,Beta[i,j],sqrt(Sig2Be[i]))
    	}
  }
  	return(BetaF)
}

BetaArray <- NULL
BetaSim2 <- BetaSim[,col1]

t11 <- Sys.time()
for(i in 1:(N-1)){
	H <- N - i
	BetaF <- Betafcast(BetaSim2,Sig2BeSim,H) #Sampling new betas
	BetaArray <- cbind(BetaArray,BetaF)
	if(dim(BetaF)[1] == 1){
		BetaSim2 <- BetaF
	} else {
		BetaSim2 <- BetaF[,-c(1)]
	}	
}
(Sys.time()-t11)

#Forecasting
t3 <- Sys.time()
for(i in 1:dim(MuSim)){
	for(j in 1:(N*N-(N*(N+1)/2))){
		Lambda <- rgamma(1,NuSim[i]/2,NuSim[i]/2)
		Claimfcast[i,j] <- exp(MuSim[i] + AlphaSim[i,row1[j]] + BetaArray[i,j] + sqrt(Sig2Sim[i]/Lambda)*rnorm(1))
	}
}
t4 <- Sys.time()
(time2 <- t4-t3)

#Fitting values, in case you want to calculate a model selection information
#criterion
Fitted <- array(NA, dim = c(dim(MuSim),(N*(N+1)/2)))

Rows <- NULL
Cols <- 1:(N*(N+1)/2)
for(i in 1:N){
	Rows <- c(Rows,seq(1,(N-i+1)))
}

t8 <- Sys.time()
for(i in 1:dim(MuSim)){
	for(j in 1:(N*(N+1)/2)){
		Fitted[i,j] <- exp(MuSim[i] + AlphaSim[i,Rows[j]] + BetaSim[i,Cols[j]] + rnorm(1,0,sqrt(Sig2Sim[i]/LambdaSim[i,Cols[j]])))
	}
}	
t9 <- Sys.time()
(time <- t9 - t8)


#############################################################################
#############################################################################
#############################################################################
#Comparing the forecasts with the values of the lower triangle 
RunOff <- array(NA, dim = c(N,N))
mediansfcast <- apply(Claimfcast,2,median)
mediansfit <- apply(Fitted,2,median)

col2 <- NULL

for(i in 1:(N-1)){
	col2 <- c(col2,seq(i+1,N))
}

col3 <- NULL
for(i in 1:N){
	col3 <- c(col3,rep(i,(N-i+1)))
}

for(i in 1:(N*(N+1)/2)){
	RunOff[Rows[i],col3[i]] <- mediansfit[i]
}
for(i in 1:(N*N-(N*(N+1)/2))){
	RunOff[row1[i],col2[i]] <- mediansfcast[i]
}

plot(data[1,],type='l')
lines(RunOff[1,],type='l',col='red')

plot(data[2,],type='l', ylim=c(1000,11000))
lines(RunOff[2,],type='l',col='red')
lines(test[2,],type='p',col='blue') #True value

plot(data[3,],type='l',ylim=c(1000,11000))
lines(RunOff[3,],type='l',col='red')
lines(test[3,],type='l',col='blue')

plot(data[4,],type='l',ylim=c(1000,15000))
lines(RunOff[4,],type='l',col='red')
lines(test[4,],type='l',col='blue')

plot(data[5,],type='l',ylim=c(1000,15000))
lines(RunOff[5,],type='l',col='red')
lines(test[5,],type='l',col='blue')

plot(data[6,],type='l',ylim=c(1000,15000))
lines(RunOff[6,],type='l',col='red')
lines(test[6,],type='l',col='blue')

plot(data[7,],type='l',ylim=c(1000,15000))
lines(RunOff[7,],type='l',col='red')
lines(test[7,],type='l',col='blue') #Structural break in the series that destroys the forecasts

plot(data[8,],type='l',ylim=c(1000,15000))
lines(RunOff[8,],type='l',col='red')
lines(test[8,],type='l',col='blue')

plot(data[9,],type='p',ylim=c(1000,15000))
lines(RunOff[9,],type='l',col='red')
lines(test[9,],type='l',col='blue')

##############################################################################
##########################Deviance Information Criterion######################
##############################################################################
#This is a function to calculete the DIC

Param <- cbind(MuSim,Sig2Sim,AlphaSim,BetaSim,LambdaSim)

logLike <- function(Y, Theta, N) {
	Mu <- Theta[1]
	Sig2 <- Theta[2]
	Alpha <- Theta[3:11]
	Beta <- Theta[12:56]
	Lambda <- Theta[57:101]

  	Mean <- Mu + VecAlpha(Alpha, N) + Beta
	sum( dnorm(Y, Mean, sqrt(Sig2/Lambda), log = TRUE) )
}

calculateDIC <- function(Y, Theta, N) {
  Theta_hat = apply(Theta, 2, mean)
  L = logLike(Y, Theta_hat, N)
  
  #Calculate P
  S = nrow(Theta) #S = number of iterations
  #Add up the log likelihoods of each iteration
  llSum = 0
  for (s in 1:S) {
    Theta_s = Theta[s,]
    llSum = llSum + logLike(Y, Theta_s, N)
  }
  Dbar <- -2 * (1 / S * llSum)
  P = 2 * L + Dbar 
  
  DIC = Dbar + P
  
  return(list(DIC=DIC, P=P, LogLik=L))
}

t8 <- Sys.time()
calculateDIC(Z, Param, N) 

t9 <- Sys.time()
(time <- t9 - t8)

#$DIC
#[1] 26.17753

#$P
#[1] 24.87375

#$LogLik
#[1] 11.78499