############################################################
# choose the causal estimand
############################################################

lW <- c("ATO", "ATT", "ATE") # labels of causal estimand
nW <- length(lW)

lN <- paste0(lW, "N") # labels of weighting Hàjek (nonaugmented) estimators
lA <- paste0(lW, "A") # labels of weighting augmented estimators
lE <- c(lN, lA) # labels of estimators
nE <- length(lE)

cE <- "ATOA" # chosen estimator


############################################################
# import the data
############################################################

# data should be imported into a matrix D, with variables in columns
# ordinal running variable should be coded as a number


############################################################
# extract individual variables from the data matrix
############################################################

tV <- 0 # threshold value of the running variable
tC <- "-1|0" # index of the cutoff vector for to the threshold value
cF <- 1 # number of column containing the running variable
cO <- 2 # number of column containing the outcome variable
cC <- c(3,4,5,6) # numbers of columns containing the pretreatment variables
cP <- c(3,4) # numbers of columns containing the predictors of treatment
cR <- c(3,4,5) # numbers of columns containing predictors of outcome
Z <- ifelse(D[,cF] < tV | is.na(D[,cF]), 0, 1) # treatment indicator
R <- as.ordered(D[,cF]) # runinng variable as an ordered factor
Y <- D[,cO] # outcome variable
X <- D[,cC] # pretreatment covariates
K <- ncol(X)


####################################################################################################
# specify and estimate the ordered probit model for the distribution of the ordinal running variable
####################################################################################################

library(MASS)

P <- D[,cP] # matrix of predictors of treatment
nP <- ncol(P)

# matrix P should be sensibly scaled to avoid having to specify starting values
oprobit <- polr(R ~ P, Hess = T, subset = complete.cases(X), method = "probit")

mu <- oprobit$zeta # vector of cutoffs
beta <- oprobit$coefficients # vector of coefficients
C <- mu[tC] - P %*% beta # linear predictor

eH <- 1 - pnorm(C) # estimated propensity scores

w <- matrix(0, nrow = length(eH), ncol = nW)
colnames(w) <- lW

w[,"ATO"] <- ifelse(Z == 0, eH, 1 - eH) # weights for ATO
w[,"ATT"] <- ifelse(Z == 0, eH / ( 1 - eH ), 1) # weights for ATT
w[,"ATE"] <- ifelse(Z == 0, ( 1 - eH ) ^ (-1), eH ^ (-1)) # weights for ATE


############################################################
# specify and estimate the outcome model
############################################################

Q <- D[,cR] # matrix of predictors of outcome

lM0 <- lm(Y ~ Q, subset = complete.cases(X) & ( Z == 0 )) # model for controls
lM1 <- lm(Y ~ Q, subset = complete.cases(X) & ( Z == 1 )) # model for treated

gamma0 <- lM0$coefficients # vector of coefficients for controls
gamma1 <- lM1$coefficients # vector of coefficient for treated

M0 <- cbind(1,Q) %*% gamma0 # predicted values for controls
M1 <- cbind(1,Q) %*% gamma1 # predicted values for treated


############################################################
# specify the subsample and estimate the treatment effect
############################################################

h <- c(0.34, 0.35, 0.36, 0.37, 0.38) # symmetric intervals
#h <- c(0.01, 0.02, 0.03, 0.04, 0.05) # asymmetric intervals
H <- length(h)

WATE <- matrix(0, nrow = H, ncol = nE) # weighted average treatment effects
colnames(WATE) <- lE
stdE <- matrix(0, nrow = H, ncol = nE) # vector of standard errors
colnames(stdE) <- lE
pVal <- matrix(0, nrow = H, ncol = nE) # vector of p-values
colnames(pVal) <- lE

nObs <- rep(0, H) # vector of numbers of observations
nC <- rep(0, H) # vector of numbers of controls
nT <- rep(0, H) # vector of numbers of treated
ASD <- array(0, dim = c(H, K, nW), dimnames = list(NULL, NULL, lW)) # matrix of absolute standardized differences

for (i in 1:H) {

iSS <- ( abs(eH - 0.5) < h[i] ) & complete.cases(X, Y, R) # index vector for symmetric interval
#iSS <- ( eH < 0.88 ) & ( eH > (0.12 - h[i]) ) & complete.cases(cbind(X, Y, R)) # index vector for asymmetric interval

nObs[i] <- sum(iSS) # number of observations in the subsample
nC[i] <- sum(Z[iSS]==0) # number of controls in the subsample
nT[i] <- sum(Z[iSS]==1) # number of treated in the subsample

wlsATO <- lm(Y ~ Z, subset = iSS, weights = w[,"ATO"])
WATE[i,"ATON"] <- unname(wlsATO$coefficients[2]) # ATO (nonaugmented estimator)

wlsATT <- lm(Y ~ Z, subset = iSS, weights = w[,"ATT"])
WATE[i,"ATTN"] <- unname(wlsATT$coefficients[2]) # ATT (nonaugmented estimator)

wlsATE <- lm(Y ~ Z, subset = iSS, weights = w[,"ATE"])
WATE[i,"ATEN"] <- unname(wlsATE$coefficients[2]) # ATE (nonagumented estimator)

WATE[i,"ATOA"] <- sum(( ( 1 - eH ) * ( Z * Y - ( Z - eH ) * M1 ) )[iSS]) /
sum((eH * ( 1 - eH ) )[iSS]) -
sum((eH * ( ( 1 - Z ) * Y + ( Z - eH ) * M0 ) )[iSS]) /
sum(( ( 1 - eH ) * eH )[iSS]) # ATO (augmented estimator)
 
WATE[i,"ATTA"] <- weighted.mean(Y[( Z == 1 ) & iSS], w[( Z == 1 ) & iSS,"ATT"]) -
( sum(( Y - M0 )[( Z == 0 ) & iSS] * w[( Z == 0 ) & iSS,"ATT"]) + 
sum(M0[( Z == 1 ) & iSS] * w[( Z == 1 ) & iSS,"ATT"]) ) /
sum(w[( Z == 1 ) & iSS,"ATT"]) # ATT (augmented estimator)

WATE[i,"ATEA"] <- ( sum(( Y - M1 )[( Z == 1 ) & iSS] * w[( Z == 1 ) & iSS,"ATE"]) -
sum(( Y - M0 )[( Z == 0 ) & iSS] * w[( Z == 0 ) & iSS,"ATE"]) ) / nObs[i] +
( sum(M1[iSS]) - sum(M0[iSS]) ) / nObs[i] # ATE (augmented estimator)


############################################################
# compute the balancing statistics
############################################################

for (j in 1:K) {
	for (l in lW) {

	WM1 <- weighted.mean(X[,j][( Z == 1 ) & iSS], w[( Z == 1 ) & iSS,l], na.rm = T)
	WM0 <- weighted.mean(X[,j][( Z == 0 ) & iSS], w[( Z == 0 ) & iSS,l], na.rm = T)
	V1 <- sd(X[,j][( Z == 1 ) & iSS], na.rm = T) ^ 2 / sum(!is.na(Y[( Z == 1 ) & iSS]))
	V0 <- sd(X[,j][( Z == 0 ) & iSS], na.rm = T) ^ 2 / sum(!is.na(Y[( Z == 0 ) & iSS]))
	ASD[i,j,l] <- ( WM1 - WM0 ) / sqrt(V1 + V0) # absolute standardized difference

	}
}


############################################################
# calculate standard errors and p-values
############################################################

# components of the treatment effect estimate
T0 <- sum(( w[,"ATO"] * Y )[( Z == 0 ) & iSS]) / sum(w[( Z == 0 ) & iSS,"ATO"]) # ATO (nonaugmented estimator)
T1 <- sum(( w[,"ATO"] * Y )[( Z == 1 ) & iSS]) / sum(w[( Z == 1 ) & iSS,"ATO"]) # ATO (nonaugmented estimator)

S0 <- sum(( w[,"ATT"] * Y )[( Z == 0 ) & iSS]) / sum(w[( Z == 0 ) & iSS,"ATT"]) # ATT (nonaugmented estimator)
S1 <- sum(( w[,"ATT"] * Y )[( Z == 1 ) & iSS]) / sum(w[( Z == 1 ) & iSS,"ATT"]) # ATT (nonaugmented estimator)

R0 <- sum(( w[,"ATE"] * Y )[( Z == 0 ) & iSS]) / sum(w[( Z == 0 ) & iSS,"ATE"]) # ATE (nonaugmented estimator)
R1 <- sum(( w[,"ATE"] * Y )[( Z == 1 ) & iSS]) / sum(w[( Z == 1 ) & iSS,"ATE"]) # ATE (nonaugmented estimator)

W0 <- sum(( eH * ( ( 1 - Z ) * Y + ( Z - eH ) * M0 ) )[iSS]) / sum(( eH * ( 1 - eH ) )[iSS]) # ATO (augmented estimator)
W1 <- sum(( ( 1 - eH ) * ( Z * Y - ( Z - eH ) * M1 ) )[iSS]) / sum(( eH * ( 1 - eH ) )[iSS]) # ATO (augmented estimator)

V1 <- sum(( w[,"ATT"] * Y )[( Z == 1 ) & iSS]) / sum(w[( Z == 1 ) & iSS,"ATT"]) # ATT (augmented estimator)
V0 <- ( sum(( w[,"ATT"] * Y )[( Z == 0 ) & iSS]) +
	    sum(( w[,"ATT"] * M0 )[( Z == 1 ) & iSS]) -
	    sum(( w[,"ATT"] * M0 )[( Z == 0 ) & iSS]) ) / sum(w[( Z == 1 ) & iSS,"ATT"]) # ATT (augmented estimator)

U1 <- sum(( w[,"ATE"] * ( Y - M1 ) )[( Z == 1 ) & iSS]) / nObs[i] + sum(M1[iSS]) / nObs[i] # ATE (augmented estimator)
U0 <- sum(( w[,"ATE"] * ( Y - M0 ) )[( Z == 0 ) & iSS]) / nObs[i] + sum(M0[iSS]) / nObs[i] # ATE (augmented estimator)

# create indicator matrix of ratings
options(na.action = "na.pass")
W <- model.matrix(~ R-1)
options(na.action = "na.omit")

u <- c(-Inf, mu, Inf) # vector of cutoffs in the ordered probit 
J <- length(u) - 1
sE <- matrix(0, nrow = nObs[i], ncol = J-1+nP) # score in the ordered probit
eX <- matrix(0, nrow = nObs[i], ncol = J-1+nP) # gradient of the propensity score
T <- which(names(mu)==tC) # indicator of the cutoff corresponding to the threshold

# first elements of sE
lB <- matrix(0, nrow = nObs[i], ncol = J)
for (j in 1:J) {
	lB[,j] <- W[iSS,j] * ( dnorm(u[j+1] - P[iSS,] %*% beta) - dnorm(u[j] - P[iSS,] %*% beta) ) /
			 ( pnorm(u[j+1] - P[iSS,] %*% beta) - pnorm(u[j] - P[iSS,] %*% beta) )
}

for (k in 1:nP) {
	sE[,k] <- rowSums(-lB) * P[iSS,k]
	eX[,k] <- -P[iSS,k] * dnorm(u[T+1] - P[iSS,] %*% beta)
}

for (j in 1:(J-1)) {
	sE[,nP+j] <- W[iSS,j] * dnorm(u[j+1] - P[iSS,] %*% beta) /
				( pnorm(u[j+1] - P[iSS,] %*% beta) - pnorm(u[j] - P[iSS,] %*% beta) ) -
				W[iSS,j+1] * dnorm(u[j+1] - P[iSS,] %*% beta) /
				( pnorm(u[j+2] - P[iSS,] %*% beta) - pnorm(u[j+1] - P[iSS,] %*% beta) )
}

eX[,nP+T] <- dnorm(u[T+1] - P[iSS,] %*% beta)

# ATO (nonaugmented estimator)
hE <- 1 / nObs[i] * ( Z[iSS] * ( Y[iSS] - T1 ) + ( 1 - Z[iSS] ) * ( Y[iSS] - T0 ) ) %*% eX
iV <- Z[iSS] * ( Y[iSS] - T1 ) * ( 1 - eH[iSS] )  - ( 1 - Z[iSS] ) * ( Y[iSS] - T0 ) * eH[iSS] - hE %*% ( nObs[i] * vcov(oprobit) ) %*% t(sE)
tH <- 1 / nObs[i] * sum(eH[iSS] * ( 1 - eH[iSS] ))
stdE[i,"ATON"] <- sqrt(sum(iV^2) / ( nObs[i] * tH ) ^ 2 )
pVal[i,"ATON"] <- 2 * pnorm(-abs(WATE[i,"ATON"]) / stdE[i,"ATON"])

# ATT (nonaugmented estimator)
hE <- 1 / nObs[i] * ( ( 1 - Z[iSS] ) * ( Y[iSS] - S0 ) / ( 1 - eH[iSS] ) ^ 2 ) %*% eX
iV <- Z[iSS] * ( Y[iSS] - S1 )  - ( 1 - Z[iSS] ) * ( Y[iSS] - S0 ) * eH[iSS] / ( 1 - eH[iSS] ) - hE %*% ( nObs[i] * vcov(oprobit) ) %*% t(sE)
tH <- 1 / nObs[i] * sum(eH[iSS])
stdE[i,"ATTN"] <- sqrt(sum(iV^2) / ( nObs[i] * tH ) ^ 2 )
pVal[i,"ATTN"] <- 2 * pnorm(-abs(WATE[i,"ATTN"]) / stdE[i,"ATTN"])

# ATE (nonaugmented estimator)
hE <- 1 / nObs[i] * ( Z[iSS] * ( Y[iSS] - R1 ) / eH[iSS] ^ 2 + ( 1 - Z[iSS] ) * ( Y[iSS] - R0 ) / ( 1 - eH[iSS] ) ^ 2 ) %*% eX
iV <- Z[iSS] * ( Y[iSS] - R1 ) / eH[iSS]  - ( 1 - Z[iSS] ) * ( Y[iSS] - R0 ) / ( 1 - eH[iSS] ) - hE %*% ( nObs[i] * vcov(oprobit) ) %*% t(sE)
tH <- 1
stdE[i,"ATEN"] <- sqrt(sum(iV^2) / ( nObs[i] * tH ) ^ 2 )
pVal[i,"ATEN"] <- 2 * pnorm(-abs(WATE[i,"ATEN"]) / stdE[i,"ATEN"])

# ATO (augmented estimator)
hE <- - ( 1 / nObs[i] ) * ( ( M1[iSS] * ( 1 - 2 * eH[iSS] + Z[iSS] ) - M0[iSS] * ( Z[iSS] - 2 * eH[iSS] ) + ( 2 * eH[iSS] - 1 ) * ( W1 - W0 ) - Y[iSS] ) %*% eX )
hG0 <- 1 / nObs[i] * ( ( Z[iSS] - eH[iSS] ) * eH[iSS] ) %*% cbind(1,Q)[iSS,] 
sG0 <- matrix(Y[iSS] - M0[iSS], nObs[i], ncol(Q)+1) * cbind(1,Q)[iSS,]
hG1 <- 1 / nObs[i] * ( ( Z[iSS] - eH[iSS] ) * ( 1 - eH[iSS] ) ) %*% cbind(1,Q)[iSS,]
sG1 <- matrix(Y[iSS] - M1[iSS], nObs[i], ncol(Q)+1) * cbind(1,Q)[iSS,]
iV <- ( 1 - eH[iSS] ) * ( Z[iSS] * Y[iSS] - ( Z[iSS] - eH[iSS] ) * M1[iSS] - eH[iSS] * W1 ) - eH[iSS] * ( ( 1 - Z[iSS] ) * Y[iSS] + ( Z[iSS] - eH[iSS] ) * M0[iSS] - ( 1 - eH[iSS] ) *W0 )  - hE %*% ( nObs[i] * vcov(oprobit) ) %*% t(sE) - hG0 %*% ( nObs[i] * summary(lM0)$cov ) %*% t(sG0) - hG1 %*% ( nObs[i] * summary(lM1)$cov ) %*% t(sG1)
tH <- 1 / nObs[i] * sum(eH[iSS] * ( 1 - eH[iSS] ))
stdE[i,"ATOA"] <- sqrt(sum(iV^2) / ( nObs[i] * tH ) ^ 2 )
pVal[i,"ATOA"] <- 2 * pnorm(-abs(WATE[i,"ATOA"]) / stdE[i,"ATOA"])

# ATT (augmented estimator)
hE <- 1 / nObs[i] * ( ( 1 - Z[iSS] ) * ( Y[iSS] - M0[iSS] ) / ( 1 - eH[iSS] ) ^ 2 ) %*% eX
hG <- 1 / nObs[i] * ( ( Z[iSS] - eH[iSS] ) / ( 1 - eH[iSS] ) ) %*% cbind(1,Q)[iSS,]
sG <- matrix(Y[iSS] - M0[iSS], nObs[i], ncol(Q)+1) * cbind(1,Q)[iSS,]
iV <- Z[iSS] * ( Y[iSS] - V1 )  - ( 1 - Z[iSS] ) * Y[iSS] * eH[iSS] / ( 1 - eH[iSS] ) + Z[iSS] * V0 - M0[iSS] * ( Z[iSS] - eH[iSS] ) / ( 1 - eH[iSS] )  - hE %*% ( nObs[i] * vcov(oprobit) ) %*% t(sE) - hG %*% ( nObs[i] * summary(lM0)$cov ) %*% t(sG)
tH <- 1 / nObs[i] * sum(eH[iSS])
stdE[i,"ATTA"] <- sqrt(sum(iV^2) / ( nObs[i] * tH ) ^ 2 )
pVal[i,"ATTA"] <- 2 * pnorm(-abs(WATE[i,"ATTA"]) / stdE[i,"ATTA"])

# ATE (augmented estimator)
hE <- 1 / nObs[i] * ( Z[iSS] * ( Y[iSS] - M1[iSS] ) / eH[iSS] ^ 2 + ( 1 - Z[iSS] ) * ( Y[iSS] - M0[iSS] ) / ( 1 - eH[iSS] ) ^ 2 ) %*% eX
hG0 <- 1 / nObs[i] * ( ( Z[iSS] - eH[iSS] ) / ( 1 - eH[iSS] ) ) %*% cbind(1,Q)[iSS,]
sG0 <- matrix(Y[iSS] - M0[iSS], nObs[i], ncol(Q)+1) * cbind(1,Q)[iSS,]
hG1 <- 1 / nObs[i] * ( ( Z[iSS] - eH[iSS] ) / eH[iSS] ) %*% cbind(1,Q)[iSS,]
sG1 <- matrix(Y[iSS] - M1[iSS], nObs[i], ncol(Q)+1) * cbind(1,Q)[iSS,]
iV <- Z[iSS] * ( Y[iSS] - M1[iSS] ) / eH[iSS] + M1[iSS] - U1  - ( 1 - Z[iSS] ) * ( Y[iSS] - M0[iSS] ) / ( 1 - eH[iSS] ) - M0[iSS] + U0 - hE %*% ( nObs[i] * vcov(oprobit) ) %*% t(sE) - hG0 %*% ( nObs[i] * summary(lM0)$cov ) %*% t(sG0) - hG1 %*% ( nObs[i] * summary(lM1)$cov ) %*% t(sG1)
tH <- 1
stdE[i,"ATEA"] <- sqrt(sum(iV^2) / ( nObs[i] * tH ) ^ 2 )
pVal[i,"ATEA"] <- 2 * pnorm(-abs(WATE[i,"ATEA"]) / stdE[i,"ATEA"])

}


############################################################
# visualize the results
############################################################

# print the estimates and the absolute standardized differences for each subsample
print(data.frame(nC, nT, WATE[,cE], stdE[,cE], pVal[,cE], ASD[,,strtrim(cE, 3)]), digits = 3)

# boxplot showing the estimated propensity scores
boxplot(eH ~ R, subset = complete.cases(X), ylab = "estimated propensity score", frame = F, las = 3)

