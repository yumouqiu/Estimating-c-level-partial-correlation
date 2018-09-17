resultprop1 = read.table("BD30360100resultpropC0.txt", head = FALSE)
resultLiu1 = read.table("BD30360100resultLiuC0.txt", head = FALSE)
resultprop2 = read.table("BD30380100resultpropC0.txt", head = FALSE)
resultLiu2 = read.table("BD30380100resultLiuC0.txt", head = FALSE)

n1 = 60; n2 = 80; p = 100; c0 = 0.25; tau = seq(0.1, 3.1, 0.01); kappa = 1
STrueNum = resultprop1[1, 1]
#STrueNum = 25 * (p / 5) - p
M = dim(resultprop1)

IndFP = 1 + (1 : (M[2] / 4) - 1) * 4 + 1
IndFN = 1 + (1 : (M[2] / 4) - 1) * 4 + 2
IndTP = 1 + (1 : (M[2] / 4) - 1) * 4 + 3
IndTN = 1 + (1 : (M[2] / 4)) * 4

alpha = c(0.05, 0.1, 0.15, 0.2)
lenalpha = length(alpha)
FDPresultpropC01 = matrix(0, 1000, 3 * lenalpha); FDPresultLiuC01 = matrix(0, 1000, 3 * lenalpha)
FDPresultpropC02 = matrix(0, 1000, 3 * lenalpha); FDPresultLiuC02 = matrix(0, 1000, 3 * lenalpha)

for (rep in 1 : M[1]){
	prop1 = resultprop1[rep, ]
	Liu1 = resultLiu1[rep, ]
	prop2 = resultprop2[rep, ]
	Liu2 = resultLiu2[rep, ]

	FPprop1 = as.numeric(prop1[IndFP]); FNprop1 = as.numeric(prop1[IndFN]); TPprop1 = as.numeric(prop1[IndTP]); TNprop1 = as.numeric(prop1[IndTN])
	FPLiu1 = as.numeric(Liu1[IndFP]); FNLiu1 = as.numeric(Liu1[IndFN]); TPLiu1 = as.numeric(Liu1[IndTP]); TNLiu1 = as.numeric(Liu1[IndTN])
	FPprop2 = as.numeric(prop2[IndFP]); FNprop2 = as.numeric(prop2[IndFN]); TPprop2 = as.numeric(prop2[IndTP]); TNprop2 = as.numeric(prop2[IndTN])
	FPLiu2 = as.numeric(Liu2[IndFP]); FNLiu2 = as.numeric(Liu2[IndFN]); TPLiu2 = as.numeric(Liu2[IndTP]); TNLiu2 = as.numeric(Liu2[IndTN])

	rejectpropC01 = pmax(FPprop1 + TPprop1, rep(1, length(FPprop1))); rejectLiuC01 = pmax(FPLiu1 + TPLiu1, rep(1, length(FPLiu1)))
	rejectpropC02 = pmax(FPprop2 + TPprop2, rep(1, length(FPprop2))); rejectLiuC02 = pmax(FPLiu2 + TPLiu2, rep(1, length(FPLiu2)))

	FDPpropC01 = 2 * ( sqrt(n1) * p * ( 1 - pnorm( tau * sqrt(log(p)) ) ) + p * (p - 1 - sqrt(n1)) * ( 1 - pnorm( sqrt(n1) * c0 / sqrt(kappa) + tau * sqrt(log(p)) ) ) ) / rejectpropC01
	FDPLiuC01 = 2 * p * (p - 1) * ( 1 - pnorm( c0 * sqrt(n1) + tau * sqrt(log(p)) ) ) / rejectLiuC01
	FDPpropC02 = 2 * ( sqrt(n2) * p * ( 1 - pnorm( tau * sqrt(log(p)) ) ) + p * (p - 1 - sqrt(n2)) * ( 1 - pnorm( sqrt(n2) * c0 / sqrt(kappa) + tau * sqrt(log(p)) ) ) ) / rejectpropC02
	FDPLiuC02 = 2 * p * (p - 1) * ( 1 - pnorm( c0 * sqrt(n2) + tau * sqrt(log(p)) ) ) / rejectLiuC02

	FDPrespropC01 = c(); FDPresLiuC01 = c()
	FDPrespropC02 = c(); FDPresLiuC02 = c()
	for (i in 1 : lenalpha){
		if (sum(FDPpropC01 <= alpha[i]) > 0) {taupropC01 = min(c(2, tau[FDPpropC01 <= alpha[i]])); kprop1 = which(tau == taupropC01)}
		else {taupropC01 = 2; kprop1 = 191}
		powerpropC01 = as.numeric(TPprop1[kprop1] / STrueNum)
		FDPalphapropC01 = FPprop1[kprop1] / max(FPprop1[kprop1] + TPprop1[kprop1], 1)
		tempFDPpropC01 = c(taupropC01, FDPalphapropC01, powerpropC01)
		FDPrespropC01 = c(FDPrespropC01, tempFDPpropC01)

		if (sum(FDPpropC02 <= alpha[i]) > 0) {taupropC02 = min(c(2, tau[FDPpropC02 <= alpha[i]])); kprop2 = which(tau == taupropC02)}
		else {taupropC02 = 2; kprop2 = 191}
		powerpropC02 = as.numeric(TPprop2[kprop2] / STrueNum)
		FDPalphapropC02 = FPprop2[kprop2] / max(FPprop2[kprop2] + TPprop2[kprop2], 1)
		tempFDPpropC02 = c(taupropC02, FDPalphapropC02, powerpropC02)
		FDPrespropC02 = c(FDPrespropC02, tempFDPpropC02)

		if (sum(FDPLiuC01 <= alpha[i]) > 0) {tauLiuC01 = min(c(2, tau[FDPLiuC01 <= alpha[i]])); kLiu1 = which(tau == tauLiuC01)}
		else {tauLiuC01 = 2; kLiu1 = 191}
		powerLiuC01 = as.numeric(TPLiu1[kLiu1] / STrueNum)
		FDPalphaLiuC01 = FPLiu1[kLiu1] / max(FPLiu1[kLiu1] + TPLiu1[kLiu1], 1)
		tempFDPLiuC01 = c(tauLiuC01, FDPalphaLiuC01, powerLiuC01)
		FDPresLiuC01 = c(FDPresLiuC01, tempFDPLiuC01)

		if (sum(FDPLiuC02 <= alpha[i]) > 0) {tauLiuC02 = min(c(2, tau[FDPLiuC02 <= alpha[i]])); kLiu2 = which(tau == tauLiuC02)}
		else {tauLiuC02 = 2; kLiu2 = 191}
		powerLiuC02 = as.numeric(TPLiu2[kLiu2] / STrueNum)
		FDPalphaLiuC02 = FPLiu2[kLiu2] / max(FPLiu2[kLiu2] + TPLiu2[kLiu2], 1)
		tempFDPLiuC02 = c(tauLiuC02, FDPalphaLiuC02, powerLiuC02)
		FDPresLiuC02 = c(FDPresLiuC02, tempFDPLiuC02)
	}

	FDPresultpropC01[rep, ] = FDPrespropC01; FDPresultpropC02[rep, ] = FDPrespropC02
	FDPresultLiuC01[rep, ] = FDPresLiuC01; FDPresultLiuC02[rep, ] = FDPresLiuC02
}

MFDP = dim(FDPresultpropC01)
FDPprop1 = colMeans(FDPresultpropC01)
FDPLiu1 = colMeans(FDPresultLiuC01)
FDPprop2 = colMeans(FDPresultpropC02)
FDPLiu2 = colMeans(FDPresultLiuC02)

SDprop1 = c(); SDprop2 = c()
SDLiu1 = c(); SDLiu2 = c()
for (i in 1 : MFDP[2]){
	SDprop1[i] = sd(FDPresultpropC01[, i])
	SDLiu1[i] = sd(FDPresultLiuC01[, i])
	SDprop2[i] = sd(FDPresultpropC02[, i])
	SDLiu2[i] = sd(FDPresultLiuC02[, i])
}

FDPIndtau = (1 : (MFDP[2] / 3) - 1) * 3 + 1
FDPInd = (1 : (MFDP[2] / 3) - 1) * 3 + 2
FDPIndPower = (1 : (MFDP[2] / 3)) * 3

#cbind(FDPprop1[FDPInd], FDPLiu1[FDPInd], FDPprop1[FDPIndPower], FDPLiu1[FDPIndPower])
#cbind(FDPprop2[FDPInd], FDPLiu2[FDPInd], FDPprop2[FDPIndPower], FDPLiu2[FDPIndPower])



FDPresultprop1 = read.table("BD30360100FDPresultpropC0.txt", head = FALSE)
FDPresultLiu1 = read.table("BD30360100FDPresultLiuC0.txt", head = FALSE)
FDPresultprop2 = read.table("BD30380100FDPresultpropC0.txt", head = FALSE)
FDPresultLiu2 = read.table("BD30380100FDPresultLiuC0.txt", head = FALSE)

FDPprop11 = colMeans(FDPresultprop1)
FDPLiu11 = colMeans(FDPresultLiu1)
FDPprop12 = colMeans(FDPresultprop2)
FDPLiu12 = colMeans(FDPresultLiu2)


c(FDPprop1[FDPInd[2]], FDPprop1[FDPIndPower[2]], FDPprop11[FDPInd[2]], FDPprop11[FDPIndPower[2]])
c(FDPprop2[FDPInd[2]], FDPprop2[FDPIndPower[2]], FDPprop12[FDPInd[2]], FDPprop12[FDPIndPower[2]])


r1alpha01 = c(FDPprop1[FDPInd[2]], FDPprop1[FDPIndPower[2]], FDPLiu11[FDPInd[2]], FDPLiu11[FDPIndPower[2]], FDPLiu1[FDPInd[2]], FDPLiu1[FDPIndPower[2]])
r2alpha01 = c(FDPprop2[FDPInd[2]], FDPprop2[FDPIndPower[2]], FDPLiu12[FDPInd[2]], FDPLiu12[FDPIndPower[2]], FDPLiu2[FDPInd[2]], FDPLiu2[FDPIndPower[2]])

round(rbind(r1alpha01, r2alpha01), 3)