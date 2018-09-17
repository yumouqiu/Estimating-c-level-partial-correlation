resultprop = read.table("BD30360100resultprop.txt", head = FALSE)
resultLiu = read.table("BD30360100resultLiu.txt", head = FALSE)
resultLSW = read.table("BD30360100resultLSW.txt", head = FALSE)
resultGL = read.table("BD30360100resultGL.txt", head = FALSE)

n = 60; p = 100; c0 = 0; tau = seq(0.1, 3.1, 0.01); kappa = 1
STrueNum = resultprop[1, 1]
#STrueNum = 25 * (p / 5) - p
M = dim(resultprop)
prop = colMeans(resultprop)
Liu = colMeans(resultLiu)
LSW = colMeans(resultLSW)
GL = colMeans(resultGL)

IndFP = 1 + (1 : (M[2] / 4) - 1) * 4 + 1
IndFN = 1 + (1 : (M[2] / 4) - 1) * 4 + 2
IndTP = 1 + (1 : (M[2] / 4) - 1) * 4 + 3
IndTN = 1 + (1 : (M[2] / 4)) * 4

FPprop = prop[IndFP]; FNprop = prop[IndFN]; TPprop = prop[IndTP]; TNprop = prop[IndTN]
FPLiu = Liu[IndFP]; FNLiu = Liu[IndFN]; TPLiu = Liu[IndTP]; TNLiu = Liu[IndTN]
FPLSW = LSW[IndFP]; FNLSW = LSW[IndFN]; TPLSW = LSW[IndTP]; TNLSW = LSW[IndTN]
FPGL = GL[IndFP]; FNGL = GL[IndFN]; TPGL = GL[IndTP]; TNGL = GL[IndTN]
rejectpropC0 = pmax(FPprop + TPprop, rep(1, length(FPprop))); rejectLiuC0 = pmax(FPLiu + TPLiu, rep(1, length(FPLiu)))

FDPpropC0 = 2 * ( sqrt(n) * p * ( 1 - pnorm( tau * sqrt(log(p)) ) ) + p * (p - 1 - sqrt(n)) * ( 1 - pnorm( sqrt(n) * c0 / sqrt(kappa) + tau * sqrt(log(p)) ) ) ) / rejectpropC0
FDPLiuC0 = 2 * p * (p - 1) * ( 1 - pnorm( tau * sqrt(log(p)) ) ) / rejectLiuC0

#plot(c(FPprop / (FPprop + TPprop), 0), c(TPprop / STrueNum, 0), xlim = c(0, 0.5), ylim = c(0, 1), type = 'l', col = 1, lty = 1, 
#xlab = "False Discovery Rate", ylab = "True Positive / Num of Signals", 
#main = expression(paste("Structure A with normal distribution,", n, "=", 80, " and ", p, "=", 100) ) )
##points(FPLiu / (FPLiu + TPLiu), TPLiu / STrueNum, col = 2, pch = 2)
#lines(FPLiu / (FPLiu + TPLiu), TPLiu / STrueNum, col = 3, lty = 3)
##points(FPLSW / (FPLSW + TPLSW), TPLSW / STrueNum, col = 3, pch = 3)


plot(c(FDPpropC0, 0), c(TPprop / STrueNum, 0), xlim = c(0, 0.5), ylim = c(0, 1), type = 'l', col = 1, lty = 1, 
xlab = "False Discovery Rate", ylab = "True Positive / Num of Signals", 
main = expression(paste("Structure A,", n, "=", 60, " and ", p, "=", 100) ) )
lines(c(FDPLiuC0, 0), c(TPLiu / STrueNum, 0), col = 2, lty = 2)
FDRGL = FPGL / (FPGL + TPGL); PowerGL = TPGL / STrueNum
FDRLSW = FPLSW / (FPLSW + TPLSW); PowerLSW = TPLSW / STrueNum
FDRLSWAdd = FDRGL[which(FDRGL[-c(181 : 201)] - FDRLSW[201] < -0.01)]
PowerLSWAdd = PowerGL[which(FDRGL[-c(181 : 201)] - FDRLSW[201] < -0.01)] + 0.02
lines(c(FDRLSW), c(PowerLSW), col = 3, lty = 3)
#points(FPGL / (FPGL + TPGL), TPGL / STrueNum, col = 4, pch = 4)
lines(FDRGL[-c(181 : 201)], PowerGL[-c(181 : 201)], col = 4, lty = 4)
legend("bottomright", legend = c("Proposed", "Liu", "LSW", "Graphic Lasso"), col = c(1, 2, 3, 4), lty = c(1, 2, 3, 4))


#plot(FPprop, TPprop, type = 'b', col = 1, lty = 1)
#lines(FPLiu, TPLiu, col = 2, lty = 2)
#lines(FPLSW, TPLSW, col = 3, lty = 3)
#lines(FPGL, TPGL, col = 4, lty = 4)


#-------FDR control------

FDPresultprop1 = read.table("BD30360100FDPresultprop.txt", head = FALSE)
FDPresultLiu1 = read.table("BD30360100FDPresultLiu.txt", head = FALSE)
FDPresultprop2 = read.table("BD30380100FDPresultprop.txt", head = FALSE)
FDPresultLiu2 = read.table("BD30380100FDPresultLiu.txt", head = FALSE)

FDPresultprop3 = read.table("BD30360200FDPresultprop.txt", head = FALSE)
FDPresultLiu3 = read.table("BD30360200FDPresultLiu.txt", head = FALSE)
FDPresultprop4 = read.table("BD30380200FDPresultprop.txt", head = FALSE)
FDPresultLiu4 = read.table("BD30380200FDPresultLiu.txt", head = FALSE)

MFDP = dim(FDPresultprop1)
FDPprop1 = colMeans(FDPresultprop1)
FDPLiu1 = colMeans(FDPresultLiu1)
FDPprop2 = colMeans(FDPresultprop2)
FDPLiu2 = colMeans(FDPresultLiu2)
FDPprop3 = colMeans(FDPresultprop3)
FDPLiu3 = colMeans(FDPresultLiu3)
FDPprop4 = colMeans(FDPresultprop4)
FDPLiu4 = colMeans(FDPresultLiu4)

#SDprop = c()
#SDLiu = c()
#for (i in 1 : MFDP[2]){
#	SDprop[i] = sd(FDPresultprop[, i])
#	SDLiu[i] = sd(FDPresultLiu[, i])
#}

FDPIndtau = (1 : (MFDP[2] / 3) - 1) * 3 + 1
FDPInd = (1 : (MFDP[2] / 3) - 1) * 3 + 2
FDPIndPower = (1 : (MFDP[2] / 3)) * 3

#cbind(FDPprop[FDPInd], FDPLiu[FDPInd], FDPprop[FDPIndPower], FDPLiu[FDPIndPower])

#cbind(SDprop[FDPInd], SDLiu[FDPInd], SDprop[FDPIndPower], SDLiu[FDPIndPower])


r1alpha01 = c(FDPprop1[FDPInd[2]], FDPprop1[FDPIndPower[2]], FDPLiu1[FDPInd[2]], FDPLiu1[FDPIndPower[2]])
r2alpha01 = c(FDPprop2[FDPInd[2]], FDPprop2[FDPIndPower[2]], FDPLiu2[FDPInd[2]], FDPLiu2[FDPIndPower[2]])
r3alpha01 = c(FDPprop3[FDPInd[2]], FDPprop3[FDPIndPower[2]], FDPLiu3[FDPInd[2]], FDPLiu3[FDPIndPower[2]])
r4alpha01 = c(FDPprop4[FDPInd[2]], FDPprop4[FDPIndPower[2]], FDPLiu4[FDPInd[2]], FDPLiu4[FDPIndPower[2]])

round(rbind(r1alpha01, r2alpha01, r3alpha01, r4alpha01), 3)

#----------------------------------------------------------------------
#-------------------- No need to run the followings -------------------
#----------------------------------------------------------------------

k = 2
DiffSize1 = t(FDPresultprop1[FDPInd[k]] - FDPLiu1[FDPInd[k]])
DiffSize2 = t(FDPresultprop2[FDPInd[k]] - FDPLiu2[FDPInd[k]])
DiffSize3 = t(FDPresultprop3[FDPInd[k]] - FDPLiu3[FDPInd[k]])
DiffSize4 = t(FDPresultprop4[FDPInd[k]] - FDPLiu4[FDPInd[k]])
sqrt(1000) * mean(DiffSize1) / sd(DiffSize1)
sqrt(1000) * mean(DiffSize2) / sd(DiffSize2)
sqrt(1000) * mean(DiffSize3) / sd(DiffSize3)
sqrt(1000) * mean(DiffSize4) / sd(DiffSize4)

DiffPower1 = t(FDPresultprop1[FDPIndPower[k]] - FDPLiu1[FDPIndPower[k]])
DiffPower2 = t(FDPresultprop2[FDPIndPower[k]] - FDPLiu2[FDPIndPower[k]])
DiffPower3 = t(FDPresultprop3[FDPIndPower[k]] - FDPLiu3[FDPIndPower[k]])
DiffPower4 = t(FDPresultprop4[FDPIndPower[k]] - FDPLiu4[FDPIndPower[k]])
sqrt(1000) * mean(DiffPower1) / sd(DiffPower1)
sqrt(1000) * mean(DiffPower2) / sd(DiffPower2)
sqrt(1000) * mean(DiffPower3) / sd(DiffPower3)
sqrt(1000) * mean(DiffPower4) / sd(DiffPower4)

#-------FPR control------

FPRresultprop = read.table("BD2035050FPRresultprop.txt", head = FALSE)
FPRresultLSW = read.table("BD2035050FPRresultLSW.txt", head = FALSE)

MFPR = dim(FPRresultprop)
FPRprop = colMeans(FPRresultprop)
FPRLSW = colMeans(FPRresultLSW)

FPRIndtau = (1 : (MFPR[2] / 3) - 1) * 3 + 1
FPRInd = (1 : (MFPR[2] / 3) - 1) * 3 + 2
FPRIndPower = (1 : (MFPR[2] / 3)) * 3

cbind(FPRprop[FPRInd], FPRLSW[FPRInd], FPRprop[FPRIndPower], FPRLSW[FPRIndPower])
