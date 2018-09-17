library(lars)
library(MASS)
library(glasso)
library(genlasso)
library(scalreg)

A = read.csv("PET.csv", head = FALSE)
AD = read.csv("AD Index.csv", head = FALSE)
NC = read.csv("NL Index.csv", head = FALSE)
#Vox = read.csv("vox mean2.csv", head = FALSE)
avoi = c(3, 4, 7, 8, 23 : 28, 31, 32, 35 : 40, 49 : 56, 59 : 62, 67, 68, 81 : 90)
group = c(rep(c(1, 2), 6), 3, 4, rep(c(7, 8), 2), rep(c(5, 6), 3), 7, 8, rep(c(3, 4), 3), rep(c(7, 8), 5))
rename = c(3, 4, 7, 8, 23 : 28, 31, 32, 59 : 62, 67, 68, 35, 36, 49 : 54, 81 : 90, 55, 56, 37 : 40)

ADdata = t(A[, AD$V1])
NCdata = t(A[, NC$V1])
MCIdata = t(A[, -union(AD$V1, NC$V1)])


X = ADdata[, rename]  # This is for AD group
                      # Use X = NCdata[, rename] for the NC group

n = dim(X)[1]; p = dim(X)[2]
t0 = 2; tau = seq(0, 3.5, 0.01); smax = n / 2; lentau = length(tau); c0 = 0.25 
IndMatrix = matrix(1, p, p) - diag(rep(1, p))
Eresidual = matrix(0, n, p)
CoefMatrix = matrix(0, p, p - 1)
meanX = colMeans(X)
X = t(t(X) - meanX)
XS = matrix(0, n, p)
for (i in 1 : p){
	XS[, i] = X[, i] / sd(X[, i])
}

for (i in 1 : p){
	#out = genlasso(X[, i], X = X[, -i], D = diag(1, p - 1))
	#Coef = coef(out,  lambda = sqrt(2 * var(X[, i]) * n * log(p * log(p) / sqrt(n))))
	#Predict = predict(out, lambda = sqrt(2 * var(X[, i]) * n * log(p * log(p) / sqrt(n))), Xnew = X[, -i])
	#CoefMatrix[i, ] = t(Coef$beta)
	#Eresidual[, i] = X[, i] - Predict$fit

	out = scalreg(X = XS[, -i], y = X[, i], lam0 = sqrt(2 * 2.01 * log(p * (log(p))^(1.5) / sqrt(n)) / n))
	Eresidual[, i] = out$residuals
	CoefMatrix[i, ] = out$coefficients / apply(X[, -i], 2, sd)
}

CovRes = t(Eresidual) %*% Eresidual / n
Est = matrix(1, p, p)

for (i in 1 : (p - 1)){
	for (j in (i + 1) : p){
		temp = Eresidual[, i] * Eresidual[, j] + Eresidual[, i]^2 * CoefMatrix[j, i] + Eresidual[, j]^2 * CoefMatrix[i, j - 1]
		Est[i, j] = mean(temp) / sqrt(diag(CovRes)[i] * diag(CovRes)[j])
		Est[j, i] = Est[i, j]
	}
}

EstThresh = Est * ( abs(Est) >= (t0 * sqrt(log(p) / n) * IndMatrix) )
kappa = (n / 3) * mean( colSums(Eresidual^4) / (colSums(Eresidual^2))^2 )
CovX = t(X) %*% X / n - matrix(colMeans(X), p, 1) %*% matrix(colMeans(X), 1, p)

GL0 = glassopath(CovX, trace = 0)
rho0 = GL0$rholist
indexGL = c()
for (i in 1 : length(rho0)){
	precision0 = GL0$wi[, , i]
	indexGL = c(indexGL, sum(precision0 != 0))
}
if (length(indexGL <= p) > 0) rhomax = min(rho0[indexGL <= p])
if (length(indexGL <= p) == 0) rhomax = max(rho0)
rho = seq(rho0[1] / 2, rhomax, length.out = lentau)
GL = glassopath(CovX, rholist = rho, trace = 0)

resprop = list(); resLiu = list(); resLSW = list(); resGL = list()
rejectprop = c(); rejectLiu = c(); rejectLSW = c(); rejectGL = c()
rejectpropC0 = c(); rejectLiuC0 = c()
for (i in 1 : lentau){
	Threshold = tau[i] * sqrt(kappa * log(p) / n) * (1 - EstThresh^2)
	SRec = 1 * (abs(Est) > Threshold); NoNSRec = 1 * (SRec == 0)
	resprop[[i]] = which(SRec == 1, arr.ind = TRUE)
	rejectprop = c(rejectprop, max(1, (sum(SRec) - p)))
	SRecC0 = 1 * (abs(Est) - c0 > Threshold)
	rejectpropC0 = c(rejectpropC0, max(1, (sum(SRecC0) - p)))

	ThresholdLiu = tau[i] * sqrt(log(p) / n) * IndMatrix
	SRecLiu = 1 * (abs(Est) > ThresholdLiu); NoNSRecLiu = 1 * (SRecLiu == 0)
	resLiu[[i]] = which(SRecLiu == 1, arr.ind = TRUE)
	rejectLiu = c(rejectLiu, max(1, (sum(SRecLiu) - p)))
	SRecLiuC0 = 1 * (abs(Est) - c0 > ThresholdLiu)
	rejectLiuC0 = c(rejectLiuC0, max(1, (sum(SRecLiuC0) - p)))

	ThresholdLSW = tau[i] * sqrt(log(p) / n) * IndMatrix
	SRecLSW = 1 * (abs(CovX) > ThresholdLSW); NoNSRecLSW = 1 * (SRecLSW == 0)
	resLSW[[i]] = which(SRecLSW == 1, arr.ind = TRUE)
	rejectLSW = c(rejectLSW, max(1, (sum(SRecLSW) - p)))

	SRecGL = 1 * (GL$wi[, , i] != 0); NoNSRecGL = 1 * (GL$wi[, , i] == 0)
	resGL[[i]] = which(SRecGL == 1, arr.ind = TRUE)
	rejectGL = c(rejectGL, max(1, (sum(SRecGL) - p)))
}

FDPprop = 2 * p * (p - 1) * ( 1 - pnorm( tau * sqrt(log(p)) ) ) / rejectprop
FDPLiu = 2 * p * (p - 1) * ( 1 - pnorm( tau * sqrt(log(p)) ) ) / rejectLiu
alpha = c(0.01, 0.05, 0.1, 0.15, 0.2)
lenalpha = length(alpha)
FDPresprop = list(); FDPresLiu = list()
for (i in 1 : lenalpha){
	if (sum(FDPprop <= alpha[i]) > 0) tauprop = min(c(2, tau[FDPprop <= alpha[i]]))
	if (sum(FDPprop <= alpha[i]) == 0) tauprop = 2
	Threshold = tauprop * sqrt(kappa * log(p) / n) * (1 - EstThresh^2)
	SRec = 1 * (abs(Est) > Threshold); NoNSRec = 1 * (SRec == 0)
	FDPresprop[[i]] = which(SRec == 1, arr.ind = TRUE)

	if (sum(FDPLiu <= alpha[i]) > 0) tauLiu = min(c(2, tau[FDPLiu <= alpha[i]]))
	else tauLiu = 2
	ThresholdLiu = tauLiu * sqrt(log(p) / n) * IndMatrix
	SRecLiu = 1 * (abs(Est) > ThresholdLiu); NoNSRecLiu = 1 * (SRecLiu == 0)
	FDPresLiu[[i]] = which(SRecLiu == 1, arr.ind = TRUE)
}

FDPpropC0 = 2 * ( sqrt(n) * p * ( 1 - pnorm( tau * sqrt(log(p)) ) ) + p * (p - 1 - sqrt(n)) * ( 1 - pnorm( sqrt(n) * c0 / sqrt(kappa) + tau * sqrt(log(p)) ) ) ) / rejectpropC0
FDPLiuC0 = 2 * p * (p - 1) * ( 1 - pnorm( tau * sqrt(log(p)) ) ) / rejectLiuC0
FDPrespropC0 = list(); FDPresLiuC0 = list()
for (i in 1 : lenalpha){
	if (sum(FDPpropC0 <= alpha[i]) > 0) taupropC0 = min(c(2, tau[FDPpropC0 <= alpha[i]]))
	else taupropC0 = 2
	ThresholdC0 = taupropC0 * sqrt(kappa * log(p) / n) * (1 - EstThresh^2)
	SRecC0 = 1 * (abs(Est) - c0 > ThresholdC0); NoNSRecC0 = 1 * (SRecC0 == 0)
	FDPrespropC0[[i]] = which(SRecC0 == 1, arr.ind = TRUE)

	if (sum(FDPLiuC0 <= alpha[i]) > 0) tauLiuC0 = min(c(2, tau[FDPLiuC0 <= alpha[i]]))
	else tauLiuC0 = 2
	ThresholdLiuC0 = tauLiuC0 * sqrt(log(p) / n) * IndMatrix
	SRecLiuC0 = 1 * (abs(Est) - c0 > ThresholdLiuC0); NoNSRecLiuC0 = 1 * (SRecLiuC0 == 0)
	FDPresLiuC0[[i]] = which(SRecLiuC0 == 1, arr.ind = TRUE)
}

#--- plot ----

#ind60 = which(rejectGL == (2 * 90))[1]
#dis = resGL[[ind60]]
#plot(1, 1, xlim = c(1, p), ylim = c(p, 1), lwd = 1, pch = 0, xlab = 'AVOIs', ylab = 'AVOIs', main = "GLasso for AD, 90 Different Edges")
#for (i in 1 : p){
#	for (j in 1 : p){
#		points(i, j, lwd = 1, pch = 0)
#	}
#}

#for (i in 1 : dim(dis)[1]){
#	if (as.numeric(dis[i, 1]) != as.numeric(dis[i, 2])) points(as.numeric(dis[i, 1]), as.numeric(dis[i, 2]), pch = 15)
#}


k = 3
dis = FDPresprop[[k]]
(dim(dis)[1] - p) / 2
plot(1, 1, xlim = c(1, p), ylim = c(p, 1), lwd = 1, pch = 0, xlab = 'AVOIs', ylab = 'AVOIs', main = "AD Brain Connectivity")
for (i in 1 : p){
	for (j in 1 : p){
		points(i, j, lwd = 1, pch = 0)
	}
}

for (i in 1 : dim(dis)[1]){
	if (as.numeric(dis[i, 1]) != as.numeric(dis[i, 2])) points(as.numeric(dis[i, 1]), as.numeric(dis[i, 2]), col = 2, pch = 15)
}

indGL = which(rejectGL == (dim(dis)[1] - p))[1]
disGL = resGL[[indGL]]
for (i in 1 : dim(disGL)[1]){
	if (as.numeric(disGL[i, 1]) != as.numeric(disGL[i, 2])) points(as.numeric(disGL[i, 1]), as.numeric(disGL[i, 2]), col = 3, pch = 20)
}

disLiu = FDPresLiu[[k]]
(dim(disLiu)[1] - p) / 2
distemp = c(); disLiutemp = c()
for (i in 1 : dim(disLiu)[1]){
	if (as.numeric(disLiu[i, 1]) != as.numeric(disLiu[i, 2])) points(as.numeric(disLiu[i, 1]), as.numeric(disLiu[i, 2]), col = 4, pch = 4)
}



for (i in 1 : dim(dis)[1]){
	distemp[i] = 100 * dis[i, 1] + dis[i, 2]
}
for (i in 1 : dim(disLiu)[1]){
	disLiutemp[i] = 100 * disLiu[i, 1] + disLiu[i, 2]
}
Difftemp = setdiff(disLiutemp, distemp)
c(length(Difftemp), dim(disLiu)[1] - dim(dis)[1])
points(floor(Difftemp / 100), Difftemp - 100 * floor(Difftemp / 100), col = 4, pch = 4)


#------------------ c-level graph -----------------------

k = 3
dis = FDPrespropC0[[k]]
(dim(dis)[1] - p) / 2
plot(1, 1, xlim = c(1, p), ylim = c(p, 1), lwd = 1, pch = 0, xlab = 'AVOIs', ylab = 'AVOIs', main = "AD Brain Connectivity, 0.25-level Graph")
for (i in 1 : p){
	for (j in 1 : p){
		points(i, j, lwd = 1, pch = 0)
	}
}

for (i in 1 : dim(dis)[1]){
	if (as.numeric(dis[i, 1]) != as.numeric(dis[i, 2])) points(as.numeric(dis[i, 1]), as.numeric(dis[i, 2]), col = 2, pch = 15)
}

disLiu = FDPresLiuC0[[k]]
(dim(disLiu)[1] - p) / 2
distemp = c(); disLiutemp = c()

for (i in 1 : dim(disLiu)[1]){
	if (as.numeric(disLiu[i, 1]) != as.numeric(disLiu[i, 2])) points(as.numeric(disLiu[i, 1]), as.numeric(disLiu[i, 2]), col = 4, pch = 4)
}


#for (i in 1 : dim(dis)[1]){
#	distemp[i] = 100 * dis[i, 1] + dis[i, 2]
#}
#for (i in 1 : dim(disLiu)[1]){
#	disLiutemp[i] = 100 * disLiu[i, 1] + disLiu[i, 2]
#}
#Difftemp = setdiff(distemp, disLiutemp)
#c(length(Difftemp), dim(disLiu)[1] - dim(dis)[1])
#points(floor(Difftemp / 100), Difftemp - 100 * floor(Difftemp / 100), col = 4, pch = 12)





#------------------ compare between correlation and partial correlation -----------------------

corr = cor(X)
plot(1, 1, xlim = c(1, p), ylim = c(p, 1), lwd = 1, pch = 0, xlab = 'AVOIs', ylab = 'AVOIs', main = "|Correlation| large than 0.3")
for (i in 1 : p){
	for (j in 1 : p){
		points(i, j, lwd = 1, pch = 0)
		if ( (abs(corr[i, j]) > 0.3) && (i != j) ) points(i, j, col = 2, pch = 15)
	}
}
lines(c(0.5, 12.5), c(12.5, 12.5), col = 3, lwd = 2)
lines(c(12.5, 12.5), c(0.5, 12.5), col = 3, lwd = 2)
lines(c(12.5, 20.5), c(12.5, 12.5), col = 5, lwd = 2)
lines(c(12.5, 12.5), c(12.5, 20.5), col = 5, lwd = 2)
lines(c(12.5, 20.5), c(20.5, 20.5), col = 5, lwd = 2)
lines(c(20.5, 20.5), c(12.5, 20.5), col = 5, lwd = 2)
lines(c(20.5, 26.5), c(20.5, 20.5), col = 6, lwd = 2)
lines(c(20.5, 20.5), c(20.5, 26.5), col = 6, lwd = 2)
lines(c(20.5, 26.5), c(26.5, 26.5), col = 6, lwd = 2)
lines(c(26.5, 26.5), c(20.5, 26.5), col = 6, lwd = 2)
lines(c(26.5, 42.5), c(26.5, 26.5), col = 7, lwd = 2)
lines(c(26.5, 26.5), c(26.5, 42.5), col = 7, lwd = 2)

plot(1, 1, xlim = c(1, p), ylim = c(p, 1), lwd = 1, pch = 0, xlab = 'AVOIs', ylab = 'AVOIs', main = "|Partial Correlation| large than 0.3")
for (i in 1 : p){
	for (j in 1 : p){
		points(i, j, lwd = 1, pch = 0)
		if ( (abs(Est[i, j]) > 0.3) && (i != j) ) points(i, j, col = 4, pch = 15)
	}
}
lines(c(0.5, 12.5), c(12.5, 12.5), col = 3, lwd = 2)
lines(c(12.5, 12.5), c(0.5, 12.5), col = 3, lwd = 2)
lines(c(12.5, 20.5), c(12.5, 12.5), col = 5, lwd = 2)
lines(c(12.5, 12.5), c(12.5, 20.5), col = 5, lwd = 2)
lines(c(12.5, 20.5), c(20.5, 20.5), col = 5, lwd = 2)
lines(c(20.5, 20.5), c(12.5, 20.5), col = 5, lwd = 2)
lines(c(20.5, 26.5), c(20.5, 20.5), col = 6, lwd = 2)
lines(c(20.5, 20.5), c(20.5, 26.5), col = 6, lwd = 2)
lines(c(20.5, 26.5), c(26.5, 26.5), col = 6, lwd = 2)
lines(c(26.5, 26.5), c(20.5, 26.5), col = 6, lwd = 2)
lines(c(26.5, 42.5), c(26.5, 26.5), col = 7, lwd = 2)
lines(c(26.5, 26.5), c(26.5, 42.5), col = 7, lwd = 2)

k = 1
corrVec = c()
EstVec = c()
for (i in 1 : (p - 1)){
	for (j in (i + 1) : p){
		corrVec[k] = corr[i, j]
		EstVec[k] = Est[i, j]
		k = k + 1
	}
}
EstVec[EstVec > 1] = 1

hist(corrVec, breaks = 15, col = rgb(1, 1, 0, 0.7), ylim = c(0, 200), xlim = c(-1, 1), main = "Histograms of correlation and partial correlation", xlab = "Correlation / Partial correlation")
par(new = TRUE)
hist(EstVec, breaks = 15, col = rgb(0, 1, 1, 0.4), ylim = c(0, 200), xlim = c(-1, 1), main = "", xlab = "", ylab = "")
legend("topright", c("Correlation", "Partial correlation"), col = c(rgb(1, 1, 0, 0.7), rgb(0, 1, 1, 0.4)), pch = c(15, 15), pt.cex = c(2, 2))
 

quan = 0.85
corr.cut = quantile(abs(corrVec), quan)
Est.cut = quantile(abs(EstVec), quan)

plot(1, 1, xlim = c(1, p), ylim = c(p, 1), lwd = 1, pch = 0, xlab = 'AVOIs', ylab = 'AVOIs', main = "Upper 15% correlations in absolute value")
for (i in 1 : p){
	for (j in 1 : p){
		points(i, j, lwd = 1, pch = 0)
		if ( (abs(corr[i, j]) > corr.cut) && (i != j) ) points(i, j, col = 2, pch = 15)
	}
}
lines(c(0.5, 12.5), c(12.5, 12.5), col = 3, lwd = 2)
lines(c(12.5, 12.5), c(0.5, 12.5), col = 3, lwd = 2)
lines(c(12.5, 20.5), c(12.5, 12.5), col = 5, lwd = 2)
lines(c(12.5, 12.5), c(12.5, 20.5), col = 5, lwd = 2)
lines(c(12.5, 20.5), c(20.5, 20.5), col = 5, lwd = 2)
lines(c(20.5, 20.5), c(12.5, 20.5), col = 5, lwd = 2)
lines(c(20.5, 26.5), c(20.5, 20.5), col = 6, lwd = 2)
lines(c(20.5, 20.5), c(20.5, 26.5), col = 6, lwd = 2)
lines(c(20.5, 26.5), c(26.5, 26.5), col = 6, lwd = 2)
lines(c(26.5, 26.5), c(20.5, 26.5), col = 6, lwd = 2)
lines(c(26.5, 42.5), c(26.5, 26.5), col = 7, lwd = 2)
lines(c(26.5, 26.5), c(26.5, 42.5), col = 7, lwd = 2)

plot(1, 1, xlim = c(1, p), ylim = c(p, 1), lwd = 1, pch = 0, xlab = 'AVOIs', ylab = 'AVOIs', main = "Upper 15% partial correlations in absolute value")
for (i in 1 : p){
	for (j in 1 : p){
		points(i, j, lwd = 1, pch = 0)
		if ( (abs(Est[i, j]) > Est.cut) && (i != j) ) points(i, j, col = 4, pch = 15)
	}
}
lines(c(0.5, 12.5), c(12.5, 12.5), col = 3, lwd = 2)
lines(c(12.5, 12.5), c(0.5, 12.5), col = 3, lwd = 2)
lines(c(12.5, 20.5), c(12.5, 12.5), col = 5, lwd = 2)
lines(c(12.5, 12.5), c(12.5, 20.5), col = 5, lwd = 2)
lines(c(12.5, 20.5), c(20.5, 20.5), col = 5, lwd = 2)
lines(c(20.5, 20.5), c(12.5, 20.5), col = 5, lwd = 2)
lines(c(20.5, 26.5), c(20.5, 20.5), col = 6, lwd = 2)
lines(c(20.5, 20.5), c(20.5, 26.5), col = 6, lwd = 2)
lines(c(20.5, 26.5), c(26.5, 26.5), col = 6, lwd = 2)
lines(c(26.5, 26.5), c(20.5, 26.5), col = 6, lwd = 2)
lines(c(26.5, 42.5), c(26.5, 26.5), col = 7, lwd = 2)
lines(c(26.5, 26.5), c(26.5, 42.5), col = 7, lwd = 2)
