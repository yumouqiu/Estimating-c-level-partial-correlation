source("PartialCorrelation1.txt")

p = 100
Sigma0 = BD2(p, 25, 1, c(0.3, 0.9))
D = diag(sqrt(runif(p, 0.1, 2)))
Sigma = D %*% Sigma0 %*% D

pc0 = solve(Sigma0)
pc = diag(sqrt(1 / diag(pc0))) %*% pc0 %*% diag(sqrt(1 / diag(pc0)))
hist(pc[pc != 0 & pc < 0.99999])

A = PC(60, p, Sigma, c0 = 0.25)

write.table(A[[1]], "BD30360100resultprop.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[2]], "BD30360100resultLiu.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[3]], "BD30360100resultLSW.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[4]], "BD30360100resultGL.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[5]], "BD30360100FDPresultprop.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[6]], "BD30360100FDPresultLiu.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[7]], "BD30360100FDPresultpropC0.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[8]], "BD30360100FDPresultLiuC0.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[9]], "BD30360100resultpropC0.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[10]], "BD30360100resultLiuC0.txt", col.names = FALSE, row.names = FALSE)


A = PC(80, p, Sigma, c0 = 0.25)

write.table(A[[1]], "BD30380100resultprop.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[2]], "BD30380100resultLiu.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[3]], "BD30380100resultLSW.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[4]], "BD30380100resultGL.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[5]], "BD30380100FDPresultprop.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[6]], "BD30380100FDPresultLiu.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[7]], "BD30380100FDPresultpropC0.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[8]], "BD30380100FDPresultLiuC0.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[9]], "BD30380100resultpropC0.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[10]], "BD30380100resultLiuC0.txt", col.names = FALSE, row.names = FALSE)




p = 200
Sigma0 = BD2(p, 50, 1, c(0.3, 0.9))
D = diag(sqrt(runif(p, 0.1, 2)))
Sigma = D %*% Sigma0 %*% D

pc0 = solve(Sigma0)
pc = diag(sqrt(1 / diag(pc0))) %*% pc0 %*% diag(sqrt(1 / diag(pc0)))
hist(pc[pc != 0 & pc < 0.99999])

A = PC(60, p, Sigma, c0 = 0.25)

write.table(A[[1]], "BD30360200resultprop.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[2]], "BD30360200resultLiu.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[3]], "BD30360200resultLSW.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[4]], "BD30360200resultGL.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[5]], "BD30360200FDPresultprop.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[6]], "BD30360200FDPresultLiu.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[7]], "BD30360200FDPresultpropC0.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[8]], "BD30360200FDPresultLiuC0.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[9]], "BD30360200resultpropC0.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[10]], "BD30360200resultLiuC0.txt", col.names = FALSE, row.names = FALSE)

A = PC(80, p, Sigma, c0 = 0.25)

write.table(A[[1]], "BD30380200resultprop.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[2]], "BD30380200resultLiu.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[3]], "BD30380200resultLSW.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[4]], "BD30380200resultGL.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[5]], "BD30380200FDPresultprop.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[6]], "BD30380200FDPresultLiu.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[7]], "BD30380200FDPresultpropC0.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[8]], "BD30380200FDPresultLiuC0.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[9]], "BD30380200resultpropC0.txt", col.names = FALSE, row.names = FALSE)
write.table(A[[10]], "BD30380200resultLiuC0.txt", col.names = FALSE, row.names = FALSE)


