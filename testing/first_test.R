library(expm)
library(spExpM)
library(tidyverse)

M <- as.matrix(1/dist(1:10000))
M[M<0.25] <- 0
# M <- M/4
diag(M) <- -rowSums(M)
M <- Matrix(M)
# v <- rbind(
#   c(1,rep(0,999)),
#   c(0,1,rep(0,998))
# )
v <- diag(10000)
# x1 <- expAtv(t(M), v)$eAtv
x2 <- SPS_v_exp_spQ(v, M, 1.0e-16)
max(abs(t(x1)-x2))

xxx <- MUSPS2_v_exp_spQ(v, M, 1.0e-16, c(1, 10, 100, 200, 300))
xxx2 <- MUSPS2r_v_exp_spQ(v, M, 1.0e-16, c(1, 10, 100, 200, 300))

x3 <- SPS_v_exp_spQ(matrix(xxx[1,],1,10000), 9*M, 1.0e-16)
max(abs(xxx-xxx2))
max(abs(xxx[1,]-x2))
max(abs(xxx[2,]-x3))

