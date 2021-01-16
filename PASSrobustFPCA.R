estimateScores <- function (fdata, tp, efun.est){
  
  est.scores <- apply(fdata, 1, function (one.fdata) {
    in.prod <- one.fdata * efun.est
    return( sum((in.prod[-length(in.prod)] + in.prod[-1])*tp[2]/2) )
  })
  
  return(est.scores)
}

trimScores <- function (pairwise.difference.scores, trim.by = 0.05) {
  
  upper.lim <- quantile(abs(pairwise.difference.scores), 1-trim.by)
  keep <- pairwise.difference.scores[abs(pairwise.difference.scores)<=upper.lim]
  
  return(keep)
}

getU <- function (fdata, efun.est, trim=0, tp = seq(0,1, length = dim(fdata)[2])
) {
  
  library(matrixStats)
  
  N <- dim(fdata)[1]
  
  Z.index = combn(c(1:N),2)
  pairwise.diff <- t(apply(Z.index, 2, function (pair.indx) {
    pa <- fdata[pair.indx[1], ] - fdata[pair.indx[2],]
    return (pa)
  }))
  
  pairwise.diff.scores <- t(apply(efun.est, 2, estimateScores, tp = tp, fdata = pairwise.diff))
  pairwise.diff.scores.trimed <- t(apply(pairwise.diff.scores, 1, trimScores, trim.by = trim))
  pairwise.diff.scores.scaled <- t(pairwise.diff.scores.trimed / rowSds(pairwise.diff.scores.trimed))
  
  return (list(pairwise.diff.scores.scaled = pairwise.diff.scores.scaled))
}

itrEigenval <- function (PASS.ratio, u.vec,
                         ratio.old = rep(1, length(PASS.ratio)),
                         tol=1e-10,
                         max.itr=1e+20) {
  
  run <- TRUE
  itr.count <- 0
  while ( run == TRUE ) {
    
    itr.count <- itr.count + 1

    u.sum <- apply(u.vec, 1, function (eachrow) {
      sum(c(1, ratio.old) * eachrow)
    })
  
    update.weight.comp <- colMeans(u.vec / u.sum, na.rm = TRUE)
    update.weight <- update.weight.comp[1]/update.weight.comp[-1]
    lambda_temp <- PASS.ratio * update.weight
    update.diff <- sqrt(sum((lambda_temp - ratio.old)^2))
    
    ratio.old <- lambda_temp
    if (update.diff < tol) {
      run <- FALSE
    }
    if (itr.count > max.itr) {
      run <- FALSE
    }
  }
  return (ratio.old)
}

smoothPASS <- function (fdata) {
  
  Z.index <- combn(c(1:dim(fdata)[1]),2)
  Y.pair <- t(apply(Z.index, 2, function(Z.idx){
    j = Z.idx[1]; k = Z.idx[2]
    Zjk = as.matrix((fdata[j,]-fdata[k,]), nrow=1)
    return(Zjk)
  }))
  
  N = dim(Y.pair)[1]  
  M = dim(Y.pair)[2]  
  w = 1/M
  t.indx <- seq(0,1,length=M)
  
  associate.index = as.matrix(expand.grid(1:M,1:M))
  Num.ind = dim(associate.index)[1]
  G_mat = matrix(NA,nrow=N,ncol=Num.ind)
  for(i in 1:N){
    G_mat[i,] = as.vector(as.matrix(Y.pair[i,],ncol=1)%*%Y.pair[i,])
  }
  diag.index = which(associate.index[,1]==associate.index[,2])
  
  time.index = associate.index[-diag.index,]
  G_mat_reg = G_mat[,-diag.index]
  
  bivar.reg.data = data.frame(Gvec = colMeans(G_mat_reg),
                              tvec = sapply(time.index[,1],function(x) t.indx[x]),
                              svec = sapply(time.index[,2],function(x) t.indx[x]) )
  
  attach(bivar.reg.data)
  fit.obj = spm(form = Gvec~f(tvec,svec))
  detach(bivar.reg.data)
  predict.timeframe = as.matrix( expand.grid(t.indx,t.indx) )
  
  G_hat_vec = predict.spm(fit.obj,data.frame(tvec=predict.timeframe[,1],
                                             svec=predict.timeframe[,2]))
  remove(fit.obj)
  
  G_hat_mat = t( matrix(G_hat_vec,nrow=M) )
  
  return ( (G_hat_mat + t(G_hat_mat))/2 )
}

PASSfpca <- function(fdata, T_span=1, threshold=0.95){
  
  N = dim(fdata)[1]
  M = dim(fdata)[2]
  w = T_span/M   
  tp = seq(0, T_span, length=M)
  
  Z.index = combn(c(1:N),2)
  Z.Mat = apply(Z.index, 2, function(Z.idx){
    j = Z.idx[1]; k = Z.idx[2]
    tempNorm = sqrt( sum( (fdata[j,]-fdata[k,])^2 )*w )
    if (tempNorm == 0) {
      Zjk = matrix(rep(0, length(fdata[j,])), nrow = 1)
    } else {
      Zjk = as.matrix((fdata[j,]-fdata[k,])/tempNorm, nrow=1)
    }
    return(Zjk)
  })
  
  K_hat = Z.Mat%*%t(Z.Mat) * 2/(N*(N-1))
  eigen_res = eigen(K_hat)
  
  ResLists = list(evalues = w*eigen_res$values, efuns = eigen_res$vectors/sqrt(w)) 
  
  return(ResLists)
}

performPASSeigenAnalysis <- function (fdata,
                                      data.cov = 'data',
                                      trim = 0,
                                      tol = 1e-10,
                                      max.itr = 1e+20) {
  
  if (data.cov == 'data') {
    pass.results <- PASSfpca(fdata)
  } else {
    w <- 1/dim(fdata)[1]
    pass.results <- eigen(fdata)
    pass.results$values <- pass.results$values*w
    pass.results$vectors <- pass.results$vectors/sqrt(w)
    names(pass.results) <- c('evalues', 'efuns')
  }
  
  pass.eigen.fun <- pass.results$efuns
  pass.evalues <- pass.results$evalues
  Q <- sum( !((cumsum(pass.evalues) / sum(pass.evalues)) >= 1) )
  pass.evalue <- pass.evalues[1:Q]
  
  pass.ratio <- pass.evalue[-1]/pass.evalue[1]
  
  U.score <- getU(fdata, pass.eigen.fun[,1:Q], trim = trim)
  PA.scores.scaled <- U.score$pairwise.diff.scores.scaled
  u.vec.Score <- PA.scores.scaled^2
  
  itr.ratio <- c(1, itrEigenval(pass.ratio, u.vec.Score,
                                ratio.old = rep(1, Q-1)))
  
  cpe <- cumsum(itr.ratio) / sum(itr.ratio)
  
  return (list(efun = pass.eigen.fun, cpe = cpe))
}

getPASSresults <- function (fdata,
                            smooth='no smoothing',
                            trim = 0,
                            tol = 1e-10,
                            max.itr = 1e+20) {
  
  library(SemiPar)
  
  if (smooth == 'no smoothing') {
    
    pass.eigen <- performPASSeigenAnalysis(fdata)
    
  } else if (smooth == 'smooth pass cov') {
    
    pass.cov <- smoothPASS(fdata)
    pass.eigen <- performPASSeigenAnalysis(pass.cov, data.cov = 'cov')
    
    
  } else if (smooth == 'pre smooth') {
    
    N <- dim(fdata)[1]
    M <- dim(fdata)[2]
    smooth.fdata = matrix(0,nrow=N,ncol=M)
    for(i in 1:N){
      tempfdata=data.frame(tempy=fdata[i,], time.p=tp)
      attach(tempfdata)
      tempobj <- spm(tempy~f(time.p),omit.missing = T)
      smooth.fdata[i,] <- predict(tempobj,newdata=data.frame(time.p=tp))
      detach(tempfdata)
      remove(tempobj)
    }
    
    pass.eigen <- performPASSeigenAnalysis(smooth.fdata)
    
  }
  
  return (pass.eigen)
}

# Example----
M <- 100
N <- 200
tp <- seq(0, 1, length=M)
eigen.functions <- list(expression(sqrt(2)*sin(2*pi*tp)),
                        expression(sqrt(2)*cos(2*pi*tp)), 
                        expression(sqrt(2)*sin(4*pi*tp)),
                        expression(sqrt(2)*cos(4*pi*tp)))

eigen.functions.eval = sapply(eigen.functions,function(x) eval(x) )

lambda <- c(1, 1/2, 1/4, 1/8)

set.seed(300)
Gcov = diag(3/5*(lambda))
library(mvtnorm)
zeta_Mat = rmvt(n = N,sigma = Gcov, df = 5)
fdata = zeta_Mat%*%t(eigen.functions.eval)

fpca.results <- getPASSresults(fdata, trim = 0.05)
cov.results <- eigen(cov(fdata))

plot(fpca.results$efun[,1], type = 'l')
fpca.results$cpe

# With smoothing
set.seed(10)
noise = matrix(rnorm(N*M,0,1),nrow=N)
fdata.noisy <- fdata + noise

fpca.results <- getPASSresults(fdata, smooth='pre smooth')
fpca.results$cpe
plot(fpca.results$efun[,1])