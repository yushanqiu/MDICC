# compute the dominate set for the matrix aff.matrix and NR.OF.KNN
"dominate.set" <- function( aff.matrix, NR.OF.KNN ) {
  
  # create the structure to save the results
  PNN.matrix = array(0,c(nrow(aff.matrix),ncol(aff.matrix)))
  
  # sort each row of aff.matrix in descending order and saves the sorted 
  # array and a collection of vectors with the original indices
  res.sort = apply(t(aff.matrix),MARGIN=2,FUN=function(x) {return(sort(x, decreasing = TRUE, index.return = TRUE))})
  sorted.aff.matrix = t(apply(as.matrix(1:length(res.sort)),MARGIN=1,function(x) { return(res.sort[[x]]$x) }))
  sorted.indices = t(apply(as.matrix(1:length(res.sort)),MARGIN=1,function(x) { return(res.sort[[x]]$ix) }))
  
  # get the first NR.OF.KNN columns of the sorted array
  res = sorted.aff.matrix[,1:NR.OF.KNN]
  
  # create a matrix of NR.OF.KNN columns by binding vectors of 
  # integers from 1 to the number of rows/columns of aff.matrix
  inds = array(0,c(nrow(aff.matrix),NR.OF.KNN))
  inds = apply(inds,MARGIN=2,FUN=function(x) {x=1:nrow(aff.matrix)})
  
  # get the first NR.OF.KNN columns of the indices of aff.matrix
  loc = sorted.indices[,1:NR.OF.KNN]
  
  # assign to PNN.matrix the sorted indices
  PNN.matrix[(as.vector(loc)-1)*nrow(aff.matrix)+as.vector(inds)] = as.vector(res)
  
  # compute the final results and return them
  PNN.matrix = (PNN.matrix + t(PNN.matrix))/2
  
  return(PNN.matrix)
  
}

# compute the transition field of the given matrix
"transition.fields" <- function( W ) {
  
  # get any index of columns with all 0s
  zero.index = which(apply(W,MARGIN=1,FUN=sum)==0)
  
  # compute the transition fields
  W = dn(W,'ave')
  
  w = sqrt(apply(abs(W),MARGIN=2,FUN=sum)+.Machine$double.eps)
  W = W / t(apply(array(0,c(nrow(W),ncol(W))),MARGIN=2,FUN=function(x) {x=w}))
  W = W %*% t(W)
  
  # set to 0 the elements of zero.index
  W[zero.index,] = 0
  W[,zero.index] = 0
  
  return(W)
  
}

# normalizes a symmetric kernel
"dn" = function( w, type ) {
  
  # compute the sum of any column
  D = apply(w,MARGIN=2,FUN=sum)
  
  # type "ave" returns D^-1*W
  if(type=="ave") {
    D = 1 / D
    D_temp = matrix(0,nrow=length(D),ncol=length(D))
    D_temp[cbind(1:length(D),1:length(D))] = D
    D = D_temp
    wn = D %*% w
  }
  # type "gph" returns D^-1/2*W*D^-1/2
  else if(type=="gph") {
    D = 1 / sqrt(D)
    D_temp = matrix(0,nrow=length(D),ncol=length(D))
    D_temp[cbind(1:length(D),1:length(D))] = D
    D = D_temp
    wn = D %*% (w %*% D)
  }
  else {
    stop("Invalid type!")
  }
  
  return(wn)
  
}


# compute the eigenvalues and eigenvectors
"eig1" <- function( A, c = NA, isMax = NA, isSym = NA ) {
  
  # set the needed parameters
  if(is.na(c)) {
    c = dim(A)[1]
  }
  if(c>dim(A)[1]) {
    c = dim(A)[1]
  }
  if(is.na(isMax)) {
    isMax = 1
  }
  if(is.na(isSym)) {
    isSym = 1
  }
  
  # compute the eigenvalues and eigenvectors of A
  if(isSym==1) {
    eigen_A = eigen(A,symmetric=TRUE)
  }
  else {
    eigen_A = eigen(A)
  }
  v = eigen_A$vectors
  d = eigen_A$values
  
  # sort the eigenvectors
  if(isMax == 0) {
    eigen_A_sorted = sort(d,index.return=TRUE)
  }
  else {
    eigen_A_sorted = sort(d,decreasing=TRUE,index.return=TRUE)
  }
  d1 = eigen_A_sorted$x
  idx = eigen_A_sorted$ix # the index after ranking
  idx1 = idx[1:c] # the index after ranking and pick top c¡®s index
  
  # compute the results
  eigval = d[idx1]
  eigvec = Re(v[,idx1])
  eigval_full = d[idx]
  
  return(list(eigval=eigval,eigvec=eigvec,eigval_full=eigval_full))
  
}

# compute the L2 distance
"L2_distance_1" <- function( a, b ) {
  
  if(dim(a)[1] == 1) {
    a = rbind(a,rep(0,dim(a)[2]))
    b = rbind(b,rep(0,dim(b)[2]))
  }
  
  aa = apply(a*a,MARGIN=2,FUN=sum)
  bb = apply(b*b,MARGIN=2,FUN=sum)
  ab = t(a) %*% b
  d1 = apply(array(0,c(length(t(aa)),length(bb))),MARGIN=2,FUN=function(x){ x = t(aa) }) # expend col
  d2 = t(apply(array(0,c(length(t(bb)),length(aa))),MARGIN=2,FUN=function(x){ x = t(bb) }))
  d = d1 + d2 - 2 * ab
  d = Re(d)
  d = matrix(mapply(d,FUN=function(x) { return(max(max(x),0)) }),nrow=nrow(d),ncol=ncol(d))
  d_eye = array(1,dim(d))
  diag(d_eye) = 0
  d = d * d_eye
  
  return(d)
  
}

# umkl function
"umkl" = function( D, beta = NA ) {
  
  # set some parameters
  if(is.na(beta)) {
    beta = 1 / length(D)
  }
  tol = 1e-4
  u = 20
  logU = log(u)
  
  # compute Hbeta
  res_hbeta = Hbeta(D, beta)
  H = res_hbeta$H
  thisP = res_hbeta$P
  
  betamin = -Inf
  betamax = Inf
  # evaluate whether the perplexity is within tolerance
  Hdiff = H - logU
  tries = 0
  while (abs(Hdiff) > tol && tries < 30) {
    #if not, increase or decrease precision
    if (Hdiff > 0) {
      betamin = beta
      if(abs(betamax)==Inf) {
        beta = beta * 2
      }
      else {
        beta = (beta + betamax) / 2
      }
    }
    else {
      betamax = beta
      if(abs(betamin)==Inf) {
        beta = beta * 2
      }
      else {
        beta = (beta + betamin) / 2
      }
    }
    # compute the new values
    res_hbeta = Hbeta(D, beta)
    H = res_hbeta$H
    thisP = res_hbeta$P
    Hdiff = H - logU
    tries = tries + 1
  }
  
  return(thisP)
  
}

"Hbeta" = function( D, beta ) {
  
  D = (D - min(D)) / (max(D) - min(D) + .Machine$double.eps)
  P = exp(-D * beta)
  sumP = sum(P)
  H = log(sumP) + beta * sum(D * P) / sumP
  P = P / sumP
  
  return(list(H=H,P=P))
  
}


"MDICC" <- function( X, c, no.dim = NA, k = 10 ) {
  
  # set any required parameter to the defaults
  if(is.na(no.dim)) {
    no.dim = c
  }
  
  # start the clock to measure the execution time
  ptm = proc.time()
  
  # set some parameters
  NITER = 30
  num = ncol(X[[1]])
  r = -1
  beta = 0.8
  # change to dgCMatrix type.
  
  for(i in length(X)){
    X[[i]] = Matrix(X[[i]], sparse = TRUE)
  }

  X2 = X

  ##Fqian
  for(i in length(X)){
  d = apply(X2[[i]],1,sum)
     d1 = diag(d)
     X2[[i]] = solve(d1^(1)) %*% X2[[i]]
     }
  
  D_Kernels = X2
  
  # set up some parameters
  alphaK = 1 / rep(length(D_Kernels),length(D_Kernels))
  distX = array(0,c(dim(D_Kernels[[1]])[1],dim(D_Kernels[[1]])[2]))

  for (i in 1:length(D_Kernels)) {
    distX = distX + D_Kernels[[i]]
  }
 
  distX = distX / length(D_Kernels)

  res = apply(distX,MARGIN=1,FUN=function(x) return(sort(x,index.return = TRUE)))
  distX1 = array(0,c(nrow(distX),ncol(distX)))
  idx = array(0,c(nrow(distX),ncol(distX)))
  for(i in 1:nrow(distX)) {
    distX1[i,] = res[[i]]$x
    idx[i,] = res[[i]]$ix
  }

  A = array(0,c(num,num)) 
  di = distX1[,2:(k+2)] 
  rr = 0.5 * (k * di[,k+1] - apply(di[,1:k],MARGIN=1,FUN=sum))
  id = idx[,2:(k+2)] 

  if(r<=0) {
    r = mean(rr)
  }
  lambda = max(mean(rr),0)
  # A[is.nan(A)] = 0
  # A0 = (A + t(A)) / 2
  S0 = max(max(distX)) - distX
  
  cat("Network fusion.\n")
  

  
  # compute dn(!!normalization)
  S0 = dn(S0,'ave') 
  S = S0
  D0 = diag(apply(S,MARGIN=2,FUN=sum)) 
  L0 = D0 - S
  
  eig1_res = eig1(L0,c,0) 
  F_eig1 = eig1_res$eigvec
  temp_eig1 = eig1_res$eigval
  evs_eig1 = eig1_res$eigval_full
  
  # perform the iterative procedure NITER times
  converge = vector()
  for(iter in 1:NITER) {
    
    cat("Iteration: ",iter,"\n")
    
    distf = L2_distance_1(t(F_eig1),t(F_eig1))
    A = array(0,c(num,num))
    b = idx[,2:dim(idx)[2]]
    a = apply(array(0,c(num,ncol(b))),MARGIN=2,FUN=function(x){ x = 1:num })
    inda = cbind(as.vector(a),as.vector(b))
    ad = (distX[inda]+lambda*distf[inda])/2/r
    dim(ad) = c(num,ncol(b))
    
    # call the c function for the optimization
    c_input = -t(ad)
    c_output = t(ad)
    ad = t(.Call("projsplx_R",c_input,c_output))# KKT
    
    A[inda] = as.vector(ad)
    A[is.nan(A)] = 0
    A = (A + t(A)) / 2
    S = (1 - beta) * S + beta * A

    
    # After updating S
    D = diag(apply(S,MARGIN=2,FUN=sum))
    L = D - S 
    F_old = F_eig1
    eig1_res = eig1(L,c,0)
    F_eig1 = eig1_res$eigvec
    temp_eig1 = eig1_res$eigval
    ev_eig1 = eig1_res$eigval_full
    evs_eig1 = cbind(evs_eig1,ev_eig1)
    DD = vector()
    # W * S
    for (i in 1:length(D_Kernels)) {
      temp1 = (.Machine$double.eps+D_Kernels[[i]]) * (S+.Machine$double.eps)
      temp2 = 0.5*(.Machine$double.eps+D_Kernels[[i]]) * (D_Kernels[[i]]+.Machine$double.eps)
      temp = temp1 - temp2
      DD[i] = mean(apply(temp,MARGIN=2,FUN=sum))
    }
    alphaK0 = umkl(DD)
    alphaK0 = alphaK0 / sum(alphaK0)
    alphaK = (1-beta) * alphaK + beta * alphaK0
    alphaK = alphaK / sum(alphaK)
    fn1 = sum(ev_eig1[1:c])
    fn2 = sum(ev_eig1[1:(c+1)])
    converge[iter] = fn2 - fn1
    if (iter<10) {
      if (ev_eig1[length(ev_eig1)] > 0.000001) {
        lambda = 1.5 * lambda
        r = r / 1.01
      }
    }
    else {
      if(converge[iter]>converge[iter-1]) {
        S = S_old
        # if(converge[iter-1] > 0.2) {
        #   warning('Maybe you should set a larger value of c.')
        # }
        break
      }
    }
    S_old = S
    
    # compute Kbeta
    distX = D_Kernels[[1]] * alphaK[1]
    for (i in 2:length(D_Kernels)) {
      distX = distX + as.matrix(D_Kernels[[i]]) * alphaK[i]
    }
    
    # sort distX for rows
    res = apply(distX,MARGIN=1,FUN=function(x) return(sort(x,index.return = TRUE)))
    distX1 = array(0,c(nrow(distX),ncol(distX)))
    idx = array(0,c(nrow(distX),ncol(distX)))
    for(i in 1:nrow(distX)) {
      distX1[i,] = res[[i]]$x
      idx[i,] = res[[i]]$ix
    }
    
  }
  LF = F_eig1
  D = diag(apply(S,MARGIN=2,FUN=sum))
  L = D - S
  
  # compute the eigenvalues and eigenvectors of P
  eigen_L = eigen(L)
  U = eigen_L$vectors
  D = eigen_L$values
  

  # compute the execution time
  execution.time = proc.time() - ptm
  
  # create the structure with the results
  results = S
  return(results)
  
}

