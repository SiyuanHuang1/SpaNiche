#' nmf_modified
#'
#' @param dat A list of data matrices to be factorized. Each matrix represents a distinct dataset.
#' @param Wzero Initial weight matrix (W). Typically used as a starting point for the factorization process. Its dimensions are determined by the data matrices and the specified number of components/topics (k).
#' @param Hzero Initial list of H matrices. Each matrix in this list corresponds to a data matrix in dat. These matrices are also used as starting points for the factorization process.
#' @param wt A weight vector. It assigns weights to the individual data matrices in dat, indicating the importance or priority of each dataset in the factorization process.
#' @param k An integer specifying the number of components or topics to be extracted. It determines the number of columns in the W matrix and the number of rows in each H matrix.
#' @param lambda lambda
#' @param sigma Gaussian kernel parameter. Used in the computation of the Laplacian matrix, influencing how spatial information is incorporated into the factorization.
#' @param maxiter The maximum number of iterations allowed for the factorization process. This acts as a stopping criterion to prevent the algorithm from running indefinitely.
#' @param st.count Convergence counter. If the change in the reconstruction error remains below epsilon for st.count consecutive iterations, the algorithm stops, assuming it has converged.
#' @param spatial.mat A matrix containing spatial information
#' @param epsilon A threshold for the change in reconstruction error. If the relative change in error is less than epsilon for st.count consecutive iterations, convergence is assumed.
#'
#' @return
#' @export
#'
#' @examples
nmf_modified = function(
  dat = dat,
  Wzero = fit$W,
  Hzero = fit$H,
  wt = rho_v,
  k=15,
  lambda = 1,
  sigma = 0.5,
  maxiter = 50,
  st.count = 20,
  spatial.mat = as.matrix(spatialdf[,c("row","col")]),
  epsilon = 1e-04 #对loss的限制，百分比可能好一点，比起具体的值
  #seed在这个过程中需要设置吗？
  ){
  # I would like to express my gratitude to TJJ and CZH,
  # for their invaluable assistance in completing the core part of this code.
  # I consulted with them on numerous mathematical issues
  # and am deeply appreciative of their insights and explanations.
  # Additionally, I would like to extend my thanks to the R package IntNMF,
  # as our preliminary concepts were inspired by it.


  # 一些基本的变量
  n = nrow(Wzero)
  n.dat = length(Hzero)
  if (n.dat != length(wt)) stop("Number of weights must match number of data")

  # 因为涉及到很多lambda和sigma，所以是一个list
  min.f.WH <- NULL
  W.list = NULL
  for (i in 1:n.dat) {assign(paste("H",i,".list",sep = ""),NULL)}
  convergence.list = NULL

  for (lambdai in lambda) {
    for (sigmai in sigma) {
      count = 0
      iter = 1
      convergence = NULL
      W=Wzero
      H=Hzero

      laplacian_res = calculate_laplacian_matrix_with_gausskernel(
        spa_mat = spatial.mat,sigma = sigmai
      )
      laplacian_matrix = laplacian_res$laplacian
      degree_matrix = laplacian_res$degree
      adj_matrix = laplacian_res$adj

      while ((iter < maxiter) & (count < st.count)) {

        #step1: H->W
        Aminus=matrix(0,nrow = n,ncol = k)
        Aplus=matrix(0,nrow = n,ncol = k)
        for (i in 1:n.dat) {
          Aminus = Aminus + 2*wt[i] * dat[[i]] %*% t(H[[i]])
          Aplus = Aplus + 2 * wt[i] * W %*% H[[i]] %*% t(H[[i]])
        }

        # right_term = 2*lambdai*laplacian_matrix %*% W
        # Bminus = matrix(0,nrow = n,ncol = k)
        # Bplus = matrix(0,nrow = n,ncol = k)
        # Bminus = ifelse(right_term < 0,-right_term,0)
        # Bplus = ifelse(right_term >= 0,right_term,0)
        Bminus = 2*lambdai* adj_matrix %*% W
        Bplus = 2*lambdai* degree_matrix %*% W
        W = W * ((Aminus+Bminus) / pmax(Aplus+Bplus,.Machine$double.eps))

        #step2: W->H
        for (i in 1:n.dat) {
          Aminus=matrix(0,nrow = k,ncol = ncol(Hzero[[i]]))
          Aplus=matrix(0,nrow = k,ncol = ncol(Hzero[[i]]))
          Aminus = 2 * wt[i] * t(W) %*% dat[[i]]
          Aplus = 2 * wt[i] * t(W) %*% W %*% H[[i]]
          Aplus = pmax(Aplus,.Machine$double.eps)
          H[[i]] = H[[i]] * (Aminus / Aplus)
        }

        # Compute the following for each iteration
        # (i) abs.diff : Change in reconstruction error
        # (ii) f.WH    : Objective function
        # (iii) difference_level
        if (iter==1) {
          for (i in 1:n.dat) {
            assign(paste("d",i,".old",sep=""),W %*% H[[i]])
          }

          f.WH <- 0
          for (i in 1:n.dat) {
            f.WH <- f.WH + wt[i]*sum(
              (dat[[i]] - eval(parse(text=paste("d",i,".old",sep=""))))^2
            )
          }
          f.WH <- f.WH + lambdai * sum(diag(t(W) %*% laplacian_matrix %*% W))

          loss.previous = f.WH
          loss.current = NA
          difference_level = NA
          abs.diff = NA
        } else {
          for (i in 1:n.dat) {
            assign(paste("d",i,".new",sep=""),W %*% H[[i]])
          }

          abs.diff <- 0
          for (i in 1:n.dat) {
            abs.diff <- abs.diff + sum(
              abs(
                eval(parse(text=paste("d",i,".new",sep=""))) - eval(parse(text=paste("d",i,".old",sep="")))
              )
            ) / sum(
              eval(parse(text=paste("d",i,".old",sep="")))
            )
          }

          for (i in 1:n.dat) {
            assign(paste("d",i,".old",sep=""),eval(parse(text=paste("d",i,".new",sep=""))))
          }

          f.WH <- 0
          for (i in 1:n.dat) {
            f.WH <- f.WH + wt[i]*sum(
              (dat[[i]] - eval(parse(text=paste("d",i,".new",sep=""))))^2
            )
          }
          f.WH <- f.WH + lambdai * sum(diag(t(W) %*% laplacian_matrix %*% W))

          loss.current = f.WH
          difference_level = abs(loss.current - loss.previous) / loss.previous
          loss.previous = loss.current
        }

        if (!is.na(difference_level) & difference_level < epsilon) {count <- count + 1} else {count <- 0}

        #保存一些指标
        convergence <- rbind(convergence,c(iter,count,ifelse(!is.na(difference_level) & difference_level < epsilon,1,0),difference_level,abs.diff,f.WH))
        iter <- iter + 1
      } #while循环结束

      colnames(convergence) <- c("iter","count","stability","difference_level","abs.diff","f.WH")
      convergence = as.data.frame(convergence)
      convergence$lambda = lambdai
      convergence$sigma = sigmai

      # make list of all W,H, convergence
      min.f.WH <- c(min.f.WH,f.WH)
      W.list <- c(W.list,list(W))
      for (i in 1:n.dat) {
        assign(paste("H",i,".list",sep=""),
               c(
                 eval(parse(text=paste("H",i,".list",sep=""))),list(H[[i]])
               )
        )
      }
      convergence.list <- c(convergence.list,list(convergence))

    }
  } #lambda和sigma循环完了
  H.list <- NULL
  for (i in 1:n.dat) H.list <- c(H.list,list(eval(parse(text=paste("H",i,".list",sep="")))))
  names(H.list) <- paste("H",1:n.dat,".list",sep="")


  # Select the factorization (W and H) that leads to the lowest objective function 'f.WH'
  convergence <- convergence.list[[which.min(min.f.WH)]]
  W <- W.list[[which.min(min.f.WH)]]
  H <- NULL
  for (i in 1:n.dat){
    assign(paste("H",i,sep=""),H.list[[i]][[which.min(min.f.WH)]])
    H <- c(H,list(eval(parse(text=paste("H",i,sep="")))))
  }
  names(H) <- paste("H",1:n.dat,sep="")

  # Compute cluster membership
  clusters <- apply(W,1,which.max)
  # Output the following results
  output <- list(W=W,H=H,convergence=convergence,min.f.WH=min.f.WH,clusters=clusters)
  return(output)
}
