# Here is some even older code (ca. 2017) which does not depend on lavaan so
# much.

library(maxLik) # so that AsymVCOVLamdaPhiTau  function works
library(lavaan) # is needed so that Ploglik and GrPloglik functions work


# output : the value of the pairwise loglikelihood written interms of lamdas, phis, taus
# input :
# free.par :  a vector with all the free parameters given in the order
#            loadings (given columnwise), phis (given columnwise) #!!!!! have to check what happens with more than two factors
#            and finally thresholds (withing each item in ascending order)
# model : an model object as should be provided in lavaan
#         make sure that it explicitly includes the variances and covariances of
#         latent variables and all thresholds for all variables (excluding
#         the redundant ones -Inf and Inf)
#         also write NA in front of the first indicator of each latent variable
#         to denote that it is free to be estimated
#         example

# model <- "
# ksi1 =~ NA*x1 + x2 + x3 + x4
# ksi2 =~ NA*x4 + x5 + x6
#
# ksi1 ~~ ksi2
# ksi1 ~~ 1*ksi1
# ksi2 ~~ 1*ksi2
#
# x1 | t1 + t2 +t3
# x2 | t1 + t2 +t3
# x3 | t1 + t2 +t3
# x4 | t1 + t2 +t3
# x5 | t1 + t2 +t3
# x6 | t1 + t2 +t3
# "

# Data : row data in matrix object  if casewise=TRUE
#       otherwise, n.xixj.ab vector  which gives all the bivariate frequencies
#       for the all pairs  xi-xj where a the categories of xi and b of xj
#       index a runs fastest of all, followed by b, followed by j and lastly by i
#       n.xixj.ab is the output of
#       unlist( lavaan:::pairwiseTables(data=data.sim, no.x=no.x)$pairTables )
# casewise : if TRUE the value of the pairwise loglikelihood for each observation
#            is computed, if FALSE the value of the pairwise loglikelihood for
#            the whole sample is computed
# NOTE: all non-redundant thresholds are assumed free in this code

Ploglik <- function(free.par, model, Data, casewise = FALSE) {
  ParTable <- lavaanify(model)
  names.x <- lavaanNames(ParTable, type = "ov")
  no.x <- length(names.x)
  no.ksi <- length(lavaanNames(ParTable, type = "lv"))

  LisrelMM <- as.data.frame(lavaan:::representation.LISREL(ParTable))

  Lx.matrix <- matrix(NA, nrow = no.x, ncol = no.ksi)
  Lx.free.index <- as.matrix(LisrelMM[LisrelMM == "lambda" & ParTable$free != 0, 2:3])
  if (nrow(Lx.free.index) > 0) {
    Lx.matrix[Lx.free.index] <- free.par[1:nrow(Lx.free.index)]
  }
  Lx.fix.index <- as.matrix(LisrelMM[LisrelMM == "lambda" & ParTable$free == 0, 2:3])
  Lx.matrix[Lx.fix.index] <- ParTable$ustart[LisrelMM == "lambda" & ParTable$free == 0]
  Lx.matrix[is.na(Lx.matrix)] <- 0

  no.free.par.used <- nrow(Lx.free.index)

  Phi.matrix <- matrix(NA, no.ksi, no.ksi)
  Phi.free.index <- as.matrix(LisrelMM[LisrelMM == "psi" & ParTable$free != 0, 2:3])
  if (nrow(Phi.free.index) > 0) {
    free.phis <- free.par[(no.free.par.used + 1):
    (no.free.par.used + nrow(Phi.free.index))]

    Phi.matrix[Phi.free.index] <- free.phis
    Phi.matrix[cbind(Phi.free.index[, 2], Phi.free.index[, 1])] <- free.phis
    no.free.par.used <- no.free.par.used + nrow(Phi.free.index)
  }
  Phi.fix.index <- as.matrix(LisrelMM[LisrelMM == "psi" & ParTable$free == 0, 2:3])
  Phi.matrix[Phi.fix.index] <- ParTable$ustart[LisrelMM == "psi" & ParTable$free == 0]


  all.thres <- free.par[(no.free.par.used + 1):length(free.par)]

  th.idx <- ParTable$lhs[ParTable$op == "|"]
  for (i in 1:no.x) {
    th.idx[th.idx == names.x[i]] <- i
  }
  th.idx <- as.numeric(th.idx)

  Lx.Phi.tLx <- Lx.matrix %*% Phi.matrix %*% t(Lx.matrix)
  rho.xixj <- Lx.Phi.tLx[lower.tri(Lx.Phi.tLx)]

  # the following copied from Yve's code so that the function can be maximized by any R version through maxLik function or nlminb
  if (any(abs(rho.xixj) > 1)) {
    return(+Inf)
  }


  LONG2 <- lavaan:::LongVecTH.Rho(
    no.x = no.x,
    all.thres = all.thres,
    index.var.of.thres = th.idx,
    rho.xixj = rho.xixj
  )

  LONG1 <- lavaan:::LongVecInd(
    no.x = no.x,
    all.thres = all.thres,
    index.var.of.thres = th.idx
  )

  PI <- lavaan:::pairwiseExpProbVec(ind.vec = LONG1, th.rho.vec = LONG2)
  if (!casewise) {
    LogLik <- sum(Data * log(PI))
    res <- LogLik
  } else {
    nobs <- nrow(Data)
    LogLik.casewise <- rep(NA, nobs)
    no.cat.x <- rep(NA, no.x)
    for (i in 1:no.x) {
      no.cat.x[i] <- sum(ParTable$lhs == names.x[i] & ParTable$op == "|") + 1
    }
    LONG3 <- vector("list", length = 4)
    names(LONG3) <- c("cat.idx.xi", "cat.idx.xj", "idx.xi", "idx.xj")
    LONG3$cat.idx.xi <-
      LONG1$index.thres.var1.of.pair[!(LONG1$index.thres.var1.of.pair == 0 |
        LONG1$index.thres.var2.of.pair == 0)]
    LONG3$cat.idx.xj <-
      LONG1$index.thres.var2.of.pair[!(LONG1$index.thres.var1.of.pair == 0 |
        LONG1$index.thres.var2.of.pair == 0)]
    LONG3$idx.xi <-
      LONG1$index.var1.of.pair[!(LONG1$index.thres.var1.of.pair == 0 |
        LONG1$index.thres.var2.of.pair == 0)]
    LONG3$idx.xj <-
      LONG1$index.var2.of.pair[!(LONG1$index.thres.var1.of.pair == 0 |
        LONG1$index.thres.var2.of.pair == 0)]

    comb.idx <- combn(no.x, 2)
    no.comb <- ncol(comb.idx)
    no.bifreq <- length(PI)
    vec.bifreq.idx <- c(1:no.bifreq)

    for (i in 1:nobs) {
      obs <- Data[i, ]
      pair.xixj.obs <- matrix(as.numeric(obs[comb.idx]), ncol = 2, byrow = TRUE)
      n.xixj.ab <- rep(0, no.bifreq)

      for (j in 1:no.comb) {
        bifreq.idx <- vec.bifreq.idx[(LONG3$idx.xi == comb.idx[1, j] &
          LONG3$idx.xj == comb.idx[2, j] &
          LONG3$cat.idx.xi == pair.xixj.obs[j, 1] &
          LONG3$cat.idx.xj == pair.xixj.obs[j, 2])]
        n.xixj.ab[bifreq.idx] <- 1
      }

      LogLik.casewise[i] <- sum(n.xixj.ab * log(PI))
    }
    res <- LogLik.casewise
  }
  res
}

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




#### Gradient function
# output : the value of the Gradient pairwise loglikelihood with respect to
#         lamdas, phis, and thresholds in this order (exactly same order as
#         given in free.par vector)
# input as above

GrPloglik <- function(free.par, model, Data, casewise = FALSE) {
  ParTable <- lavaanify(model)
  names.x <- lavaanNames(ParTable, type = "ov")
  no.x <- length(names.x)
  no.ksi <- length(lavaanNames(ParTable, type = "lv"))

  LisrelMM <- as.data.frame(lavaan:::representation.LISREL(ParTable))

  Lx.matrix <- matrix(NA, nrow = no.x, ncol = no.ksi)
  Lx.free.index <- as.matrix(LisrelMM[LisrelMM == "lambda" & ParTable$free != 0, 2:3])
  if (nrow(Lx.free.index) > 0) {
    Lx.matrix[Lx.free.index] <- free.par[1:nrow(Lx.free.index)]
    idx.free.lamdas <- c(1:length(Lx.matrix))[!is.na(Lx.matrix)] # is needed further below
  }
  Lx.fix.index <- as.matrix(LisrelMM[LisrelMM == "lambda" & ParTable$free == 0, 2:3])
  Lx.matrix[Lx.fix.index] <- ParTable$ustart[LisrelMM == "lambda" & ParTable$free == 0]
  Lx.matrix[is.na(Lx.matrix)] <- 0

  no.free.par.used <- nrow(Lx.free.index)

  Phi.matrix <- matrix(NA, no.ksi, no.ksi)
  Phi.free.index <- as.matrix(LisrelMM[LisrelMM == "psi" & ParTable$free != 0, 2:3])
  if (nrow(Phi.free.index) > 0) {
    free.phis <- free.par[(no.free.par.used + 1):
    (no.free.par.used + nrow(Phi.free.index))]
    Phi.matrix[cbind(Phi.free.index[, 2], Phi.free.index[, 1])] <- free.phis
    idx.free.phis <- c(1:length(Phi.matrix))[!is.na(Phi.matrix)] # is needed further below
    Phi.matrix[Phi.free.index] <- free.phis
    no.free.par.used <- no.free.par.used + nrow(Phi.free.index)
  }
  Phi.fix.index <- as.matrix(LisrelMM[LisrelMM == "psi" & ParTable$free == 0, 2:3])
  Phi.matrix[Phi.fix.index] <- ParTable$ustart[LisrelMM == "psi" & ParTable$free == 0]

  all.thres <- free.par[(no.free.par.used + 1):length(free.par)]

  th.idx <- ParTable$lhs[ParTable$op == "|"]
  for (i in 1:no.x) {
    th.idx[th.idx == names.x[i]] <- i
  }
  th.idx <- as.numeric(th.idx)

  Lx.Phi.tLx <- Lx.matrix %*% Phi.matrix %*% t(Lx.matrix)
  rho.xixj <- Lx.Phi.tLx[lower.tri(Lx.Phi.tLx)]

  # the following copied from Yve's code so that the function can be maximized by any R version through maxLik function or nlminb
  if (any(abs(rho.xixj) > 1)) {
    return(+Inf)
  }

  LONG2 <- lavaan:::LongVecTH.Rho(
    no.x = no.x,
    all.thres = all.thres,
    index.var.of.thres = th.idx,
    rho.xixj = rho.xixj
  )

  LONG1 <- lavaan:::LongVecInd(
    no.x = no.x,
    all.thres = all.thres,
    index.var.of.thres = th.idx
  )

  PI <- lavaan:::pairwiseExpProbVec(ind.vec = LONG1, th.rho.vec = LONG2)

  .MLIST <- list(lambda = Lx.matrix, psi = Phi.matrix)

  if (nrow(Lx.free.index) > 0) {
    der.sigmavec.to.lamda <- lavaan:::derivative.sigma.LISREL(
      m = "lambda",
      idx = idx.free.lamdas, MLIST = .MLIST
    )
  } else {
    der.sigmavec.to.lamda <- c()
  }

  if (nrow(Phi.free.index) > 0) {
    der.sigmavec.to.phi <- lavaan:::derivative.sigma.LISREL(
      m = "psi",
      idx = idx.free.phis, MLIST = .MLIST
    )
  } else {
    der.sigmavec.to.phi <- c()
  }


  if (nrow(Lx.free.index) > 0 | nrow(Phi.free.index) > 0) {
    der.sigmavec.to.lamda.and.phi <- cbind(
      der.sigmavec.to.lamda,
      der.sigmavec.to.phi
    )

    cov.idx <- lavaan:::vech.idx(no.x)
    covd.idx <- lavaan:::vech.idx(no.x, diag = FALSE)
    cor.idx <- match(covd.idx, cov.idx)
    derRhotoLamdaAndPhi <- der.sigmavec.to.lamda.and.phi[cor.idx, , drop = FALSE]
  }

  out.derPItoRho <- derPItoRho(ind.vec = LONG1, th.rho.vec = LONG2)

  derPItoRho.times.PI.inv <- (1 / PI) * out.derPItoRho

  out.der.pi.xixj.to.tau.xi <- der.pi.xixj.to.tau.xi(ind.vec = LONG1, th.rho.vec = LONG2)

  out.der.pi.xixj.to.tau.xj <- der.pi.xixj.to.tau.xj(ind.vec = LONG1, th.rho.vec = LONG2)

  comb.idx <- combn(no.x, 2)


  if (!casewise) {
    if (nrow(Lx.free.index) > 0 | nrow(Phi.free.index) > 0) {
      prod.terms.for.derPLtoRho <- Data * derPItoRho.times.PI.inv
      xnew <- lapply(LONG1[c("index.pairs.extended")], function(y) {
        y[LONG1$index.thres.var1.of.pair != 0 &
          LONG1$index.thres.var2.of.pair != 0]
      })
      derPLtoRho <- tapply(prod.terms.for.derPLtoRho, xnew$index.pairs.extended, sum)

      derPLtoLamdaPhi <- apply(derRhotoLamdaAndPhi, 2, function(x) {
        matrix(derPLtoRho, nrow = 1) %*%
          matrix(x, ncol = 1)
      })
    }
    n.over.pi <- Data / PI

    out.derPLtoTau <- derPLtoTau(
      ind.vec = LONG1, n.over.pi = n.over.pi, no.x = no.x,
      comb.idx = comb.idx,
      der.pi.xixj.to.tau.xi = out.der.pi.xixj.to.tau.xi,
      der.pi.xixj.to.tau.xj = out.der.pi.xixj.to.tau.xj
    )

    if (nrow(Lx.free.index) > 0 | nrow(Phi.free.index) > 0) {
      res <- c(derPLtoLamdaPhi, out.derPLtoTau) # this is the order of lavaan
    } else {
      res <- out.derPLtoTau
    }
  } else { # for the casewise=TRUE option

    nobs <- nrow(Data)
    no.comb <- ncol(comb.idx)

    derPLtoRho.casewise <- matrix(NA, nrow = nobs, ncol = no.comb)
    derPLtoLamdaPhi.casewise <- matrix(NA,
      nrow = nobs,
      ncol = (nrow(Lx.free.index) + nrow(Phi.free.index))
    )

    out.derPLtoTau.casewise <- matrix(NA, nrow = nobs, ncol = length(all.thres))

    no.cat.x <- rep(NA, no.x)
    for (i in 1:no.x) {
      no.cat.x[i] <- sum(ParTable$lhs == names.x[i] & ParTable$op == "|") + 1
    }
    LONG3 <- vector("list", length = 4)
    names(LONG3) <- c("cat.idx.xi", "cat.idx.xj", "idx.xi", "idx.xj")
    LONG3$cat.idx.xi <-
      LONG1$index.thres.var1.of.pair[!(LONG1$index.thres.var1.of.pair == 0 |
        LONG1$index.thres.var2.of.pair == 0)]
    LONG3$cat.idx.xj <-
      LONG1$index.thres.var2.of.pair[!(LONG1$index.thres.var1.of.pair == 0 |
        LONG1$index.thres.var2.of.pair == 0)]
    LONG3$idx.xi <-
      LONG1$index.var1.of.pair[!(LONG1$index.thres.var1.of.pair == 0 |
        LONG1$index.thres.var2.of.pair == 0)]
    LONG3$idx.xj <-
      LONG1$index.var2.of.pair[!(LONG1$index.thres.var1.of.pair == 0 |
        LONG1$index.thres.var2.of.pair == 0)]


    no.bifreq <- length(PI)
    vec.bifreq.idx <- c(1:no.bifreq)

    for (i in 1:nobs) {
      obs <- Data[i, ]
      pair.xixj.obs <- matrix(as.numeric(obs[comb.idx]), ncol = 2, byrow = TRUE)
      n.xixj.ab <- rep(0, no.bifreq)

      if (nrow(Lx.free.index) > 0 | nrow(Phi.free.index) > 0) {
        derPLtoRho.casei <- rep(NA, no.comb)
      }

      for (j in 1:no.comb) {
        bifreq.idx <- vec.bifreq.idx[(
          LONG3$idx.xi == comb.idx[1, j] &
            LONG3$idx.xj == comb.idx[2, j] &
            LONG3$cat.idx.xi == pair.xixj.obs[j, 1] &
            LONG3$cat.idx.xj == pair.xixj.obs[j, 2])]

        if (nrow(Lx.free.index) > 0 | nrow(Phi.free.index) > 0) {
          derPLtoRho.casei[j] <- derPItoRho.times.PI.inv[bifreq.idx]
        }
        n.xixj.ab[bifreq.idx] <- 1
      }

      if (nrow(Lx.free.index) > 0 | nrow(Phi.free.index) > 0) {
        derPLtoRho.casewise[i, ] <- derPLtoRho.casei

        derPLtoLamdaPhi.casewise[i, ] <- apply(
          derRhotoLamdaAndPhi, 2,
          function(x) {
            matrix(derPLtoRho.casei, nrow = 1) %*%
              matrix(x, ncol = 1)
          }
        )
      }

      n.over.pi <- n.xixj.ab / PI

      out.derPLtoTau.casewise[i, ] <- derPLtoTau(
        ind.vec = LONG1,
        n.over.pi = n.over.pi,
        no.x = no.x,
        comb.idx = comb.idx,
        der.pi.xixj.to.tau.xi = out.der.pi.xixj.to.tau.xi,
        der.pi.xixj.to.tau.xj = out.der.pi.xixj.to.tau.xj
      )
    }

    if (nrow(Lx.free.index) > 0 | nrow(Phi.free.index) > 0) {
      res <- cbind(derPLtoLamdaPhi.casewise, out.derPLtoTau.casewise)
    } else {
      res <- out.derPLtoTau.casewise
    }
  }

  res # output
} # end of function




####### internal functions of GrPloglik
derPItoRho <- function(ind.vec, th.rho.vec) {
  prob.vec <- rep(NA, length(ind.vec$index.thres.var1.of.pair))

  prob.vec[ind.vec$index.thres.var1.of.pair == 0 |
    ind.vec$index.thres.var2.of.pair == 0 |
    ind.vec$last.thres.var1.of.pair |
    ind.vec$last.thres.var2.of.pair] <- 0

  prob.vec[is.na(prob.vec)] <- lavaan:::dbinorm(th.rho.vec$thres.var1.of.pair,
    th.rho.vec$thres.var2.of.pair,
    rho = th.rho.vec$rho.vector
  )

  den.term1 <- prob.vec[ind.vec$index.thres.var1.of.pair != 0 &
    ind.vec$index.thres.var2.of.pair != 0]

  den.term2 <- prob.vec[ind.vec$index.thres.var1.of.pair != 0 &
    !ind.vec$last.thres.var2.of.pair]

  den.term3 <- prob.vec[ind.vec$index.thres.var2.of.pair != 0 &
    !ind.vec$last.thres.var1.of.pair]

  den.term4 <- prob.vec[!ind.vec$last.thres.var1.of.pair &
    !ind.vec$last.thres.var2.of.pair]

  den.term1 - den.term2 - den.term3 + den.term4
}



#######
der.pi.xixj.to.tau.xi <- function(ind.vec, th.rho.vec) {
  xi <- lapply(ind.vec[c(
    "index.thres.var2.of.pair",
    "last.thres.var2.of.pair"
  )], function(y) {
    y[!(ind.vec$index.thres.var1.of.pair == 0 |
      ind.vec$last.thres.var1.of.pair)]
  })

  cum.prob.vec <- rep(NA, length(xi$index.thres.var2.of.pair))

  cum.prob.vec[xi$index.thres.var2.of.pair == 0] <- 0

  cum.prob.vec[xi$last.thres.var2.of.pair] <- 1

  denom <- (1 - th.rho.vec$rho.vector^2)^0.5

  cum.prob.vec[is.na(cum.prob.vec)] <- pnorm((th.rho.vec$thres.var2.of.pair -
    th.rho.vec$rho.vector * th.rho.vec$thres.var1.of.pair) / denom)

  den.prob.vec <- dnorm(th.rho.vec$thres.var1.for.dnorm.in.der.pi.to.tau.xi)

  den.prob.vec * (cum.prob.vec[xi$index.thres.var2.of.pair != 0] -
    cum.prob.vec[!xi$last.thres.var2.of.pair])
}


##############
der.pi.xixj.to.tau.xj <- function(ind.vec, th.rho.vec) {
  xj <- lapply(
    ind.vec[c(
      "index.thres.var1.of.pair",
      "last.thres.var1.of.pair"
    )],
    function(y) {
      y[!(ind.vec$index.thres.var2.of.pair == 0 |
        ind.vec$last.thres.var2.of.pair)]
    }
  )

  cum.prob.vec <- rep(NA, length(xj$index.thres.var1.of.pair))

  cum.prob.vec[xj$index.thres.var1.of.pair == 0] <- 0

  cum.prob.vec[xj$last.thres.var1.of.pair] <- 1

  denom <- (1 - th.rho.vec$rho.vector^2)^0.5

  cum.prob.vec[is.na(cum.prob.vec)] <- pnorm((th.rho.vec$thres.var1.of.pair -
    th.rho.vec$rho.vector * th.rho.vec$thres.var2.of.pair) / denom)

  den.prob.vec <- dnorm(th.rho.vec$thres.var2.for.dnorm.in.der.pi.to.tau.xj)

  den.prob.vec * (cum.prob.vec[xj$index.thres.var1.of.pair != 0] -
    cum.prob.vec[!xj$last.thres.var1.of.pair])
}

########



derPLtoTau <- function(ind.vec, n.over.pi, no.x, comb.idx,
                       der.pi.xixj.to.tau.xi, der.pi.xixj.to.tau.xj) {
  x3a <- lapply(ind.vec, function(y) {
    y[!(ind.vec$index.thres.var1.of.pair == 0 |
      ind.vec$index.thres.var2.of.pair == 0)]
  })

  diff.n.over.pi.to.xi <- n.over.pi[!x3a$last.thres.var1.of.pair] -
    n.over.pi[x3a$index.thres.var1.of.pair != 1]

  diff.n.over.pi.to.xj <- n.over.pi[!x3a$last.thres.var2.of.pair] -
    n.over.pi[x3a$index.thres.var2.of.pair != 1]

  terms.der.Lxixj.to.tau.xi <- diff.n.over.pi.to.xi * der.pi.xixj.to.tau.xi

  terms.der.Lxixj.to.tau.xj <- diff.n.over.pi.to.xj * der.pi.xixj.to.tau.xj

  x3b <- lapply(
    ind.vec[c(
      "index.pairs.extended",
      "index.thres.var1.of.pair"
    )],
    function(y) {
      y[!(ind.vec$index.thres.var1.of.pair == 0 |
        ind.vec$last.thres.var1.of.pair |
        ind.vec$index.thres.var2.of.pair == 0)]
    }
  )
  x4b <- lapply(
    ind.vec[c(
      "index.pairs.extended",
      "index.thres.var2.of.pair"
    )],
    function(y) {
      y[!(ind.vec$index.thres.var2.of.pair == 0 |
        ind.vec$last.thres.var2.of.pair |
        ind.vec$index.thres.var1.of.pair == 0)]
    }
  )

  # ind.pairs <- utils::combn(no.x, 2) # comb.idx

  der.Lxixj.to.tau.xi <- tapply(
    terms.der.Lxixj.to.tau.xi,
    list(
      x3b$index.pairs.extended,
      x3b$index.thres.var1.of.pair
    ), sum
  )

  der.Lxixj.to.tau.xj <- tapply(
    terms.der.Lxixj.to.tau.xj,
    list(
      x4b$index.pairs.extended,
      x4b$index.thres.var2.of.pair
    ), sum
  )

  split.der.Lxixj.to.tau.xi <- split(
    as.data.frame(der.Lxixj.to.tau.xi),
    comb.idx[1, ]
  )

  sums.der.Lxixj.to.tau.xi <- lapply(
    split.der.Lxixj.to.tau.xi,
    function(x) {
      y <- apply(x, 2, sum)
      y[!is.na(y)]
    }
  )

  split.der.Lxixj.to.tau.xj <- split(
    as.data.frame(der.Lxixj.to.tau.xj),
    comb.idx[2, ]
  )

  sums.der.Lxixj.to.tau.xj <- lapply(
    split.der.Lxixj.to.tau.xj,
    function(x) {
      y <- apply(x, 2, sum)
      y[!is.na(y)]
    }
  )

  c(
    sums.der.Lxixj.to.tau.xi[[1]],
    c(unlist(sums.der.Lxixj.to.tau.xi[2:(no.x - 1)]) +
      unlist(sums.der.Lxixj.to.tau.xj[1:(no.x - 2)])),
    sums.der.Lxixj.to.tau.xj[[no.x - 1]]
  )
}




#######################################################################


# Estimation of the Asymptotic Covariance matrix of the estimated parameters
# output : the estimated asymptotic variance-covariance matrix of the estimated
#          lamdas, phis, taus

# input :
# par.est : the parameter estimates given in the order: lamdas columnwise,
#           phis columnwise of the lower triangular of matrix Phi
# MODEL : an model object as should be provided in lavaan
#         make sure that it explicitly includes the variances and covariances of
#         latent variables and all thresholds for all variables (excluding
#         the redundant ones -Inf and Inf)
#         also write NA in front of the first indicator of each latent variable
#         to denote that it is free to be estimated
# Data : row data in matrix object
# GrPloglik.fun : it is the GrPloglik function


AsymVCOVLamdaPhiTau <- function(par.est, MODEL, DATA) {
  ParTable <- lavaanify(MODEL)
  names.x <- lavaanNames(ParTable, type = "ov")
  no.x <- length(names.x)

  GrPloglik.casewise <- GrPloglik(
    free.par = par.est, model = MODEL,
    Data = DATA, casewise = TRUE
  )

  no.par <- length(par.est)
  # for case i the elements of J matrix is saved in column i
  # the elements of the J matrix are listed columnwise in the ith column
  J.casewise <- apply(GrPloglik.casewise, 1, function(x) {
    outer(x, x)
  })
  J.mat <- matrix(rowSums(J.casewise), nrow = no.par, ncol = no.par)

  n.xixj.ab <- unlist(lavaan:::pairwiseTables(data = DATA, no.x = no.x)$pairTables)

  H.mat <- numericGradient(
    f = GrPloglik, t0 = par.est,
    model = MODEL, Data = n.xixj.ab, casewise = FALSE
  )

  minusHinv <- (-1) * solve(H.mat)

  minusHinv %*% J.mat %*% minusHinv
}


##### function to give J and H, H is the inverted second derivative *(-1)
J.InvMinusHinv.timesn <- function(par.est, MODEL, DATA) {
  ParTable <- lavaanify(MODEL)
  names.x <- lavaanNames(ParTable, type = "ov")
  no.x <- length(names.x)

  GrPloglik.casewise <- GrPloglik(
    free.par = par.est, model = MODEL,
    Data = DATA, casewise = TRUE
  )

  no.par <- length(par.est)
  # for case i the elements of J matrix is saved in column i
  # the elements of the J matrix are listed columnwise in the ith column
  J.casewise <- apply(GrPloglik.casewise, 1, function(x) {
    outer(x, x)
  })
  J.mat <- matrix(rowSums(J.casewise), nrow = no.par, ncol = no.par)

  n.xixj.ab <- unlist(lavaan:::pairwiseTables(data = DATA, no.x = no.x)$pairTables)

  H.mat <- numericGradient(
    f = GrPloglik, t0 = par.est,
    model = MODEL, Data = n.xixj.ab, casewise = FALSE
  )

  minusHinv <- (-1) * solve(H.mat)

  J.mat %*% minusHinv
}

##### function to give Hinv which is the inverted second derivative *(-1)
MinusHinv.timesn <- function(par.est, MODEL, DATA, GrPloglik.fun = GrPloglik) {
  ParTable <- lavaanify(MODEL)
  names.x <- lavaanNames(ParTable, type = "ov")
  no.x <- length(names.x)

  n.xixj.ab <- unlist(lavaan:::pairwiseTables(data = DATA, no.x = no.x)$pairTables)

  H.mat <- numericGradient(
    f = GrPloglik, t0 = par.est,
    model = MODEL, Data = n.xixj.ab, casewise = FALSE
  )

  (-1) * solve(H.mat)
}

# Commands for SRS and no weights ----------------------------------------------
library(lavaan)
library(mvtnorm)
library(fCopulae)
library(Matrix)
library(dr)
library(mnormt)
# library(raster)
# simulation set up
no.x <- 5
no.ksi <- 1
n.obs <- 1000
no.sim <- 10
# quantities to be saved in the simulations
# L_2 saves the values of the test statistics per simulation and two types of degrees of freedom
power <- matrix(NA, no.sim, 3)
L_2 <- matrix(NA, no.sim, 3)
Chi_2 <- rep(NA, no.sim)
df_eigen_omega <- matrix(NA, no.sim, (no.x * (no.x - 1) / 2 + no.x))
df_eigen_pearson <- matrix(NA, no.sim, (no.x * (no.x - 1) / 2 + no.x))
Chi_2_Rao <- rep(NA, no.sim)
Chi_2_Rao2 <- matrix(NA, no.sim, 2)
df_rank <- matrix(NA, no.sim, 2)
Sav.Omega.diag <- matrix(NA, (no.x * (no.x - 1) / 2 + no.x), no.sim)
Sav.P2 <- matrix(NA, (no.x * (no.x - 1) / 2 + no.x), no.sim)
Sav.P2.theta <- matrix(NA, (no.x * (no.x - 1) / 2 + no.x), no.sim)
Sav.P2.t <- matrix(NA, (no.x * (no.x - 1) / 2 + no.x), no.sim)
#
# the true loadings and the thresholds have been set equal to the lsat results
# one-facro
Lx.matrix <- matrix(data = c(0.8, 0.7, 0.47, 0.38, 0.34), nrow = no.x, ncol = no.ksi)
Phi.matrix <- matrix(c(1), no.ksi, no.ksi)
# Two factot model
# Lx.matrix <- matrix(data=c(0.7, 0.6, 0.5, 0.4, 0.3, 0.3, 0.4, 0.5, 0.6, 0.7), nrow=no.x,ncol=no.ksi )
# Phi.matrix <- matrix(c(1,0,0,1), no.ksi, no.ksi)
Lx.Phi.tLx <- Lx.matrix %*% Phi.matrix %*% t(Lx.matrix)
Theta.d.matrix <- diag(1 - diag(Lx.Phi.tLx))
Sigma <- Lx.Phi.tLx + Theta.d.matrix
cors <- Sigma[lower.tri(Sigma)]
free.tau.x <- c(-1.4325, -0.5504, -0.1332, -0.71599, -1.1263)
# TH <-  rep(free.tau.x , no.x)  # same thresholds for all binary variables
TH <- c(free.tau.x) # different thresholds for all binary variables
# th.idx  <- rep(1:no.x, each=length(free.tau.x ))
# one-factor
free.par.values <- c(0.6, 0.5, 0.47, 0.38, 0.34, TH)
# specify the model in lavaan syntax
model <- "
 ksi1 =~ NA*x1 + x2 + x3 + x4 +x5

 ksi1 ~~ 1*ksi1

 x1 | t1
 x2 | t1
 x3 | t1
 x4 | t1
 x5 | t1
"
# Generate data
no.ksi <- 1
# data.sim.array <- array(NA, dim=c(n.obs, no.x, no.sim))
#
for (l in (1:no.sim)) {
  data <- array(NA, dim = c(n.obs, no.x))
  # GENERATE DATA   # library(mnormt)  included in lavaan
  under.data.sim.x <- rmnorm(n.obs, mean = rep(0, no.x), varcov = Sigma)
  # binary variables
  for (i in (1:no.x)) {
    data[, i] <- (under.data.sim.x[, i] <= free.tau.x[i]) * 1 +
      (under.data.sim.x[, i] > free.tau.x[i]) * 2
  }
  #   data <-  (under.data.sim.x <= free.tau.x )*1 +
  #                (under.data.sim.x >  free.tau.x ) *2
  #  the data.sim.array will keep all the simulated data to be re-used in the future
  #   data.sim.array[,,l] <- data
  #
  data <- as.data.frame(data)
  data[, ] <- lapply(data, ordered)
  names(data) <- paste("x", 1:no.x, sep = "")
  #
  n.xixj.ab <- unlist(lavaan:::pairwiseTables(
    data = data,
    no.x = no.x
  )$pairTables)
  #
  # maximise the logL using maxLik
  start.values <- c(0.4, 0.5, 0.7, 0.6, 0.8, -0.5, -0.7, -0.6, -0.8, -0.4)
  # the maxLik uses the derivatives
  pml <- maxLik(
    logLik = Ploglik, grad = GrPloglik, hess = NULL,
    start = start.values, method = "BFGS",
    #              start=c( rep(0.5,5), rep(0,5)), method="BFGS",
    model = model, Data = n.xixj.ab, casewise = FALSE
  )

  pml$est # gives the vector with the estimates

  AsymVCOV <- AsymVCOVLamdaPhiTau(
    par.est = pml$est, model = model,
    Data = data, GrPloglik.fun = GrPloglik, J.mat
  )

  # AsymVCOV #it gives the asymptotic variance-covariance matrix of the estimates

  # compute the test elements
  par.est <- coef(pml)

  # compute unstandardised residuals (n.xixj.ab/n-pi.xixj.ab)
  # For that I will need LONG1, LONG2 and PI
  # we also have to recode the data to 0 and 1.
  data.r <- (data == 1) * 0 + (data == 2) * 1
  data.r <- as.data.frame(data.r)
  #
  Lx.matrix <- coef(pml)[1:(no.x * no.ksi)]
  len_Lx.matrix <- length(Lx.matrix)
  all.thres <- coef(pml)[(len_Lx.matrix + 1):(no.ksi * no.x + no.x)]
  th.idx <- seq(1:no.x)
  Phi.matrix <- matrix(1, no.ksi, no.ksi)
  Lx.Phi.tLx <- Lx.matrix %*% Phi.matrix %*% t(Lx.matrix)
  # the matrix below contains all the lower triangular rho's without the diagonal, elements are stored columnwise
  rho.xixj <- Lx.Phi.tLx[lower.tri(Lx.Phi.tLx)]
  LONG2 <- lavaan:::LongVecTH.Rho(
    no.x = no.x,
    all.thres = all.thres,
    index.var.of.thres = th.idx,
    rho.xixj = rho.xixj
  )

  LONG1 <- lavaan:::LongVecInd(
    no.x = no.x,
    all.thres = all.thres,
    index.var.of.thres = th.idx
  )
  #
  # PI computes the estimated bivariate probabilities for all combinations of items and categories
  PI <- lavaan:::pairwiseExpProbVec(ind.vec = LONG1, th.rho.vec = LONG2)
  # n.xixj.ab/n.obs contains the corresponding observed proportions
  # we need to select only the (1,1) responses for all combinations of items.
  # to do that we need to create a vector that contains the place of the elements to be selected
  to.no.pairs <- no.x * (no.x - 1) / 2 * 4
  # in PI the elements are given in the order of P00,P10,P01,P11 for each pair of items
  ind.pi <- seq(4, to.no.pairs, 4)
  P_11 <- PI[ind.pi]
  n_11 <- n.xixj.ab[ind.pi]
  # observed proportions of positive responses only
  p_11 <- n_11 / n.obs
  # compute the univariate estimated positive responses
  P_1 <- pnorm(all.thres, lower.tail = F)
  # observed univariate proportions of positive responses
  p_1 <- colMeans(data.r, na.rm = FALSE, dims = 1)
  #
  # compute trivariate estimated positive proportions for all combinations of items.
  tri.pairs <- combn(no.x, 3)
  no.tri.pairs <- ncol(tri.pairs)
  P_111 <- array(NA, no.tri.pairs)
  # P_111<-array(NA,no.tri.pairs,1)
  m <- 3
  S <- diag(m)
  for (i in 1:no.tri.pairs) {
    ind1.1 <- c(tri.pairs[1, i])
    ind1.2 <- c(tri.pairs[2, i])
    ind2.1 <- c(tri.pairs[1, i])
    ind2.2 <- c(tri.pairs[3, i])
    ind3.1 <- c(tri.pairs[2, i])
    ind3.2 <- c(tri.pairs[3, i])
    S[2, 1] <- S[1, 2] <- Lx.Phi.tLx[ind1.1, ind1.2]
    S[3, 1] <- S[1, 3] <- Lx.Phi.tLx[ind2.1, ind2.2]
    S[3, 2] <- S[2, 3] <- Lx.Phi.tLx[ind3.1, ind3.2]
    diag(S) <- 1
    S
    prob.temp <- pmvnorm(
      lower = c(all.thres[tri.pairs[1, i]], all.thres[tri.pairs[2, i]], all.thres[tri.pairs[3, i]]),
      upper = c(Inf, Inf, Inf),
      mean = rep(0, m), sigma = S,
      algorithm = Miwa()
    )
    P_111[i] <- prob.temp[[1]]
  }
  # compute four-way estimated positive proportions.
  four.pairs <- combn(no.x, 4)
  no.four.pairs <- ncol(four.pairs)
  P_1111 <- array(NA, no.four.pairs)
  # P_1111<-array(NA,no.four.pairs,1)
  m <- 4
  S <- diag(m)
  for (i in 1:no.four.pairs) {
    ind1.1 <- c(four.pairs[1, i])
    ind1.2 <- c(four.pairs[2, i])
    ind2.1 <- c(four.pairs[1, i])
    ind2.2 <- c(four.pairs[3, i])
    ind3.1 <- c(four.pairs[1, i])
    ind3.2 <- c(four.pairs[4, i])
    ind4.1 <- c(four.pairs[2, i])
    ind4.2 <- c(four.pairs[3, i])
    ind5.1 <- c(four.pairs[2, i])
    ind5.2 <- c(four.pairs[4, i])
    ind6.1 <- c(four.pairs[3, i])
    ind6.2 <- c(four.pairs[4, i])
    S[2, 1] <- S[1, 2] <- Lx.Phi.tLx[ind1.1, ind1.2]
    S[3, 1] <- S[1, 3] <- Lx.Phi.tLx[ind2.1, ind2.2]
    S[4, 1] <- S[1, 4] <- Lx.Phi.tLx[ind3.1, ind3.2]
    S[3, 2] <- S[2, 3] <- Lx.Phi.tLx[ind4.1, ind4.2]
    S[4, 2] <- S[2, 4] <- Lx.Phi.tLx[ind5.1, ind5.2]
    S[4, 3] <- S[3, 4] <- Lx.Phi.tLx[ind6.1, ind6.2]
    S
    prob.temp <- pmvnorm(
      lower = c(
        all.thres[four.pairs[1, i]], all.thres[four.pairs[2, i]],
        all.thres[four.pairs[3, i]], all.thres[four.pairs[4, i]]
      ),
      upper = c(Inf, Inf, Inf, Inf),
      mean = rep(0, m), sigma = S,
      algorithm = Miwa()
    )
    P_1111[i] <- prob.temp[[1]]
  }
  # unstandardised univariate residuals
  res.1 <- p_1 - P_1
  # unstandardised bivariate residuals
  res.11 <- p_11 - P_11
  # unstandardised residuals up to bivariate all put in one vector
  res.2 <- c(res.1, res.11)
  # Now I need to compute the components of the Omega matrix needed for computing
  # the M2 statistic
  #
  # the commands below compute the Sigma_2 matrix under the multinomial, SRS
  # the following three lines create the indices for the rows of the covariance matrix
  a1 <- seq(1:no.x)
  a2 <- combn(no.x, 2, simplify = FALSE)
  a <- c(a1, a2)
  #
  # the following three lines create the indices for the elements inside the covariance matrix
  a3 <- combn(no.x, 3, simplify = FALSE)
  if (no.x > 3) {
    a4 <- combn(no.x, 4, simplify = FALSE)
  } else {
    a4 <- NULL
  }
  b <- c(a1, a2, a3, a4)
  # the y vector contains the estimated probabilities given in the same order as the list b
  y <- c(P_1, P_11, P_111, P_1111)
  # k1 is the total number of univarite and bivariate probabilities
  k1 <- no.x * (no.x - 1) / 2 + no.x
  M <- matrix(0, k1, k1)
  for (i in 1:k1) {
    for (j in 1:k1) {
      index <- union(a[[i]], a[[j]])
      index <- sort(index)
      for (k in (1:length(b))) {
        if (all(index == b[[k]])) {
          g <- k
        }
      }
      M[i, j] <- y[g] - y[i] * y[j]
    }
  }
  # Compute Beta(theta)
  Beta <- Beta_design(no.x)
  PI_inv <- PI^{
    -1
  }
  y_1 <- biv_der(par.est, no.x, no.ksi)
  Delta_theta <- y_1$Delta_all
  Beta_theta <- matrix(NA, (no.x + no.x * no.ksi), (no.x + length(a2)))
  for (j in 1:(no.x + no.x * no.ksi)) {
    y1 <- PI_inv * Delta_theta[j] * Beta
    Beta_theta[j, ] <- colSums(y1)
  }
  # compute the L2 test statistic
  # the function below computes the derivatives of the logL
  free.par <- pml$est
  # der.pl<-GrPloglik(free.par, model, Data=n.xixj.ab, casewise=FALSE)
  #####

  # Compute the H(theta)^{-1} matrix
  H.mat <- numericGradient(
    f = GrPloglik, t0 = par.est,
    model = model, Data = n.xixj.ab, casewise = FALSE
  )
  minusHinv <- (-1) * solve(H.mat)

  # Omega matrix
  y_1$Delta
  I_s <- diag(no.x + length(a2))
  Omega <- (I_s - y_1$Delta %*% minusHinv %*% Beta_theta) %*% M %*% t((I_s - y_1$Delta %*% minusHinv %*% Beta_theta))
  Sav.Omega.diag[, l] <- diag(Omega)
  Omega_inv <- ginv(Omega)
  df_eigen_omega[l, ] <- eigen(Omega)$values
  # compute the rank of the matrix using the Matrix package
  d1 <- rankMatrix(Omega, tol = 0.0001)[1]
  d2 <- rankMatrix(Omega, tol = 0.0001, "useGrad")[1]
  df_rank[l, 1] <- d1
  df_rank[l, 2] <- d2
  #
  # the wald test statistic
  L_2[l, 1] <- n.obs * (t(res.2)) %*% Omega_inv %*% res.2
  L_2[l, 2] <- d1 # enter the d.f in the L_2 test computed by the two different methods
  L_2[l, 3] <- d2
  # the Pearson test statistic
  # the pi_2_theta is needed for the Pearson test
  pi_2 <- c(P_1, P_11) # estimated proportions
  pi_2t <- c(p_1, p_11) # observed proportions
  # Sav.P2[,l]<-pi_2 #estimated univariate and bivariate proportions
  # Sav.P2.t[,l]=(pi_2t-pi_2) #residuals
  Sav.P2.theta[, l] <- n.obs * (pi_2t - pi_2)^2 / pi_2
  # t4<-rowMeans(Sav.P2.theta) #sum all the components across simulations
  # sum(t4) #that should be the chi^2 Pearson test
  #
  pi_2_inv <- (pi_2)^{
    -1
  }
  pi_2_theta_inv <- diag(pi_2_inv)
  #
  pi_2_inv1 <- (pi_2)^{
    -0.5
  }
  pi_2_theta_inv1 <- diag(pi_2_inv1)
  Chi_2[l] <- n.obs * ((t(res.2)) %*% pi_2_theta_inv %*% res.2)
  #
  # distribution of the Chi_2 test statistic
  # compute the eigen values of D_{-0.5}OmegaD_{-0.5}
  mat_dist <- (pi_2_theta_inv1 %*% Omega %*% pi_2_theta_inv1)
  eig <- eigen(mat_dist)$values
  #
  # Apply the first and second order correction.
  # The sum of the weighted independent chi-square variables of 1 df each should go from 1 to s-t
  # So that all the computations below we only use the s-t eigen values
  # mean_eigen<-mean(eig[1:(length(pi_2)-2*no.x)])
  # if we use all s values the
  mean_eigen <- mean(eigen(mat_dist)$values)
  #
  # First order Rao-Scott Chi square test
  Chi_2_Rao[l] <- Chi_2[l] / mean_eigen
  #
  # Second order Rao-Scott Chi square and Wald tests
  # compute the coefficient of variation of the eigen values
  # use only s-t
  # sd_eigen<-sd(eig[1:(length(pi_2)-2*no.x)])
  # use all s values
  sd_eigen <- sd(eig)
  # coefficient of variation
  cv_eigen <- sd_eigen / mean_eigen
  #
  Chi_2_Rao2[l, 1] <- Chi_2_Rao[l] / (1 + cv_eigen^2)
  # save also the degrees of freedom
  # Chi_2_Rao2[l,2]<-(sum(eig[1:(length(pi_2)-2*no.x)])^2)/(sum(eig[1:(length(pi_2)-2*no.x)]^2))
  # or when all s values are used
  Chi_2_Rao2[l, 2] <- (no.x * (no.x - 1) / 2 + no.x) / (1 + cv_eigen^2)
  #
  # estimate the power of the test
  power[l, 1] <- L_2[l, 1] > qchisq(0.99, df = mean(L_2[l, 3]), ncp = 0, lower.tail = TRUE, log.p = FALSE)
  summary(x1)
  power[l, 2] <- L_2[l, 1] > qchisq(0.95, df = mean(L_2[l, 3]), ncp = 0, lower.tail = TRUE, log.p = FALSE)
  summary(x1)
  power[l, 3] <- L_2[l, 1] > qchisq(0.90, df = mean(L_2[l, 3]), ncp = 0, lower.tail = TRUE, log.p = FALSE)
  summary(x1)

  cat(l, "\n")
} # the end of the simulation (for command)

# Check some results for the Chi^2
# t1<-rowMeans(Sav.Omega.diag)
# t2<-rowMeans(Sav.P2)
# t1/t2

# plots
mean(L_2[, 1])
var(L_2[, 1])
sd(L_2[, 1])
hist(L_2[, 1], prob = TRUE, main = "Histogram of Wald test, SRS", xlab = "Wald test")
lines(density(L_2[, 1]), lwd = 3, col = "orange")
curve(dchisq(x, df = 5), lwd = 3, col = "green", add = TRUE)
curve(dchisq(x, df = 6), col = "green", add = TRUE)
qqplot(rchisq(1000, df = 5), L_2[, 1],
  main = "L_2, df=5 Q-Q Plot",
  ylab = "Sample Quantiles"
)
qqplot(rchisq(1000, df = 6), L_2[, 1],
  main = "L_2, df=6 Q-Q Plot",
  ylab = "Sample Quantiles"
)
# Test the hypothesis that L_2 comes from a chi-square distribution with df=x
ks.test(L_2[, 1], pchisq,
  df = mean(L_2[, 3]), alternative = c("two.sided"),
  exact = NULL
)
#
# obtain empirical p-values
x1 <- L_2[, 1] > qchisq(0.99, df = mean(L_2[, 3]), ncp = 0, lower.tail = TRUE, log.p = FALSE)
summary(x1)
x1 <- L_2[, 1] > qchisq(0.95, df = mean(L_2[, 3]), ncp = 0, lower.tail = TRUE, log.p = FALSE)
summary(x1)
x1 <- L_2[, 1] > qchisq(0.90, df = mean(L_2[, 3]), ncp = 0, lower.tail = TRUE, log.p = FALSE)
summary(x1)
#
#
# Pearson chi-square
mean(Chi_2)
var(Chi_2)
hist(Chi_2, prob = TRUE)
lines(density(Chi_2[]), col = "orange")
curve(dchisq(x, df = 3), col = "green", add = TRUE)
curve(dchisq(x, df = 5), col = "green", add = TRUE)
# test statistics with first order correctioin
mean(Chi_2_Rao)
var(Chi_2_Rao)
hist(Chi_2_Rao, prob = TRUE)
lines(density(Chi_2_Rao[]), col = "orange")
curve(dchisq(x, df = 3), col = "green", add = TRUE)
curve(dchisq(x, df = 4), col = "green", add = TRUE)
curve(dchisq(x, df = 5), col = "green", add = TRUE)
#
# test statistics with second order correctioin
mean(Chi_2_Rao2[, 1])
var(Chi_2_Rao2[, 1])
# degrees of freedom
mean(Chi_2_Rao2[, 2])
sd(Chi_2_Rao2[, 2])
#
hist(Chi_2_Rao2[, 1], prob = TRUE)
lines(density(Chi_2_Rao2[, 1]), col = "orange")
curve(dchisq(x, df = 1.3), col = "green", add = TRUE)
curve(dchisq(x, df = 2), col = "green", add = TRUE)
curve(dchisq(x, df = 3), col = "red", add = TRUE)
# obtain empirical p-values
x1 <- Chi_2_Rao2[, 1] > qchisq(0.99, df = mean(Chi_2_Rao2[, 2]), ncp = 0, lower.tail = TRUE, log.p = FALSE)
summary(x1)
x1 <- Chi_2_Rao2[, 1] > qchisq(0.95, df = mean(Chi_2_Rao2[, 2]), ncp = 0, lower.tail = TRUE, log.p = FALSE)
summary(x1)
x1 <- Chi_2_Rao2[, 1] > qchisq(0.90, df = mean(Chi_2_Rao2[, 2]), ncp = 0, lower.tail = TRUE, log.p = FALSE)
summary(x1)
#

hist(Chi_2[], prob = TRUE)
x <- rchisq(1000, 3)
hist(x, prob = TRUE)
curve(dchisq(x, df = 3), col = "green", add = TRUE)
dr.pvalue(df_eigen_pearson, Chi_2)

L2 <- princomp(covmat = Omega)
screeplot(L2)
