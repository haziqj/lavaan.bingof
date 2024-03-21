# Legacy code written by Myrsini Katsikatsou. Used for testing purposes.

# ------------------------------------------------------------------------------
#the derivatives of the UNIVARIATE & BIVARIATE probabilities (for all pairs
#and all response patterns) with respect to model parameter vector theta
#
#The computation of derModelBivProbToTheta is split into three parts:
# der.pi.xixj.ab.to.theta = der.pi.xixj.ab.to.Rho %*% der.Rho.to.Theta +
#                           der.pi.xixj.ab.to.tau.xi %*% der.Tau.to.Theta +
#                           der.pi.xixj.ab.to.tau.xj %*% der.Tau.to.Theta +

# NOTE: in lavaan 0.6-18.2002 the names of the lavcache have been lowercased!
# Convert $long -> $long

derModelUnivBivProbToTheta <- function(nvar = nvar,
                                         TH = TH,
                                     th.idx = th.idx,
                                  Sigma.hat = Sigmahat,
                                   lavcache = lavcache,
                                   lavmodel = lavmodel){

  idxTH_var1 <- lavcache$long$index.thres.var1.of.pair
  idxTH_var2 <- lavcache$long$index.thres.var2.of.pair
  idxLastTH_var1 <- lavcache$long$last.thres.var1.of.pair
  idxLastTH_var2 <- lavcache$long$last.thres.var2.of.pair

  idx_pairs_ext <- lavcache$long$index.pairs.extended
  cdtn <- (idxTH_var1!=0) & (idxTH_var2!=0)
  idx_pairs_S <- idx_pairs_ext[cdtn]
  idxVar1_S   <- lavcache$long$index.var1.of.pair[cdtn]
  idxVar2_S   <- lavcache$long$index.var2.of.pair[cdtn]
  idxCat_var1 <- idxTH_var1[cdtn]
  idxCat_var2 <- idxTH_var2[cdtn]

  cors <- Sigma.hat[lower.tri(Sigma.hat)]
  th.rho.vec <- LongVecTH.Rho(no.x = nvar,
                                  all.thres = TH,
                         index.var.of.thres = th.idx,
                                   rho.xixj = cors)

  #Delta below gives the derivatives of standardised thresholds and polychoric
  #correlations (rows of Delta in this order) with respect to model parameter vector
  #theta. In lavaan in a factor analysis model the order of the individual
  #parameters is: loadings, unstandardized thresholds, factor correlations
  Delta <- computeDelta(lavmodel = lavmodel)[[1]]
  noTH <- length(TH)
  derRhoToTheta <- Delta[-c(1:noTH), ]
  ##derRhoToTheta provides the derivatives of rhoij wrt theta
  derTauToTheta <- Delta[c(1:noTH), ]
  ##derivative of Tau Unstandardised To Theta


  ########### SECTION 1: der.pi.xixj.ab.to.Rho %*% der.Rho.to.Theta
  #######################################################################
  prob.vec <- rep(NA, length(idxTH_var1))

  prob.vec[idxTH_var1==0 | idxTH_var2==0 |
           idxLastTH_var1 | idxLastTH_var2] <- 0

  prob.vec[is.na(prob.vec)] <- dbinorm(th.rho.vec$thres.var1.of.pair,
                                       th.rho.vec$thres.var2.of.pair,
                                       rho = th.rho.vec$rho.vector)

  den.term1 <- prob.vec[idxTH_var1!=0 & idxTH_var2!=0]

  den.term2 <- prob.vec[ idxTH_var1!=0 & !idxLastTH_var2]

  den.term3 <- prob.vec[ idxTH_var2!=0 & !idxLastTH_var1]

  den.term4 <- prob.vec[!idxLastTH_var1 & !idxLastTH_var2]

  der.pi.xixj.ab.to.rho.xixj <- den.term1 - den.term2 - den.term3 + den.term4


  #Note that each row of derRhoToTheta corresponds to one rhoij. Thus it
  #needs to be repeated as many times as the cells of the bivariate frequency
  #table of the corresponding pair of variables (xi,xj) so that it can then be
  #mulitplied with der.pi.xixj.to.rho.xixj.
  #That's why we do the following:
  der.pi.xixj.ab.to.theta_part1 <- derRhoToTheta[idx_pairs_S,] *
                                   der.pi.xixj.ab.to.rho.xixj
  ##der.pi.xixj.ab.to.theta_part1 has as many rows as the bivariate probabilities
  ##pi.xixj.ab and as many columns as the dimension of theta.
  ##Note that in pi.xixj.ab, index "a" (index of category for var 1)
  #runs the fastest, index "b" (index of the vategory of variable 2 of the pair)
  #is the next fastest, next it is "j" (index of variable 2 of the pair), and
  #last it is "i" (index of variable 1 of the pair)


  ################### SECTION 2: der.pi.xixj.ab.to.tau.xi %*% der.Tau.to.Theta
  ############################################################################
  # to compute der.pi.xixj.1b.to.tau.xi
  #note that der.pi.xixj.to.tau.xi include only those pi.xixj.ab where the
  #category "a" corresponds to a threshold that is free to be estimated.
  #In the case of binary data it means the derivatives of pi.xixj.1b.
  #The derivatives of pi.xixj.2b are (-1)*der.pi.xixj.1b.
  xi <- lapply( lavcache$long[c("index.thres.var2.of.pair",
                                            "last.thres.var2.of.pair")],
                function(y){ y[!(idxTH_var1==0 | idxLastTH_var1)] } )

  cum.prob.vec <- rep(NA, length(xi$index.thres.var2.of.pair) )
  cum.prob.vec[xi$index.thres.var2.of.pair==0] <- 0
  cum.prob.vec[xi$last.thres.var2.of.pair] <- 1
  denom <- sqrt(1-(th.rho.vec$rho.vector*th.rho.vec$rho.vector))
  cum.prob.vec[is.na(cum.prob.vec)] <-
      pnorm( (th.rho.vec$thres.var2.of.pair -
              th.rho.vec$rho.vector* th.rho.vec$thres.var1.of.pair) /
             denom)
  den.prob.vec <- dnorm(th.rho.vec$thres.var1.for.dnorm.in.der.pi.to.tau.xi)
  der.pi.xixj.1b.to.tau.xi <-  den.prob.vec *
                            (cum.prob.vec[ xi$index.thres.var2.of.pair!=0] -
                             cum.prob.vec[!xi$last.thres.var2.of.pair] )

  #Since the derivatives of pi.xixj.2b are (-1)*der.pi.xixj.1b ,
  #the complete vector der.pi.xixj.ab.to.tau.xi is
  idc_der.pi.xixj.1b.to.tau.xi <- rep(TRUE, length(idx_pairs_S))
  idc_der.pi.xixj.1b.to.tau.xi[idxCat_var1==2] <- FALSE
  der.pi.xixj.ab.to.tau.xi <- rep(NA, length(idx_pairs_S))
  der.pi.xixj.ab.to.tau.xi[idc_der.pi.xixj.1b.to.tau.xi] <-
                                                    der.pi.xixj.1b.to.tau.xi
  der.pi.xixj.ab.to.tau.xi[!idc_der.pi.xixj.1b.to.tau.xi] <-
                                                (-1)*der.pi.xixj.1b.to.tau.xi

  ##Note that for binary data each row of derTauToTheta corresponds to one tau.i
  ##where i the vairable index, while der.pi.xixj.ab.to.tau.xi is a vector of
  ##length equal to as many as bivariate variables pi.xixj.ab. Thus each row of
  ##derTauToTheta needs to be repeated as many times as the cells of the bivariate
  ##frequency table of the corresponding pair of variables (xi,xj) (4 cells for
  ##binary data). Then, it can then be mulitplied with der.pi.xixj.ab.to.tau.xi.
  #That's why we do the following:
  der.pi.xixj.ab.to.theta_part2 <- derTauToTheta[idxVar1_S,] *
                                   der.pi.xixj.ab.to.tau.xi


  ################### SECTION 3: der.pi.xixj.ab.to.tau.xj %*% der.Tau.to.Theta
  ##########################################################################
  # to compute der.pi.xixj.to.tau.xj
  xj <- lapply(lavcache$long[c("index.thres.var1.of.pair",
                                           "last.thres.var1.of.pair")],
                function(y){ y[!(idxTH_var2==0 | idxLastTH_var2)] } )

  cum.prob.vec <- rep(NA, length(xj$index.thres.var1.of.pair) )
  cum.prob.vec[xj$index.thres.var1.of.pair==0] <- 0
  cum.prob.vec[xj$last.thres.var1.of.pair] <- 1
  #denom computed above
  #denom <- sqrt(1-(th.rho.vec$rho.vector * th.rho.vec$rho.vector))
  cum.prob.vec[is.na(cum.prob.vec)] <-
      pnorm( (th.rho.vec$thres.var1.of.pair -
              th.rho.vec$rho.vector* th.rho.vec$thres.var2.of.pair) /
              denom)
  den.prob.vec <- dnorm(th.rho.vec$thres.var2.for.dnorm.in.der.pi.to.tau.xj)
  der.pi.xixj.a1.to.tau.xj <-  den.prob.vec *
                            (cum.prob.vec[ xj$index.thres.var1.of.pair!=0] -
                             cum.prob.vec[!xj$last.thres.var1.of.pair] )
  #note that der.pi.xixj.to.tau.xj include only those pi.xixj.ab where the
  #category "b" corresponds to a threshold that is free to be estimated.
  #In the case of binary data only the derivative of pi.xixj.a1 is computed.
  #We use the fact that der.pi.xixj.a2 = - der.pi.xixj.a1.to.tau.xj to
  #complete the vector  der.pi.xixj.ab.to.tau.xj.
  idc_der.pi.xixj.a1.to.tau.xj <- rep(TRUE, length(idx_pairs_S))
  idc_der.pi.xixj.a1.to.tau.xj[idxCat_var2==2] <- FALSE
  der.pi.xixj.ab.to.tau.xj <- rep(NA, length(idx_pairs_S))
  der.pi.xixj.ab.to.tau.xj[idc_der.pi.xixj.a1.to.tau.xj] <-
                                             der.pi.xixj.a1.to.tau.xj
  der.pi.xixj.ab.to.tau.xj[!idc_der.pi.xixj.a1.to.tau.xj] <-
                                        (-1)*der.pi.xixj.a1.to.tau.xj

  ##As said, for binary data, each row of derTauToTheta corresponds to one tau.i
  ##where i the variable index, while der.pi.xixj.ab.to.tau.xj is a vector of
  ##length equal to as many as bivariate variables pi.xixj.ab. Thus each row of
  ##derTauToTheta needs to be repeated as many times as the cells of the bivariate
  ##frequency table of the corresponding pair of variables (xi,xj) (4 cells for
  ##binary data). Then, it can then be mulitplied with der.pi.xixj.ab.to.tau.xj.
  #That's why we do the following:
  der.pi.xixj.ab.to.theta_part3 <- derTauToTheta[idxVar2_S,] *
                                   der.pi.xixj.ab.to.tau.xj


  ################### SECTION 4: der.pi.xixj.ab.to.theta  and
  ###################            der.pi.xixj.22.to.theta (positive bivariate probs)
  ##########################################################################
  der.pi.xixj.ab.to.theta_F <- der.pi.xixj.ab.to.theta_part1 +
                               der.pi.xixj.ab.to.theta_part2 +
                               der.pi.xixj.ab.to.theta_part3

  der.pi.xixj.22.to.theta_F <- der.pi.xixj.ab.to.theta_F[(idxCat_var1==2 &
                                                          idxCat_var2==2  ), ]
  ################### SECTION 5: der.pi.xi.1.to.theta
  ##########################################################################
  ##Note that for our problem we need ONLY the derivatives of the POSITIVE
  ##univariate probabilities w.r.t to theta and NOT of all all univariate
  ##probabilities. The univariate probabilities are only functions of the
  ##standardised thresholds.
  ##The command below is correct only in the case of binary data.
  derUnivProb1wrtTau <- (-1)*dnorm(TH)
  derUnivProb1wrtTheta <- derTauToTheta * derUnivProb1wrtTau

  # function's output
  list(derBivProbwrtTheta = der.pi.xixj.ab.to.theta_F,
       derBivProb11wrtTheta = der.pi.xixj.22.to.theta_F,
       derUnivProb1wrtTheta = derUnivProb1wrtTheta )
}
############################################################


######### COMPUTE THE POSITIVE UNIVARIATE, BIVARIATE, TRIVARIATE, TETRAVARIATE
######### MODEL-IMPLIED PROBABILITIES as well as
######### THE POSITIVE UNIVARIATE and BIVARIATE PROPORTIONS AND RESIDUALS

positiveUniBivTrivTetravModelProb_ObsUnivBivPropResid <- function(
                                      lavobject = lavobject,
                                        lavdata = lavdata,
                                             TH = TH,
                                         th.idx = th.idx,
                                           nvar = nvar,
                                       Sigmahat = Sigmahat,
                                        Meanhat = Meanhat) {

 #### SECTION 1 : Compute the univariate quantities #######################
 #observed univariate frequencies
 univfreq_table <- lav_tables_oneway (lavobject = lavobject,
                                                 lavdata = lavdata)
 prop_1 <- univfreq_table$obs.prop[univfreq_table$rhs==1]
 ##Note that in lavaan has recoded the data from 0 and 1 to 1 and 2 so the
 ##positive answer is coded as 2 rather than 1 that is in the theory.

 univModelProb_PI <- univariateExpProbVec(TH = TH, th.idx = th.idx)
 modelProb_1 <- univModelProb_PI[univfreq_table$rhs==1]


 #unstandardised univariate residuals
 residuals_1 <- prop_1 - modelProb_1


 #### SECTION 2 : Compute the bivariate quantities #######################
 #the observed bivariate frequencies and proportions
 bivfreq_table <- lav_tables_pairwise_freq_cell(lavdata = lavdata,
                                                  as.data.frame. = TRUE)
 prop_11 <- (bivfreq_table$obs.freq / bivfreq_table$nobs)[bivfreq_table$row==2 &                                                               bivfreq_table$col==2]
 ##note that in lavaan the positive bivariate response pattern is coded as (2,2)

 #the model implied positive bivariate probabilities
 bivModelProb_PI <-
          unlist(lav_tables_pairwise_model_pi(lavobject = lavobject) )
 ##the function above gives the same result as
 ## pairwiseExpProbVec(ind.vec = LONG1, th.rho.vec=LONG2)
 modelProb_11 <- bivModelProb_PI[bivfreq_table$row==2 & bivfreq_table$col==2]

 #unstandardised bivariate probability residuals
 residuals_11 <- prop_11 - modelProb_11


 #### SECTION 3 : Compute the positive trivariate model probabilities, i.e.
 ####             Pr(yi=1, yj=1, yk=1) #######################
 triGroups <- combn(nvar,3)
 no_triGroups <- ncol (triGroups)
 modelProb_111 <- sapply(1:no_triGroups, function(i) {
                         idx_group <- triGroups[,i]
                        Sigmahat.r <- Sigmahat[idx_group, idx_group, drop=FALSE]
                              TH.r <- TH[idx_group]
                         Meanhat.r <- Meanhat[idx_group]
                  mnormt::sadmvn(lower=TH.r, upper=c(Inf,Inf,Inf),
                                 mean=Meanhat.r, varcov=Sigmahat.r) })

 ##to check if modelProb_111 computed correctly, for a simulated data, compute
 ##the observed trivariate positive proportions. The latter should be close to
 ##modelProb_111 and to obtain them do the following:
 ## obsData <- lavdata@X[[1]]
 ## triv_prop <- sapply(1:no_triGroups, function(i) {
 ##                       idx_group <- triGroups[,i]
 ##                       trivData <- obsData[,idx_group]
 ##                       table(trivData[,1], trivData[,2], trivData[,3])  })
 ## prop_111 <- triv_prop[2^3,] /nsize
 ## For this example prop_111 is very close to modelProb_111.


 #### SECTION 4: Compute the positive tetra-variate model probabilities, i.e.
 ####             Pr(yi=1, yj=1, yk=1, yl=1) #######################
 tetrGroups <- combn(nvar,4)
 no_tetrGroups <- ncol(tetrGroups)
 modelProb_1111 <- sapply(1:no_tetrGroups, function(i) {
                        idx_group <- tetrGroups[,i]
                       Sigmahat.r <- Sigmahat[idx_group, idx_group, drop=FALSE]
                             TH.r <- TH[idx_group]
                        Meanhat.r <- Meanhat[idx_group]
                 mnormt::sadmvn(lower=TH.r, upper=c(Inf,Inf,Inf, Inf),
                                mean=Meanhat.r, varcov=Sigmahat.r)     })

 ##to check if modelProb_1111 computed correctly, for a simulated data, compute
 ##the observed tetra-variate positive proportions. The latter should be close to
 ##modelProb_1111 and to obtain them do the following:
 ## obsData <- lavdata@X[[1]]
 ## tetra_prop <- sapply(1:no_tetrGroups, function(i) {
 ##                       idx_group <- tetrGroups[,i]
 ##                       tetrData <- obsData[,idx_group]
 ##                       table(tetrData[,1], tetrData[,2],
 ##                             tetrData[,3], tetrData[,4])  })
 ## prop_1111 <- tetra_prop[2^4,] /nsize
 ## For this example prop_1111 is very close to modelProb_1111.


 list(prop_1    = prop_1,
      prop_11   = prop_11,

      modelProb_1    = modelProb_1,
      modelProb_11   = modelProb_11,
      modelProb_111  = modelProb_111,
      modelProb_1111 = modelProb_1111,

      ALLbivModelProb = bivModelProb_PI,

      residuals_1    = residuals_1,
      residuals_11   = residuals_11,

      UnstdResUniv1Biv11 = c(residuals_1, residuals_11) )
}
####################################################################


######### SIGMA_2  ##############################
#Cov_pXi1_pXiXj11 is computed in 5 steps: first Var(pXi1) and then Cov(pXi1, pXj1)
#to fill in the upper left part of the matrix; Cov(pXi1, pXiXj11) to fill in
#the lower left part of the matrix ; Var(pXiXj11) and then Cov(pXiXj11, pXkXl11)
#to fill in the lower right part of the matrix.
#Note that first we compute ALL the elements of the lower tringular of the matrix
#and then we use them to provide the complete matrix.

Cov_pXi1_pXiXj11 <- function(modelProb_1 = modelProb_1,
                            modelProb_11 = modelProb_11,
                           modelProb_111 = modelProb_111,
                          modelProb_1111 = modelProb_1111) {

  cov_dim <- length( c(modelProb_1, modelProb_11))
  Cov_pXi1_pXiXj11 <- matrix(NA, nrow=cov_dim , ncol=cov_dim)

  no_univ <-  length(modelProb_1)
  no_biv  <- length(modelProb_11)
  idx_univ <- 1:no_univ
  idx_biv  <- (no_univ+1):(no_univ+no_biv)
  idx_pairs_tbl <- combn(no_univ, 2)
  no_pairs <- ncol(idx_pairs_tbl)
  idx_triples_tbl <- combn(no_univ,3)
  idx_tetrads_tbl <- combn(no_univ,4)
  #add names to the elements of modelProb_11, modelProb_111 and modelProb_1111;
  #will be used later
  names(modelProb_11)  <- apply(idx_pairs_tbl, 2, paste, collapse="")
  names(modelProb_111)  <- apply(idx_triples_tbl, 2, paste, collapse="")
  names(modelProb_1111) <- apply(idx_tetrads_tbl, 2, paste, collapse="")


  #compute Var(pXi1) and allocate them at the right place in Cov_pXi1_pXiXj11
  diag(Cov_pXi1_pXiXj11[ idx_univ, idx_univ]) <- modelProb_1 * (1-modelProb_1)

  #compute the Cov(pXi1, pXj1) and allocate them  at the right place in the lower
  #trigular part of Cov_pXi1_pXiXj11
  cov_pXi1_pXj1 <- modelProb_11 - modelProb_1[idx_pairs_tbl[1,]] *
                                  modelProb_1[idx_pairs_tbl[2,]]
  tmp <- Cov_pXi1_pXiXj11[ idx_univ, idx_univ]
  tmp[lower.tri(tmp)] <- cov_pXi1_pXj1
  Cov_pXi1_pXiXj11[ idx_univ, idx_univ] <- tmp


  #compute the Cov(pXi1, pXjXk11)
  ##If i=j then Cov(pXi1, pXiXk11)= pXiXk11    - pXi1 * pXk1
  ##Otherwise,  Cov(pXi1, pXjXk11)= pXiXjXk111 - pXi1 * pXjXk11
  #To compute term1 in the above two equations we need to keep track of the
  #indices i,j,k and thus we do the following:
  idx_cov_pXi1_pXjXk11 <- cbind( rep(idx_univ, each=no_pairs),
                                 rep(idx_pairs_tbl[1,], times=no_univ),
                                 rep(idx_pairs_tbl[2,], times=no_univ) )
  idx_cov_pXi1_pXjXk11_Char <- apply(idx_cov_pXi1_pXjXk11, 1, function(x) {
                                paste( unique(sort(x)), collapse="") })
  BivTrivProb <- c(modelProb_11, modelProb_111)
  Cov_pXi1_pXiXj11[idx_biv, idx_univ] <-
                       BivTrivProb[idx_cov_pXi1_pXjXk11_Char]  -
                       ( rep(modelProb_1,   each = no_pairs)*
                         rep(modelProb_11, times = no_univ)   )

  #compute Var(pXiXj11) and allocate them in the Cov matrix
  diag(Cov_pXi1_pXiXj11[idx_biv, idx_biv]) <- modelProb_11 * (1-modelProb_11)

  #compute Cov(pXiXj11, pXlXk11)
  ##If i=l then Cov(pXiXj11, pXiXk11) = pXiXjXk111    - pXiXj11 * pXiXk11
  ##Otherwise   Cov(pXiXj11, pXlXk11) = pXiXjXlXk1111 - pXiXj11 * pXlXk11
  #To compute the first term of the above equations correctly we need to keep track
  #of the indices i,j,l,k. Thus we do the following:
  idx_cov_pXiXj11_pXlXk11_list <- sapply(1:(no_pairs-1), function(i){
                        cbind( rep(idx_pairs_tbl[1,i], times=(no_pairs-i) ),
                               rep(idx_pairs_tbl[2,i], times=(no_pairs-i) ),
                               t(idx_pairs_tbl[,-(1:i)]) )  } )
  idx_cov_pXiXj11_pXlXk11 <- Reduce(rbind, idx_cov_pXiXj11_pXlXk11_list)
  idx_cov_pXiXj11_pXlXk11_Char <- apply(idx_cov_pXiXj11_pXlXk11, 1,
                        function(x) { paste( unique(sort(x)), collapse="") })
  TrivTetraProb <- c(modelProb_111, modelProb_1111)
  term1_cov_pXiXj11_pXlXk11 <- TrivTetraProb[idx_cov_pXiXj11_pXlXk11_Char]

  #the second term of the covariances
  term2_cov_pXiXj11_pXlXk11 <- unlist( sapply(1:(no_pairs-1), function(i) {                                  modelProb_11[i] * modelProb_11[(i+1):no_pairs]  } ) )

  #do the subtraction between the two terms and allocate
  part_Cov <- Cov_pXi1_pXiXj11[idx_biv, idx_biv]
  part_Cov[lower.tri(part_Cov)] <- term1_cov_pXiXj11_pXlXk11 -                                                                   term2_cov_pXiXj11_pXlXk11
  Cov_pXi1_pXiXj11[idx_biv, idx_biv] <- part_Cov

  #so far the lower tringular of the matrix is computed
  #produce the complete matrix
  lower_tri_elements <- Cov_pXi1_pXiXj11[lower.tri(Cov_pXi1_pXiXj11)]
  Cov_pXi1_pXiXj11 <- t(Cov_pXi1_pXiXj11)
  Cov_pXi1_pXiXj11[lower.tri(Cov_pXi1_pXiXj11)] <- lower_tri_elements

  Cov_pXi1_pXiXj11
}
##################################################################





# Beta matrix  - applies only to binary data
# Beta_mat_design <- function(nvar) {
#
#   no_biv  <- nvar * (nvar-1) /2
#   idx_univ <- 1:nvar
#   idx_biv <- (nvar+1):(nvar+no_biv)
#   idx_cat_var1 <- rep(c(0,1,0,1), times=no_biv)
#   idx_cat_var2 <- rep(c(0,0,1,1), times=no_biv)
#   idx_pairs <- combn(nvar,2)
#   idx_var1 <- rep( idx_pairs[1,], each=4)
#   idx_var2 <- rep( idx_pairs[2,], each=4)
#
#   nrow_B <- 4*no_biv
#   idx_row <- 1:nrow_B
#   B_mat <- matrix(0, nrow=nrow_B, ncol=(nvar + no_biv) )
#
#   #determine the pxixj11
#   diag( B_mat[ (idx_cat_var1==1)&(idx_cat_var2==1) , idx_biv]) <- 1
#
#   #determines the pxixj10
#   cdtn10 <- (idx_cat_var1==1)&(idx_cat_var2==0)
#   idc_row10 <- idx_row[cdtn10]
#   B_mat[ cbind(idc_row10, idx_univ[idx_var1[cdtn10]]) ] <- 1
#   diag( B_mat[idc_row10, idx_biv] ) <- (-1)
#
#   #determines the pxixj01
#   cdtn01 <- (idx_cat_var1==0)&(idx_cat_var2==1)
#   idc_row01 <- idx_row[cdtn01]
#   B_mat[ cbind(idc_row01, idx_univ[idx_var2[cdtn01]]) ] <- 1
#   diag( B_mat[idc_row01, idx_biv] ) <- (-1)
#
#   #determines the pxixj00
#   cdtn00 <- (idx_cat_var1==0)&(idx_cat_var2==0)
#   cdtn11 <- (idx_cat_var1==1)&(idx_cat_var2==1)
#   B_mat[cdtn00, ] <- (-1)* (B_mat[cdtn10, ] +
#                             B_mat[cdtn01, ] +
#                             B_mat[cdtn11, ] )
#   B_mat
# }
###############################################################################


#The function below applies only to Binary data and single-group analysis
Wald_Pearson_test_function <- function(lavobject = fit){
  lavdata        <- lavobject@Data
  lavmodel       <- lavobject@Model
  lavsamplestats <- lavobject@SampleStats
  lavoptions     <- lavobject@Options
  lavcache       <- lavobject@Cache
  lavpartable    <- lavobject@ParTable

  TH      <- lavobject@Fit@TH[[1]]
  th.idx  <- lavobject@Model@th.idx[[1]]
  nvar <- if(lavobject@Model@nexo==0){
    lavobject@Model@nvar
  } else {
    cat("Case with covariates has not been done yet.")
  }
  #model-implied Sigma matrix and Thresholds
  Sigmahat <- lavobject@Fit@Sigma.hat[[1]] #g for multigroup analysis
  Meanhat <- lavobject@implied$mean[[1]]  #g for multigroup analysis
  nsize <- lavobject@Data@nobs[[1]]


  #### SECTION 1 : Compute the positive univariate and bivariate proportions, model-implied probabilities and residuals as well as the positive trivariate and tetra-variate model-implied probabilities ###
  Prob_Prop_Resid <- positiveUniBivTrivTetravModelProb_ObsUnivBivPropResid (
                     lavobject = lavobject,
                       lavdata = lavdata,
                            TH = TH,
                        th.idx = th.idx,
                          nvar = nvar,
                      Sigmahat = Sigmahat,
                       Meanhat = Meanhat)

  prop_1             <- Prob_Prop_Resid$prop_1
  prop_11            <- Prob_Prop_Resid$prop_11
  modelProb_1        <- Prob_Prop_Resid$modelProb_1
  modelProb_11       <- Prob_Prop_Resid$modelProb_11
  modelProb_111      <- Prob_Prop_Resid$modelProb_111
  modelProb_1111     <- Prob_Prop_Resid$modelProb_1111
  residuals_1        <- Prob_Prop_Resid$residuals_1
  residuals_11       <- Prob_Prop_Resid$residuals_11
  UnstdResUniv1Biv11 <- Prob_Prop_Resid$UnstdResUniv1Biv11
  ALLbivModelProb    <- Prob_Prop_Resid$ALLbivModelProb

  #### SECTION 2 : Compute Sigma 2 - the covariance matrix of Univ1 and Biv11 probabilities###
  ################################################################
  Sigma2 <- Cov_pXi1_pXiXj11(modelProb_1 = modelProb_1,
                            modelProb_11 = modelProb_11,
                           modelProb_111 = modelProb_111,
                          modelProb_1111 = modelProb_1111)

  #### SECTION 3 : Compute Delta's - derivatives of Univ1 and Biv11 probabilities wrt theta###
  ################################################################
 derUnivBivProbToTheta <-
 derModelUnivBivProbToTheta (nvar = nvar,
                              TH = TH,
                          th.idx = th.idx,
                       Sigma.hat = Sigmahat,
                        lavcache = lavcache[[1]], #g for multigroup
                        lavmodel = lavmodel)

 Delta2 <- rbind( derUnivBivProbToTheta$derUnivProb1wrtTheta,
                  derUnivBivProbToTheta$derBivProb11wrtTheta)
 derALLBivProbwrtTheta <- derUnivBivProbToTheta$derBivProbwrtTheta

  #### SECTION 4 : Compute Beta indicator matrix###
  ################################################################
  B_mat <- Beta_mat_design(nvar = nvar)


  #### SECTION 5 : Compute the Hessian matrix ###
  ################################################################
  VCOV <-lav_model_vcov(lavmodel = lavmodel,
                           lavsamplestats = lavsamplestats,
                           lavoptions     = lavoptions,
                           lavdata        = lavdata,
                           lavpartable    = lavpartable,
                           lavcache       = lavcache)
  InvHattheta0 <- attr(VCOV, "E.inv")


  #### SECTION 6 : Compute Beta(theta) matrix###
  ################################################################
  B_theta_mat <- t(derALLBivProbwrtTheta * (1/ALLbivModelProb)) %*% B_mat


  #### SECTION 7 : Compute Omega matrix###
  ################################################################
  left_mat <- diag(nrow(Delta2)) - Delta2 %*% InvHattheta0 %*% B_theta_mat
  Omega_mat <- left_mat %*% Sigma2 %*% t(left_mat)


  #### SECTION 8 : Compute the Wald test ###
  ################################################################
  Omega_mat_geninv <- MASS::ginv(Omega_mat)
  Wald_test <- nsize * ( matrix(UnstdResUniv1Biv11, nrow=1) %*%
                         Omega_mat_geninv %*%
                         matrix(UnstdResUniv1Biv11, ncol=1)   )
  AsyMean_Wald <- sum( diag(Omega_mat) )
  AsyVar_Wald <- 2*sum( diag(Omega_mat %*% Omega_mat) )
  FSMadj_Wald <- (AsyMean_Wald / (AsyVar_Wald/2) )* Wald_test
  FSMadj_df_Wald <- (AsyMean_Wald*AsyMean_Wald) / (AsyVar_Wald/2)

  # Version-2 Wald test
  Wald_test_v2 <- nsize * (
    matrix(UnstdResUniv1Biv11, nrow=1) %*%
      diag(1 / diag(Omega_mat)) %*%
      matrix(UnstdResUniv1Biv11, ncol=1)
  )

  #### SECTION 9 : Compute the Pearson test ###
  ################################################################
  D2_inv_mat <- diag( 1/c(modelProb_1, modelProb_11) )
  Pearson_test <- nsize * (matrix(UnstdResUniv1Biv11, nrow=1)   %*%
                           D2_inv_mat %*%
                           matrix(UnstdResUniv1Biv11, ncol=1)       )

  D2_inv_times_Omega <- D2_inv_mat %*% Omega_mat
  AsymMean_Pearson <-  sum( diag(D2_inv_times_Omega) )
  AsymVar_Pearson <-  2*sum( diag(D2_inv_times_Omega %*% D2_inv_times_Omega) )
  FSMadj_Pearson <- (AsymMean_Pearson / (AsymVar_Pearson/2) )* Pearson_test
  FSMadj_df_Pearson <- (AsymMean_Pearson*AsymMean_Pearson) / (AsymVar_Pearson/2)


  list(Wald_test      = c(Wald_test),
       Wald_test_v2   = c(Wald_test_v2),
       FSMadj_Wald    = FSMadj_Wald,
       FSMadj_df_Wald = FSMadj_df_Wald,
       AsyMean_Wald   = AsyMean_Wald,
       AsyVar_Wald    = AsyVar_Wald,

       Pearson_test      = c(Pearson_test),
       FSMadj_Pearson    = FSMadj_Pearson,
       FSMadj_df_Pearson = FSMadj_df_Pearson,
       AsymMean_Pearson  = AsymMean_Pearson,
       AsymVar_Pearson   = AsymVar_Pearson,

       Omega        = Omega_mat,
       Omega_geninv = Omega_mat_geninv,
       Sigma2       = Sigma2,
       Delta2       = Delta2,
       B_theta      = B_theta_mat,
       InvH         = InvHattheta0,
       derALLBivProbwrtTheta = derALLBivProbwrtTheta,

       prop_1  = prop_1,
       prop_11 = prop_11,
       modelProb_1    = modelProb_1,
       modelProb_11   = modelProb_11,
       modelProb_111  = modelProb_111,
       modelProb_1111 = modelProb_1111,
       residuals_1  = residuals_1,
       residuals_11 = residuals_11,
       UnstdResUniv1Biv11 = UnstdResUniv1Biv11,
       ALLbivModelProb = ALLbivModelProb )
}
