Data_fine_unlabelled <- prepData(Data_fine[Data_fine$divenum %in% test_dives,
                                      c("ID","divenum","stime","ad",names(dist))],
                            coordNames=NULL)

Data_fine_unlabelled$true_label <- Data_fine_unlabelled$knownState
Data_fine_unlabelled$knownState <- 7

Par0 <- getPar0(hmm)

hmm0 <- fitHMM(data=Data_fine_unlabelled,
               nbStates=N,
               dist=hmm$conditions$dist,
               fixPar=list(knownState = hmm$conditions$fixPar$knownState),
               DM=hmm$conditions$DM,
               beta0=Par0$beta,
               delta0=(1-eps)*(Par0$delta)+eps*rep(1/N,N),
               Par0=Par0$Par,
               nlmPar = list('stepmax'=1e-100,
                             'iterlim'=1))

Data_fine_unlabelled$viterbi <- viterbi(hmm0)
Data_fine_unlabelled$p_catch <- stateProbs(hmm0)[,4]
Data_fine_unlabelled$p_no_fish <- stateProbs(hmm0)[,5]
Data_fine_unlabelled$p_fish <- stateProbs(hmm0)[,6]

probs <- c()
labs <- c()

conf_matrices_PHMM[[k]] = matrix(0, nrow = 2, ncol = 2)
rownames(conf_matrices_PHMM[[k]]) <- c("True Positive", "True Negative")
colnames(conf_matrices_PHMM[[k]]) <- c("Pred Positive", "Pred Negative")

for(divenum in test_dives){
  
  if(divenum %in% pos_dives){
    rownum <- 1
    labs <- c(labs,TRUE)
  } else if (divenum %in% neg_dives){
    rownum <- 2
    labs <- c(labs,FALSE)
  }
  
  p_fish <- tail(Data_fine_unlabelled$p_fish[Data_fine_unlabelled$ID == divenum],1)
  print(p_fish)
  
  probs <- c(probs,p_fish)
  conf_matrices_PHMM[[k]][rownum,1] <- conf_matrices_PHMM[[k]][rownum,1] + p_fish
  conf_matrices_PHMM[[k]][rownum,2] <- conf_matrices_PHMM[[k]][rownum,2] + 1 - p_fish
  
}

AUCs_PHMM[k] <- roc(response = labs, predictor=probs)$auc
