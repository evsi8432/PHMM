dives <- unique(df$ID)

pos_dives <- df$ID[df$knownState %in% 6]
neg_dives <- df$ID[df$knownState %in% 5]

rand_pos <- sample(pos_dives)
rand_neg <- sample(neg_dives)

N_pos <- length(pos_dives) 
N_neg <- length(neg_dives)

test_sets <- list()
train_sets <- list()

for(k in 1:K){
  test_sets[[k]] <- c(rand_pos[floor(N_pos*(k-1)/K + 1):floor(N_pos*k/K)],
                      rand_neg[floor(N_neg*(k-1)/K + 1):floor(N_neg*k/K)])
  train_sets[[k]] <- setdiff(dives, test_sets[[k]])
}

if(K == 1){
  train_sets[[1]] <- dives
  test_sets[[1]] <- c(pos_dives,neg_dives)
}

