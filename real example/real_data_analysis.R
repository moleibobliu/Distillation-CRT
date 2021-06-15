
##### Read file and pre-process ######

# The gene list file 'ncomms11479-s4.csv' was downloaded from https://www.nature.com/articles/ncomms11479#Sec32
# The data files were from https://www.cbioportal.org/study/summary?id=brca_metabric


######## Read the data ########

# list of genes 
gene_lst <- read.csv('ncomms11479-s4.csv', head = T)
gene_lst <- as.character(gene_lst$X.Supplementary.Dataset.3...Matrix.for.mutations.across.all.genes.and.samples...NA.no.coding.mutation...For.inframe.indels.and.missense.SNVs..the.distinction.between.recurrent.and.non.recurrent.events.is.made.)

# CNA and RNA

cna_data <- read.csv('data_CNA.csv', head = T)
RNA_data <- fread('data_expression.txt')
gene_lst <- intersect(intersect(RNA_data$Hugo_Symbol, gene_lst),  cna_data[,1])

# We use the intersection of the gene_lst and the measured genes in cna_data and RNA data as the covariates

gene_indx_cna <- which(as.character(cna_data[,1]) %in% gene_lst)
gene_indx_rna <- which(as.character(RNA_data$Hugo_Symbol) %in% gene_lst)

cna_data <- cna_data[gene_indx_cna, ]
cna_data <- t(cna_data)
RNA_data <- RNA_data[gene_indx_rna,]
RNA_data <- t(RNA_data)
RNA_data <- RNA_data[,match(cna_data[1,], RNA_data[1,])]

for (j in 2:length(rownames(cna_data))){
  string <- strsplit(rownames(cna_data)[j], split = '[.]')
  rownames(cna_data)[j] <- paste(string[[1]][1], '-', string[[1]][2], sep = '')
}

# Clinical outcome

clinical_data <- read.table('data_clinical_sample.csv', head = T,
                            sep = ',')

#### Merge the data by patient ID #####

pat_set <- intersect(intersect(as.character(clinical_data$PATIENT_ID),
                               rownames(RNA_data)), rownames(cna_data))

X_rna <- c()
X_cna <- c()
ER_status <- c()
Y_clinical <- c()

for (i in 1:length(pat_set)) {
  pat_id <- pat_set[i]
  indx_rna <- which(rownames(RNA_data) == pat_id)
  indx_cna <- which(rownames(cna_data) == pat_id)
  indx_clinic <- which(as.character(clinical_data$PATIENT_ID) == pat_id)
  ER_status <- c(ER_status, clinical_data$ER_STATUS[indx_clinic])
  
  ## choosing the outcome
  Y_clinical <- rbind(Y_clinical, 
                      c(clinical_data$TUMOR_SIZE[indx_clinic], clinical_data$GRADE[indx_clinic],
                        clinical_data$TUMOR_STAGE[indx_clinic]))
  
  X_cna <- rbind(X_cna, as.vector(as.numeric(cna_data[indx_cna,])))
  X_rna <- rbind(X_rna, as.vector(as.numeric(RNA_data[indx_rna,])))
}


X_rna_pos <- X_rna
X_cna_pos <- X_cna


##### Adjust the RNA data using the CNA data #####

piecewise_linear <- function(X_rna, X_cna){
  p <- length(X_rna[1,])
  R_sq_lst <- c()
  for (j in 1:p) {
    R_tot <- var(X_rna[,j])
    Y1 <- mean(X_rna[which(X_cna[,j] <= -1), j])
    Y2 <- mean(X_rna[which(X_cna[,j] == 0), j])
    Y3 <- mean(X_rna[which(X_cna[,j] >= 1), j])
    X_rna[which(X_cna[,j] <= -1), j] <- X_rna[which(X_cna[,j] <= -1), j] - Y1
    X_rna[which(X_cna[,j] == 0), j] <- X_rna[which(X_cna[,j] == 0), j] - Y2
    X_rna[which(X_cna[,j] >= 1), j] <- X_rna[which(X_cna[,j] >= 1), j] - Y3
    R_res <- var(X_rna[,j])
    R_sq_lst <- c(R_sq_lst, (R_tot - R_res) / R_tot)
  }
  return(list(R = R_sq_lst, X_rna = X_rna))
}

adj_result <- piecewise_linear(X_rna_pos, X_cna_pos)
X_rna_new <- adj_result$X_rna

# The processed data:
data_use <- as.data.frame(cbind(X_rna_new, ER_status, Y_clinical))

N <- length(data_use[,1])
p <- length(data_use[1,])
data_use <- data_use[,c(1:165,167)]
data_use <- data_use[which(complete.cases(data_use)),]




########### Extract the final data set ###########

ER_lst <- data_use$ER_status
X <- as.matrix(data_use[,1:length(gene_use)])
Y <- data_use[,length(data_use[1,])]

# We only use patients with postive ER status

ER_pos_set <- which(ER_lst == 2)
X_pos <- X[ER_pos_set,]
X_pos <- normalize(X_pos)
p <- length(X_pos[1,])

# Merge two categories in Y since one of them is of very small size:

Y_pos <- Y[ER_pos_set]
Y_binary <- ifelse(Y_pos == 3, 1, 0)


# We finally use X_pos and Y_binary

#save(X_pos, Y_binary, file = 'real_example.rda')


#### Model X as joint gaussian using graphic lasso ####

library(CVglasso)
prec_est <- CVglasso(X = X_pos, lam.min.ratio = 1e-6)
length(which(prec_est$Omega != 0))
Sigma_est <- prec_est$Sigma

var_lst <- c()
for (indx in 1:p){
  beta_x <- solve(Sigma_est[-indx, -indx], Sigma_est[-indx, indx])
  X_bar <- X_pos[ ,-indx] %*% beta_x
  var_lst <- c(var_lst, mean((X_pos[,indx] - X_bar)^2))
  print(indx)
}
theta_lst <- diag(solve(Sigma_est))
d <- sqrt(theta_lst * var_lst)
Sigma_est <- as.matrix(diag(d)) %*% as.matrix(Sigma_est) %*% as.matrix(diag(d))



##### Run analysis with our proposed method and the benchmark methods #####
#### We ran the following codes for 300 replications in our paper. Strongly recommand using server for parallel.

for (seed.id in 1:300){
  
  time_lst <- c()
  rej_FDR_mat <- c()
  rej_FWER_mat <- c()
  pvl_mat <- c()
  
  set.seed(seed.id)
  
  ################## d0CRT ################## 
  
  
  ptm <- proc.time()[1]
  rCRT_results <- dCRT(X_pos, Y_binary, Sigma_X = Sigma_est, FDR = 0.1, model = 'Binomial_lasso')
  time_lst <- c(time_lst, proc.time()[1] - ptm)
  
  rej_FDR <- rep(0, ncol(X_pos))
  rej_FDR[rCRT_results$select_set] <- 1
  
  rej_FWER <- rep(0, ncol(X_pos))
  rej_FWER[which(rCRT_results$p_values < 0.1 / ncol(X_pos))] <- 1
  
  rej_FWER_mat <- cbind(rej_FWER_mat, rej_FWER)
  rej_FDR_mat <- cbind(rej_FDR_mat, rej_FDR)
  pvl_mat <- cbind(pvl_mat, rCRT_results$p_values)
  
  
  
  ################## dkCRT ################## 
  
  ptm <- proc.time()[1]
  dICRT_results <- dCRT(X_pos, Y_binary, Sigma_X = Sigma_est, FDR = 0.1, d.interaction = T,
                        model = 'Binomial_lasso')
  time_lst <- c(time_lst, proc.time()[1] - ptm)
  
  rej_FDR <- rep(0, ncol(X_pos))
  rej_FDR[dICRT_results$select_set] <- 1
  
  rej_FWER <- rep(0, ncol(X_pos))
  rej_FWER[which(dICRT_results$p_values < 0.1 / ncol(X_pos))] <- 1
  
  rej_FWER_mat <- cbind(rej_FWER_mat, rej_FWER)
  rej_FDR_mat <- cbind(rej_FDR_mat, rej_FDR)
  pvl_mat <- cbind(pvl_mat, dICRT_results$p_values)
  
  
  ################### HRT ####################
  
  ptm <- proc.time()[1]
  HRT_results <- HRT(Y_binary, X_pos, Sigma = Sigma_est, FDR = 0.1, N = 25000, 
                     model = 'binomial')
  time_lst <- c(time_lst, proc.time()[1] - ptm)
  
  
  rej_FDR <- rep(0, ncol(X_pos))
  rej_FDR[HRT_results$select_set] <- 1
  
  rej_FWER <- rep(0, ncol(X_pos))
  rej_FWER[HRT_results$select_FWER] <- 1
  
  rej_FWER_mat <- cbind(rej_FWER_mat, rej_FWER)
  rej_FDR_mat <- cbind(rej_FDR_mat, rej_FDR)
  pvl_mat <- cbind(pvl_mat, HRT_results$p_values)
  
  
  ###############  Original CRT ###################
  ptm <- proc.time()[1]
  o_CRT <- CRT_sMC(Y_binary, X_pos, Sigma_est, m = 25000, type = 'LASSO', 
                   FDR = 0.1, model = 'binomial')
  time_lst <- c(time_lst, proc.time()[1] - ptm)
  
  
  rej_FDR <- rep(0, ncol(X_pos))
  rej_FDR[o_CRT$select_set] <- 1
  
  rej_FWER <- rep(0, ncol(X_pos))
  rej_FWER[which(o_CRT$pvl < 0.1 / ncol(X_pos))] <- 1
  
  rej_FWER_mat <- cbind(rej_FWER_mat, rej_FWER)
  rej_FDR_mat <- cbind(rej_FDR_mat, rej_FDR)
  pvl_mat <- cbind(pvl_mat, o_CRT$pvl)
  
  
  ################## Knockoff #################
  
  ptm <- proc.time()[1]
  knockoffs = function(X) create.gaussian(X_pos, rep(0, p), Sigma_est, method = 'sdp',
                                          diag_s =  create.solve_sdp(Sigma_est))
  foo = stat.glmnet_coefdiff
  k_stat = function(X, X_k, y) foo(X, X_k, y, family = 'binomial')
  Knockoff_Result <- knockoff.filter(X_pos, Y_binary, statistic = k_stat, 
                                     knockoffs = knockoffs, fdr = 0.1, offset = 1)
  time_lst <- c(time_lst, proc.time()[1] - ptm)
  
  rej_FDR <- rep(0, ncol(X_pos))
  rej_FDR[Knockoff_Result$selected] <- 1
  rej_FDR_mat <- cbind(rej_FDR_mat, rej_FDR)
  
  
  filename <- paste('real_exp_seed', seed.id, '.rda', sep = '')
  
  save(rej_FDR_mat, rej_FWER_mat, pvl_mat, time_lst, file = filename)
  
}



#### Summarize the results and making plots ####

## Read the results from the .rda files:

num_selected_FDR <- c()
num_selected_FWER <- c()
pvl_all <- 0
time_all <- 0
names <- c('dCRT', 'HRT', 'CRT', 'Knockoff')

select_FDR_mat <- 0
for (j in 1:300){
  filename <- paste('real_exp_seed', j, '.rda', sep = '')
  load(filename)
  num_selected_FDR <- rbind(num_selected_FDR, colSums(rej_FDR_mat))
  num_selected_FWER <- rbind(num_selected_FWER, colSums(rej_FWER_mat))
  select_FDR_mat <- select_FDR_mat + rej_FDR_mat
  pvl_all <- pvl_all + pvl_mat
  time_all <- time_all + time_lst
}
pvl_all <- pvl_all / 300
time_all <- time_all / 300


## Make mass point plots for the number of detection of the methods

############# FDR ###############

set_num <- sort(unique(as.vector(num_selected_FDR)))
num_mat <- matrix(0, ncol(num_selected_FDR), length(set_num))

for (i in 1:length(set_num)) {
  num <- set_num[i]
  for (j in 1:ncol(num_selected_FDR)) {
    num_mat[j,i] <- length(which(num_selected_FDR[,j] == num))
  }
}


dat1 <- as.data.frame(cbind(set_num, rep('d[0]CRT', length(set_num))))
colnames(dat1) <- c('x', 'y')  
dat1$x <- factor(as.character(set_num), levels = as.character(set_num))

dat2 <- as.data.frame(cbind(set_num, rep('d[I]CRT', length(set_num))))
colnames(dat2) <- c('x', 'y')
dat2$x <- factor(as.character(set_num), levels = as.character(set_num))

dat3 <- as.data.frame(cbind(set_num, rep('HRT', length(set_num))))
colnames(dat3) <- c('x', 'y')  
dat3$x <- factor(as.character(set_num), levels = as.character(set_num))

dat4 <- as.data.frame(cbind(set_num, rep('oCRT', length(set_num))))
colnames(dat4) <- c('x', 'y')  
dat4$x <- factor(as.character(set_num), levels = as.character(set_num))

dat5 <- as.data.frame(cbind(set_num, rep('knockoffs', length(set_num))))
colnames(dat5) <- c('x', 'y')  
dat5$x <- factor(as.character(set_num), levels = as.character(set_num))



vec1 <- sqrt((num_mat[1,] + ifelse(num_mat[1,] > 0, 0, 0)) / 1.5)
vec2 <- sqrt((num_mat[2,] + ifelse(num_mat[2,] > 0, 0, 0)) / 1.5)
vec3 <- sqrt((num_mat[3,] + ifelse(num_mat[3,] > 0, 0, 0)) / 1.5)
vec4 <- sqrt((num_mat[4,] + ifelse(num_mat[4,] > 0, 0, 0)) / 1.5)
vec5 <- sqrt((num_mat[5,] + ifelse(num_mat[5,] > 0, 0, 0)) / 1.5)


p <- ggplot(dat1, aes(x=x, y=y)) + 
  geom_point(data = dat1, size = vec1, color = ifelse(num_mat[1,] > 0, 'black', 'white')) +
  geom_point(data = dat2, size = vec2, color = ifelse(num_mat[2,] > 0, 'black', 'white')) +
  geom_point(data = dat3, size = vec3, color = ifelse(num_mat[3,] > 0, 'black', 'white')) + 
  geom_point(data = dat4, size = vec4, color = ifelse(num_mat[4,] > 0, 'black', 'white')) +
  geom_point(data = dat5, size = vec5, color = ifelse(num_mat[5,] > 0, 'black', 'white')) + 
  labs(title = "Multiple Testing with FDR Control", y = "Methods",
       x = 'Number of detections') +
  theme(axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 25),
        axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 25),
        plot.title = element_text(size = 35, hjust = 0.5))
p
ggsave(filename = '../plots/FDR_real_example.pdf', 
       plot = p, width = 10, height =  7.5)





############# FWER ###############


set_num <- sort(unique(as.vector(num_selected_FWER)))
num_mat <- matrix(0, ncol(num_selected_FWER), length(set_num))

for (i in 1:length(set_num)) {
  num <- set_num[i]
  for (j in 1:ncol(num_selected_FWER)) {
    num_mat[j,i] <- length(which(num_selected_FWER[,j] == num))
  }
}


dat1 <- as.data.frame(cbind(set_num, rep('d0CRT', length(set_num))))
colnames(dat1) <- c('x', 'y')  
dat1$x <- factor(as.character(set_num), levels = as.character(set_num))

dat2 <- as.data.frame(cbind(set_num, rep('d[I]CRT', length(set_num))))
colnames(dat2) <- c('x', 'y')
dat2$x <- factor(as.character(set_num), levels = as.character(set_num))

dat3 <- as.data.frame(cbind(set_num, rep('HRT', length(set_num))))
colnames(dat3) <- c('x', 'y')  
dat3$x <- factor(as.character(set_num), levels = as.character(set_num))

dat4 <- as.data.frame(cbind(set_num, rep('oCRT', length(set_num))))
colnames(dat4) <- c('x', 'y')  
dat4$x <- factor(as.character(set_num), levels = as.character(set_num))

vec1 <- sqrt((num_mat[1,] + ifelse(num_mat[1,] > 0, 0, 0)) / 1.5)
vec2 <- sqrt((num_mat[2,] + ifelse(num_mat[2,] > 0, 0, 0)) / 1.5)
vec3 <- sqrt((num_mat[3,] + ifelse(num_mat[3,] > 0, 0, 0)) / 1.5)
vec4 <- sqrt((num_mat[4,] + ifelse(num_mat[4,] > 0, 0, 0)) / 1.5)

p <- ggplot(dat1, aes(x=x, y=y)) + 
  geom_point(data = dat1, size = vec1, color = ifelse(num_mat[1,] > 0, 'black', 'white')) +
  geom_point(data = dat2, size = vec2, color = ifelse(num_mat[2,] > 0, 'black', 'white')) +
  geom_point(data = dat3, size = vec3, color = ifelse(num_mat[3,] > 0, 'black', 'white')) + 
  geom_point(data = dat4, size = vec4, color = ifelse(num_mat[4,] > 0, 'black', 'white')) +
  labs(title = "Multiple Testing with FWER Control", y = "Methods",
       x = 'Number of detections') +
  theme(axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 25),
        axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 25),
        plot.title = element_text(size = 35, hjust = 0.5))
p
ggsave(filename = '../plots/FWER_real_example.pdf', 
       plot = p, width = 10, height =  7.5)






