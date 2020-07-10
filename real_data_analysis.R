
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


######## run d0CRT ########

d0CRT_results <- dCRT(X_pos, Y_binary, Sigma_X = Sigma_est, FDR = 0.1, model = 'Binomial_lasso')

######## run dICRT ########

dICRT_results <- dCRT(X_pos, Y_binary, Sigma_X = Sigma_est, FDR = 0.1, d.interaction = T,
                     model = 'Binomial_lasso')

######## HRT #######

set.seed(1)

HRT_results <- HRT(Y_binary, X_pos, Sigma = Sigma_est, FDR = 0.1, N = 25000, 
                   model = 'binomial')

######## original CRT with lasso #######

original_CRT <- CRT_sMC(Y_binary, X_pos, Sigma_est, m = 25000,  
                        type = 'LASSO', FDR = 0.1, model = 'binomial')

######## knockoffs ######## 

set.seed(666)

knockoffs = function(X) create.gaussian(X_pos, rep(0, p), Sigma_est, method = 'sdp',
                                        diag_s =  create.solve_sdp(Sigma_est))
foo = stat.glmnet_coefdiff
k_stat = function(X, X_k, y) foo(X, X_k, y, family = 'binomial')
Knockoff_Result <- knockoff.filter(X_pos, Y_binary, statistic = k_stat, knockoffs = knockoffs,
                                   fdr = 0.1, offset = 1)
Knockoff_Result$selected







