Scripts real_data_analysis.R include several parts:

1. Load and preprocess the CNA, mRNA and phenotype (caner) data;

2. Adjust mRNA with the CNA data and estimate the covariance structure of the adjusted expression level using graphi lasso;

3. Implement the algorithms with the processed data and estimated covariance structure (for 300 replications);

4. Summarize the results: obtain the detected variables, their p-values and the point mass plot for the number of detections of the algorithms shown in the paper.
