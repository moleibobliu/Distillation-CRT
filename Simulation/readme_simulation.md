Scripts simulations.R contains the setups and instructions for generating data and implementation of the methods in our simulation studies. We recommand using multiple clusters for parallel and replicate each setting for many times (says 300), and summarizing the paralleled results to estimate the average power, FDR and consumed time of the included algorithms. 

The setups inlcude: 

1. Moderate size simulations;
2. Large size simulations;
3. Power of dICRT in the presence of interaction (linear model and random forest);
4. Robustness evaluation (In-sample estimate of X's model, and known first and second moments);
5. Sensitiviy of dICRT to the choice of k (dimensionality of the distilled data);
6. dCRT implementation with no screening for dimension reduction;
7. Compare oCRT with dCRT and variation of oCRT (studies in Section 4.1 of our paper);
8. Assess algorithmic variability;
9. dCRT: resample-free v.s. non-resample-free (gamma and binary X);
10. Compare with DML/GCM with small size and laplacian tailed data.
