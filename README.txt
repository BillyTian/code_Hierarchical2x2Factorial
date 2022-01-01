This repository includes R code for reproducing Table 3 to 5, and Figure 1 in the article "Sample size calculation in hierarchical 2x2 factorial trials with unequal cluster sizes" (Stat Med accepted), and Web Table 1-27 in its supplementary material.

For questions or comments about the code, please contact Zizhong Tian at <zizhong.tian@yale.edu>.

I. Supporting Files: These supporting files are sourced in the corresponding main files that reproduce the simulation results in the main manuscript as well as the supplementary web appendix.

1) gendata.R = functions to generate the simulated data in each iteration based on three types of CV distributions;
2) CE_X.R = functions to calculate the empirical power and empirical type I error rate about the test for the controlled effect of the cluster-level treatment with and without finite-sample correction;
3) CE_Z.R = function to calculate the empirical power and empirical type I error rate about the test for the controlled effect of the individual-level treatment;
4) NE_X.R = functions to calculate the empirical power and empirical type I error rate about the test for the natural effect of the cluster-level treatment with and without finite-sample correction;
5) NE_Z.R = function to calculate the empirical power and empirical type I error rate about the test for the natural effect of the individual-level treatment;
6) interaction.R = function to calculate the empirical power and empirical type I error rate about the interaction test;
7) CE_joint.R = functions to calculate the empirical power and empirical type I error rate about the joint test of the two controlled effects with and without finite-sample correction;
8) NE_joint.R = functions to calculate the empirical power and empirical type I error rate about the joint test of the two natural effects with and without finite-sample correction;
9) CE_IU.R = functions to calculate the empirical power and empirical type I error rate about the intersection-union (I-U) test of the two controlled effects with and without finite-sample correction;
10) NE_IU.R = functions to calculate the empirical power and empirical type I error rate about the I-U test of the two natural effects with and without finite-sample correction.

II. Main Files: These main files are used to reproduce the simulation results and the illustrative plot in the main manuscript as well as the web appendix.

11) Simulation_T3_WT1_WT14_WT15_WT21_WT22 (CE_X).R = reproduce simulation results in Table 3, Web Table 1, 14, 15, 21, and 22;
12) Simulation_T4_WT3_WT17_WT18_WT24_WT25 (NE_X).R = reproduce simulation results in Table 4, Web Table 3, 17, 18, 24, and 25;
13) Simulation_T5_WT20_WT27 (IE_XZ).R = reproduce simulation results in Table 5, Web Table 20 and 27;
14) Simulation_WT2_WT16_WT23 (CE_Z).R = reproduce simulation results in Web Table 2, 16, and 23;
15) Simulation_WT4_WT19_WT26 (NE_Z).R = reproduce simulation results in Web Table 4, 19, and 26;
16) Simulation_WT6_WT7 (CE_joint).R = reproduce simulation results in Web Table 6 and 7;
17) Simulation_WT8_WT9 (CE_IU).R = reproduce simulation results in Web Table 8 and 9;
18) Simulation_WT10_WT11 (NE_joint).R = reproduce simulation results in Web Table 10 and 11;
19) Simulation_WT12_WT13 (NE_IU).R = reproduce simulation results in Web Table 12 and 13;
20) Application_Figure1.R = reproduce the illustrative plot in the application part, Figure 1;

III. Software 

Analyses were conducted with R, version 4.0.3 (https://www.r-project.org/)
The calculations used R packages nlme (version 3.1-143) and Matrix (version 1.2-18).

IV. R commands for the installation of R packages 

install.packages(c("nlme", "Matrix", "mvtnorm")) 
