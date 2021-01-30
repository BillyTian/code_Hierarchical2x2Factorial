# code_Hierarchical2x2Factorial
This repository includes R code for reproducing Table 3, 4, and Figure 1 in the article "Sample size calculation in hierarchical 2x2 factorial trials with unequal cluster sizes" (under review), and Web Table 1-6 in its supplementary material.

For questions or comments about the code, please contact Zizhong Tian at <zizhong.tian@yale.edu>.

I. Supporting Files: These supporting files are sourced in the main files that reproduce the simulation results in the manuscript.

1) marginal_cluster.R = function to calculate the empirical power and empirical type I error rate about the test for marginal cluster-level treatment effect without finite-sample correction;
2) marginal_cluster_corrected.R = function to calculate the empirical power and empirical type I error rate about the test for marginal cluster-level treatment effect with finite-sample correction;
3) marginal_ind.R = function to calculate the empirical power and empirical type I error rate about the test for marginal individual-level treatment effect;
4) interaction.R = function to calculate the empirical power and empirical type I error rate about the interaction test of the two treatments;
5) joint.R = function to calculate the empirical power and empirical type I error rate about the joint test without finite-sample correction;
6) joint_corrected.R = function to calculate the empirical power and empirical type I error rate about the joint test with finite-sample correction;
7) IU.R = function to calculate the empirical power and empirical type I error rate about the intersection-union (I-U) test without finite-sample correction;
8) IU_corrected.R = function to calculate the empirical power and empirical type I error rate about the I-U test with finite-sample correction.

II. Main Files: These main files are used to reproduce the simulation results and the illustrative plot in the manuscript.

9) Simulation_WebTable1_and_2.R = reproduce simulation results in Web Table 1 and 2 (corresponds to the marginal cluster-level test);
10) Simulation_WebTable3.R = reproduce simulation results in Web Table 3 (corresponds to the marginal individual-level test);
11) Simulation_Table3.R = reproduce simulation results in Table 3 (corresponds to the interaction test);
12) Simulation_WebTable4_and_5.R = reproduce simulation results in Web Table 4 and 5 (corresponds to the joint test);
13) Simulation_WebTable6_and_Table4.R = reproduce simulation results in Web Table 6 and Table 4 (corresponds to the I-U test);
14) Application_Figure1.R = reproduce the illustrative plot in the application part, Figure 1;

III. Software 

Analyses were conducted with R, version 3.6.2 (https://www.r-project.org/)
The calculations used R packages nlme (version 3.1-143) and Matrix (version 1.2-18).

IV. R commands for the installation of R packages 

install.packages(c("nlme", "Matrix")) 
