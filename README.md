# CABS
Conditional Adaptive Bayesian Spectrum Analysis (CABS)  
"Adaptive Bayesian Spectral Analysis of Nonstationary Biomedical Time Series"  
by Bruce, Hall, Buysse, and Krafty (2017)  
  
Author: Scott A. Bruce  
  
## Description:
Instructions for implementing CABS estimation from "Adaptive Bayesian  
Spectral Analysis of Nonstationary Biomedical Time Series"  
by Bruce, Hall, Buysse, and Krafty (2017)  
  
## Dependencies:  
Code was developed using MATLAB 2016b (version 1.0.0) and has since been updated
to be compatible with MATLAB 2018a (version 2.0.0), so code may not function properly  
on other versions of MATLAB.  If choose parallel computing option (parcomp=1), the Parallel   
Computing Toolbox is required.  
  
Code has been tested and developed for data of the following dimensions:  
<= 5000 time points, <= 100 subjects, <= 50000 iterations.  
Larger datasets may require more RAM for processing.  Reduce batch sizes  
in MCMC settings to reduce required RAM by outputting more frequently  
to dataset stored in ROM.    
      
## Quick start guide:  
Follow the steps below to simulate data, run the CABS estimation procedure,  
and create data visualizations to assess the convergence and fit of the   
method to simulated slowly varying and piecewise AR processes detailed in   
the paper.  
  
1) Download the .zip file and extract the contents to a folder of your  
choosing.  In what follows, we will assume the folder is saved with the   
following path 'C:\CABSdemo'.  
  
2) Open MATLAB (version 2016b or newer is recommended) and change your   
current working directory to 'C:\CABSdemo\programs'.  This   
directory contains all of the MATLAB functions needed to run the MCMC sampler.  
  
3) Open the file 'C:\CABSdemo\demo.m'.  This script contains   
all necessary code to reproduce the data and MCMC sampler results for the   
piecewise and slowly varying AR processes.   

4) Follow the instructions in the comments of the demo file to simulate  
data from the piecewise and slowly varying AR processes described in the  
paper and apply the CABS procedure to obtain spectral estimates and summary  
plots of the MCMC sampler fit. There is an additional section which   
demonstrates how to obtain the convergence diagnostic measures and plots   
used to assess convergence.  
  
## Using CABS on other data:  
Three inputs are necessary to run the MCMC sampler for the CABS estimation  
procedure.  
  
1) '**x**' is a matrix whereby each column contains a realization of the time  
series for a given subject.  The rows are indexed by time   
(t=1,...,T) and the columns are indexed by subject (l=1,...,L).  For   
example, {x}_11 contains the first time series instance for the first   
subject.    
  
2) '**u**' is a vector containing the corresponding covariate values for each  
subject.  The columns are indexed by subject, which is the same as for  
the time series matrix 'x'.  For example, {u}_1 is the covariate value  
for the first subject.  
  
3) '**opt**' is a structure containing the options required for the MCMC   
sampler to run.  See 'setMCMCOptions.m' for descriptions of all of the  
options as well as the default values.  Once you have determined the  
appropriate settings for your data, use this function to create the  
structure needed as in the demo file.  
  
Once you have created these three inputs based on your data, you can pass  
them into the 'CABS' function just as in the demo file for estimation.  
