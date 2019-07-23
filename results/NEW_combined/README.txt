Description of the new result files.

GENERAL COMMENTS
* For all the studies -- except different threshold exercises in 
iid_MSE_dif_thresholds.xlsx -- DM test statistics are reported. 
* For cell (i,j) the loss differentials for DM statistics 
are calculated as value_of_method_i - value_of_method_j [row method - column method].
* For simulation studies RMSE are considered, for empirical studies
CSR of Diks et al. (2011).
* For convenience colours are used to distinguish significance levels.
For "good" results shades of green are used, for "bad" results shared of red. The darker the shade the higher the significance (the darkest for 1%, medium for 5%, the lightest for 10%). Gray coloured values are neither good 
nor bad (differences between methods of the same type),
they simply indicate that some results are significant.


A. SIMULATION STUDIES
* All the results are for 99.5%, 99% and 95% VaR.
* For iid studies 100 MC replications were considered.
For AR(1) exercises with T=100 and T=1000 -- 50 MC replications, for T=10000 -- 20 MC replications.
* For iid studies only the regular posterior and fully censored posterior (CP) we used.
For AR(1) studies also partially censored posterior (PCP) was considered.
* There are the following files:

mix_iid_DM_lowerLambdas.xlsx -- 
-- iid data for the DGP, generated from a mixture distribution
-- mixture as the one used for GARCH errors by Ausin & Galeano 
-- inverse of lambda amplifies volatility of the "bad state"
-- probability of the "good state" rho = 0.75
-- lambdas used from 0.05 by 0.05 to 4.5, then from 0.5 by 0.1 to 0.9
-- three spreadsheets, for T=100, T=1000, T=10000

skt_iid_DM_moreLambda.xlsx -- 
-- iid data for the DGP, generated from a skewed t distribution
-- Hansen's skewed t distribution
-- domain of lambda (the skewness parameter) is (-1,1)
-- lambdas considered from -0.6 by 0.1 to 0.6
-- three spreadsheets, for T=100, T=1000, T=10000

iid_MSE_dif_thresholds.xlsx -- 
-- iid data for the DGP 
-- two spreadsheets:
one for split normal distribution with different sigma2, from 1.5 by 0.5 to 3.5;
one for skewed t distribution, with lambdas from -0.5 by 0.1 to 0.5
-- thresholds for CP: fixed at 0, at the 5th, 10th, 15th, 20th percentile of the in-sample data.
-- only MSEs reported! (DM tables would be not too convenient to read, large and many of them)

skt_ar1_DM.xlsx -- 
-- AR(1) data for the GDP, specification as in the paper
except that Hansen's skewed t distribution used for the error terms
-- lambdas used from -0.5 by 0.1 to 0.5
-- three spreadsheets, for T=100, T=1000, T=10000
-- estimation methods as in the paper: regular posterior, fully CP, PCP;
thresholds for the latter two at 0 and at the 10th data percentile


B. EMPIRICAL STUDIES
* All the results are for CSR with threshold at 99.5%, 99% and 95% 
(which is called "threshold for evaluation"). 
* Thresholds for evaluation are either time constant (based on the given 
percentile of the in-sample data) 
or time varying (based on the given percentile of the MLE-implied predictive distribution).
* Two  types of threshold for estimation were used: 
time-constant (based on the 10th, 20th, 30th, 40th and 50th percentile of the in-sample data)
and time-varying (based on the xth percentile of the MLE-implied predictive distribution, 
with x={20,30,40,50} for AGARCH-skt 
and x={40,50} for GAS-t, for which x=30 MitISEM approximation did not succeed).
* IBM data from the paper plus 507 new observations from 2017--2018,
T = 1000 as before, H=1500 (previous) + 507 (new) = 2007
[hope was that a longer out of sample would give us more significance]. 
* three types of spreadsheets: 
for time-constant evaluation, 
for time-varying evaluation 
and estimation results (means and standard deviations).
* Additional notation:
underlined are the result comparing the CP and PCP with the same $x$ threshold;
gray fields show cased for which time-varying threshold for estimation performs better than the time-constant threshold with the same x
* There are the following files:

DM_skt_agarch_H2007out_of_sample.xlsx -- ARARCH-skt

DM_t_gas_H2007out_of_sample -- GAS-t