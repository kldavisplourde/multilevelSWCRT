The ZIP folder contains files for implementing the sample size and power calculations under a generalized linear mixed model (GLMM) framework as shown in "Sample size considerations for stepped wedge designs with subclusters" by Davis-Plourde, Taljaard, and Li.

List of Files:
1) VARd.r = R file for generating the variance of the intervention effect for Gaussian or binary (with logit link) outcomes under a SW-CRT with subclusters and equal cluster sizes.
2) VARd_CV.r = R file for generating the variance of the intervention effect for Gaussian or binary (with logit link) outcomes taking into account coefficients of variation (CVs) in terms of the number of subclusters and subcluster sizes under a SW-CRT with subclusters.
3) study_power.r = R file for generating power predictions under the noncentral t-distribution.
4) designEffect.r = R file for generating the design effect comparing a longitudinal CRT with subclusters to individual randomization (no assumption made on randomization schedule).
5) designEffect_SWCRT.r = R file for generating the design effect comparing a SW-CRT with subclusters to individual randomization.
6) Simulation_gaussian.r = R file for simulating and fitting SW-CRTs with subclusters under a Gaussian outcome with equal cluster sizes.
7) Simulation_binomial.r = R file for simulating and fitting SW-CRTs with subclusters under a binary outcome with logit link and equal cluster sizes.
8) Simulation_unequalCS_gaussian.r = R file for simulating and fitting SW-CRTs with subclusters under a Gaussian outcome and taking into account CVs in terms of the number of subclusters and subcluster sizes.
9) Simulation_unequalCS_binomial.r = R file for simulating and fitting SW-CRTs with subclusters under a binary outcome with logit link and taking into account CVs in terms of the number of subclusters and subcluster sizes.
10) Application_Study.r = R file for producing power and sample size estimates for two SW-CRTs with subclusters, the LIRE trial and the Washington State EPT trial.

NOTES:  1) This program requires the MASS package (comes preloaded into R so does not require installation). Simulation files additionally require the installation of the lme4, foreach, doMC, and doRNG packages.
	2) You will need to change path names before running the programs.
	3) Latest version of all files are available on GitHub: https://github.com/kldavisplourde/multilevelSWCRT
	4) R Shiny App for the SW-CRT design effect: https://kendra-davis-plourde.shinyapps.io/SWCRT_3Level_DesignEffect/

This work is supported by the National Institute of Aging (NIA) of the National Institutes of Health (NIH) under Award Number U54AG063546, which funds NIA Imbedded Pragmatic Alzheimer's Disease and AD-Related Dementias Clinical Trials Collaboratory (NIA IMPACT Collaboratory). The content is solely the responsibility of the authors and does not necessarily represent the official views of the NIH. The authors are grateful to Professor Jim Hughes for providing data and information from the Washington State EPT trial. We also thank the associate editor and an anonymous referee for their valuable suggestions, which greatly improved the exposition of this work.