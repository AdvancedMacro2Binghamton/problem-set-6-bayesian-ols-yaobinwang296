# Problem Set 6: Bayesian OLS

Yaobin's Comments:

The data from "card.csv" is orgainzed and saved in "data.xlsx", which will be used in the main code.

The main matlab code for this problem set is the .m file named "PS6_YW".

For question 1, I run the OLS regression in the main matlab code file and get the estimates.

For question 2a, I find the posterior dist'n of (beta, sigma_sq) by using the MH algorith
while flat prior is assumed for all parameters.
   -In order to find the posterior, this part of the main code calls for a function named
    "MHstep1_YW". The "MHstep1_YW" function contains the code for the Metropolis-Hastings 
    algorithm. It generates a new proposal and then uses the MH ratio to determine whether 
    the new proposal is kept or not.
       -In order to get the MH ratio, the "MHstep1_YW" function calls for another function
        named "targetdist1_YW". This function returns "log-likelihood" 
        (i.e. log[L(theta)] = log[P(theta)]+log[P(Y|theta)]). In the case of using flat 
        prior for all parameters, this function simply returns log[P(Y|theta)].

For question 2b, I find the posterior dist'n of (beta, sigma_sq) by using the MH algorith
while assuming beta_edu~N(0.06, 0.0007) and flat priors for other parameters.
   -In order to find the posterior, this part of the main code calls for a function named
    "MHstep2_YW". The "MHstep2_YW" function has exactly the same syntax as the 
    "MHstep1_YW" function except it calls for the "targetdist2_YW" function rather than 
    the "targetdist1_YW" function.
       -The only difference between the "targetdist1_YW" and the "targetdist2_YW" 
        function is that "targetdist2" accounts for the prior for beta_edu. So the
        "log-likelihood" that this function returns is slightly different.

To be more clear, the structure of the code is listed below:

Main code "PS6_YW"
    1. OLS
    2. Posterior
        a. flat prior
            -calls for "MHstep1_YW"
                -calls for "targetdist1_YW"
        b. beta_edu~N and flat prior
            -calls for "MHstep2_YW"
                -calls for "targetdist2_YW"

In the end of parts 2a. and 2b., I also report the acceptance rates, which are adjusted 
to be between 20%-25%. And then I draw the histograms for each posterior.
                