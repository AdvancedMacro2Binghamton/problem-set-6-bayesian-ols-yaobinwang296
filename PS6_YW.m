%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ECON 634 Macro II
%%% Problem Set 6
%%% Bayesian OLS
%%% Yaobin Wang
%%% 11/21/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Question 1 OLS: log wage on education, experience, and dummy variables
%%%                 for SMSA status, race ("black"), and region ("south").
clear all; close all; clc;
cd 'D:\Binghamton University--SUNY\5. Fall 2017\ECON 634\Github\problem-set-6-bayesian-ols-yaobinwang296'
data = xlsread('data.xlsx');
[n,m] = size(data); % n = number of observation, m = number of regressors
Y = data(:,1); % ln_wage as dependent variable
X = [ones(n,1) data(:,2:end)]; % matrix of intercept and independent variables
% OLS results
beta = (X'*X)\(X'*Y); % estimated coefficients
sigma_sq = (1/(n-m-1))*(Y-X*beta)'*(Y-X*beta); % estimated variance of the residuals
var_beta = sigma_sq*diag(inv(X'*X));
se_beta = sqrt(var_beta); % standard errors of beta


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Question 2 find the posterior dist'n of (beta, sigma_sq) by usiing the
%%%            Metropolis-Hastings algorithm.
k = 0.11; % ajust k such that the acceptance rate is between 20%-25%
sigma_sq_v = k*[[bsxfun(@times,var_beta,eye(m)) zeros(m,1)];...
    [zeros(1,m) ((2/(n-m-1))*(sigma_sq)^2)]]; 
nsamp = 1e4; % number of samples (iterations)
burnin = 0.1*nsamp; % number of burn-in iterations
lag = 1; % could allow some iterations between successive samples


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% a. flat prior for all parameters;
theta = [beta;sigma_sq]; % initial points
% storage
THETA1 = zeros(m+1,nsamp);
acc1_b = [0,0];
acc1 = [0,0];
% MH routine
for i = 1:burnin
       [theta,a] = MHstep1_YW(theta,sigma_sq_v,Y,X,n,m);
       acc1_b = acc1_b+[a 1]; % track accept-reject status
end
for i = 1:nsamp
    for j = 1:lag
        [theta,a] = MHstep1_YW(theta,sigma_sq_v,Y,X,n,m);
        acc1 = acc1 + [a 1]; % track accept-reject status
    end
    THETA1(:,i) = theta; % store the i-th sample
    disp(i);
end
acc1_rt = (acc1(1)/acc1(2))*100; % acceptance rate
disp(['The acceptance rate is ', num2str(acc1_rt),'%']);

% summary stat of posteriors
mean11 = mean(THETA1(1,:));
var11 = var(THETA1(1,:));
mean11_label=['Mean = ',num2str(mean11)];
var11_label=['Variance = ',num2str(var11)];

mean12 = mean(THETA1(2,:));
var12 = var(THETA1(2,:));
mean12_label=['Mean = ',num2str(mean12)];
var12_label=['Variance = ',num2str(var12)];

mean13 = mean(THETA1(3,:));
var13 = var(THETA1(3,:));
mean13_label=['Mean = ',num2str(mean13)];
var13_label=['Variance = ',num2str(var13)];

mean14 = mean(THETA1(4,:));
var14 = var(THETA1(4,:));
mean14_label=['Mean = ',num2str(mean14)];
var14_label=['Variance = ',num2str(var14)];

mean15 = mean(THETA1(5,:));
var15 = var(THETA1(5,:));
mean15_label=['Mean = ',num2str(mean15)];
var15_label=['Variance = ',num2str(var15)];

mean16 = mean(THETA1(6,:));
var16 = var(THETA1(6,:));
mean16_label=['Mean = ',num2str(mean16)];
var16_label=['Variance = ',num2str(var16)];

mean17 = mean(THETA1(7,:));
var17 = var(THETA1(7,:));
mean17_label=['Mean = ',num2str(mean17)];
var17_label=['Variance = ',num2str(var17)];

% Plot histograms of the posterior approximations
dim = [.15 .6 .1 .3];
figure
histfit(THETA1(1,:),50,'kernel')
annotation('textbox',dim,'String',{mean11_label, var11_label},'FitBoxToText','on')
ylabel('Density')
title('\beta_{cons}')

figure
histfit(THETA1(2,:),50,'kernel')
annotation('textbox',dim,'String',{mean12_label, var12_label},'FitBoxToText','on')
ylabel('Density')
title('\beta_{edu}')

figure
histfit(THETA1(3,:),50,'kernel')
annotation('textbox',dim,'String',{mean13_label, var13_label},'FitBoxToText','on')
ylabel('Density')
title('\beta_{exp}')

figure
histfit(THETA1(4,:),50,'kernel')
annotation('textbox',dim,'String',{mean14_label, var14_label},'FitBoxToText','on')
ylabel('Density')
title('\beta_{SMSA}')

figure
histfit(THETA1(5,:),50,'kernel')
annotation('textbox',dim,'String',{mean15_label, var15_label},'FitBoxToText','on')
ylabel('Density')
title('\beta_{race}')

figure
histfit(THETA1(6,:),50,'kernel')
annotation('textbox',dim,'String',{mean16_label, var16_label},'FitBoxToText','on')
ylabel('Density')
title('\beta_{south}')

figure
histfit(THETA1(7,:),50,'kernel')
annotation('textbox',dim,'String',{mean17_label, var17_label},'FitBoxToText','on')
ylabel('Density')
title('\sigma^2_{\epsilon}')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% b. prior: beta_edu~N[0.06, 2*se_beta(2)], flat prior for others;
theta = [beta;sigma_sq]; % initial points
% storage
THETA2 = zeros(m+1,nsamp);
acc2_b = [0,0];
acc2 = [0,0];
% MH routine
for i = 1:burnin
       [theta,a] = MHstep2_YW(theta,sigma_sq_v,Y,X,n,m);
       acc2_b = acc2_b+[a 1]; % track accept-reject status
end
for i = 1:nsamp
    for j = 1:lag
        [theta,a] = MHstep2_YW(theta,sigma_sq_v,Y,X,n,m);
        acc2 = acc2 + [a 1]; % track accept-reject status
    end
    THETA2(:,i) = theta; % store the i-th sample
    disp(i);
end
acc2_rt = (acc2(1)/acc2(2))*100; % acceptance rate
disp(['The acceptance rate is ', num2str(acc2_rt),'%']);

% summary stat of posteriors
mean21 = mean(THETA2(1,:));
var21 = var(THETA2(1,:));
mean21_label=['Mean = ',num2str(mean21)];
var21_label=['Variance = ',num2str(var21)];

mean22 = mean(THETA2(2,:));
var22 = var(THETA2(2,:));
mean22_label=['Mean = ',num2str(mean22)];
var22_label=['Variance = ',num2str(var22)];

mean23 = mean(THETA2(3,:));
var23 = var(THETA2(3,:));
mean23_label=['Mean = ',num2str(mean23)];
var23_label=['Variance = ',num2str(var23)];

mean24 = mean(THETA2(4,:));
var24 = var(THETA2(4,:));
mean24_label=['Mean = ',num2str(mean24)];
var24_label=['Variance = ',num2str(var24)];

mean25 = mean(THETA2(5,:));
var25 = var(THETA2(5,:));
mean25_label=['Mean = ',num2str(mean25)];
var25_label=['Variance = ',num2str(var25)];

mean26 = mean(THETA2(6,:));
var26 = var(THETA2(6,:));
mean26_label=['Mean = ',num2str(mean26)];
var26_label=['Variance = ',num2str(var26)];

mean27 = mean(THETA2(7,:));
var27 = var(THETA2(7,:));
mean27_label=['Mean = ',num2str(mean27)];
var27_label=['Variance = ',num2str(var27)];

% Plot histograms of the posterior approximations
figure
histfit(THETA2(1,:),50,'kernel')
annotation('textbox',dim,'String',{mean21_label, var21_label},'FitBoxToText','on')
ylabel('Density')
title('\beta_{cons}')

figure
histfit(THETA2(2,:),50,'kernel')
annotation('textbox',dim,'String',{mean22_label, var22_label},'FitBoxToText','on')
ylabel('Density')
title('\beta_{edu}')

figure
histfit(THETA2(3,:),50,'kernel')
annotation('textbox',dim,'String',{mean23_label, var23_label},'FitBoxToText','on')
ylabel('Density')
title('\beta_{exp}')

figure
histfit(THETA2(4,:),50,'kernel')
annotation('textbox',dim,'String',{mean24_label, var24_label},'FitBoxToText','on')
ylabel('Density')
title('\beta_{SMSA}')

figure
histfit(THETA2(5,:),50,'kernel')
annotation('textbox',dim,'String',{mean25_label, var25_label},'FitBoxToText','on')
ylabel('Density')
title('\beta_{race}')

figure
histfit(THETA2(6,:),50,'kernel')
annotation('textbox',dim,'String',{mean26_label, var26_label},'FitBoxToText','on')
ylabel('Density')
title('\beta_{south}')

figure
histfit(THETA2(7,:),50,'kernel')
annotation('textbox',dim,'String',{mean27_label, var27_label},'FitBoxToText','on')
ylabel('Density')
title('\sigma^2_{\epsilon}')
