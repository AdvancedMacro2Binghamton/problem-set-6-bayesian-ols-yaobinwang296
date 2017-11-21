function [theta_1,a] = MHstep2_YW(theta_0,sig,Y,X,n,m)
theta_p = mvnrnd(theta_0,sig)'; % generate candidate from Gaussian
accprob = exp(targetdist2_YW(theta_p,Y,X,n,m)-targetdist2_YW(theta_0,Y,X,n,m)); % acceptance probability
u = rand; % uniform random number
if u <= accprob % if accepted
    theta_1 = theta_p; % new point is the candidate
    a = 1; % note the acceptance
else % if rejected
    theta_1 = theta_0; % new point is the same as the old one
    a = 0; % note the rejection
end
end