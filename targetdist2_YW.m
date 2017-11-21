function logLL = targetdist2_YW(theta,Y,X,n,m)
beta_edu = theta(2);
beta_edu_prior = normpdf(beta_edu, 0.06, 0.007);
Y_pred = X*theta(1:m);
single_PY = zeros(n,1);
for i=1:n
    if theta(m+1) > 0
        single_PY(i) = normpdf(Y(i),Y_pred(i),sqrt(theta(m+1)));
    else
        single_PY(i) = 0;
    end
end
logLL = log(beta_edu_prior)+sum(log(single_PY));
end