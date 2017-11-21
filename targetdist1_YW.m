function logLL = targetdist1_YW(theta,Y,X,n,m)
Y_pred = X*theta(1:m);
single_LL = zeros(n,1);
for i=1:n
    if theta(m+1) > 0
        single_LL(i) = normpdf(Y(i),Y_pred(i),sqrt(theta(m+1)));
    else
        single_LL(i) = 0;
    end
end
logLL = sum(log(single_LL));
end