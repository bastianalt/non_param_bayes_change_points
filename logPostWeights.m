function [ log_w ] = logPostWeights(hyper,stat,a,b)
%LOGPOSTWEIGHTS Calculates the log of a weight function that is
%proportional to the posterior (exp(log_w)/Z=p(a,b| data))
%
% Input: hyper - struct containing the hyper parameters 
%        stat - struct containing the sufficient statistics
%        a - shape parameter of the proposal parameter
%        b - scale parameter of the proposal parameter
%
% Output: log_w -  Log of the posterior weights
%


log_w=-stat.N*gammaln(a)-stat.N*a*log(b)+(a-1)*stat.T_1-stat.T_2/b+...
    (hyper.nu_a-1)*log(a)-a/hyper.chi_a+(hyper.nu_b-1)*log(b)-b/hyper.chi_b;

end

