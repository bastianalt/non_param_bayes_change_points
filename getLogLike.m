function LL=getLogLike(path)
%% getLogLike Calculates the log-likelihood of a path
%
% Input: path - the path
%
% Output: LL -  Likelihood
%

a=path.params(:,1);
b=path.params(:,2);

LL=-path.stat.n_y_i'*(gammaln(a)+log(b).*a)...
   +path.stat.sum_log_y_i'*(a-1)...
   - path.stat.sum_y_i' *(1./b);

end