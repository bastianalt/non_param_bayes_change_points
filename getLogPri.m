function logPri=getLogPri(hyper,a,b)
%% getLogPostProp Calculates the log-prior probability for the 
%                 parameters a and b
%
% Input: hyper - struct containing the hyper parameters 
%        a - shape parameter of the prior parameter
%        b - scale parameter of the prior parameter
%
% Output: logPri -  Log prior probability
nu_a=hyper.nu_a;
chi_a=hyper.chi_a;
nu_b=hyper.nu_b;
chi_b=hyper.chi_b;

logPri=-gammaln(nu_a)-nu_a*log(chi_a)+(nu_a-1)*log(a)-a/chi_a+ ...
    -gammaln(nu_b)-nu_b*log(chi_b)+(nu_b-1)*log(b)-b/chi_b;

end