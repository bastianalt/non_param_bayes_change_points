function R=priorRnd(hyper)
%POSTRND Samples new values for the shape and scale parameters from the
%        prior
%
% Input: hyper - struct containing the hyper parameters 
%
% Output: R -  new set of sampled parameters [a,b]
%
a=gamrnd(hyper.nu_a,hyper.chi_a);
b=gamrnd(hyper.nu_b,hyper.chi_b);
R=[a,b];