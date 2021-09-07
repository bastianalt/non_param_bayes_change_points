function [R] = PostPropRnd(hyper,stat)
%PostPropRnd Samples a random variable R from the posterior probability 
%            for the proposal parameters a and b
%
% Input: hyper - struct containing the hyper parameters of
%        stat - struct containing the sufficient statistics
%
% Output: R -  Sampled random variable
%

lambda=gamrnd(hyper.a_lambda+stat.N,hyper.b_lambda./(hyper.b_lambda*stat.T_2+1));
a=1;
b=1/lambda;
R=[a,b];



