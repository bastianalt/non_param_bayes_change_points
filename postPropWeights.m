function [ w ] = postPropWeights(hyper,stat,a,b)
%postPropWeights Calculates the logarithm of the weights that are 
%                proportional to the posterior probability for the proposal
%                 parameters a and b
%
% Input: hyper - struct containing the hyper parameters 
%        stat - struct containing the sufficient statistics
%        a - shape parameter of the proposal parameter
%        b - scale parameter of the proposal parameter
%
% Output: w -  weights for a and b given hyper and stat
%
lambdas=1./a.*b;

a_post=hyper.a_lambda+stat.N;
b_post=hyper.b_lambda./(stat.T_2*hyper.b_lambda+1);

w_log=(a_post-1)*log(lambdas)-lambdas./b_post;
w=exp(w_log-max(w_log));
end