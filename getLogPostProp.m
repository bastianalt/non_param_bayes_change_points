function   LP=getLogPostProp(hyper,stat,a,b)
%% getLogPostProp Calculates the log-posterior probability for the proposal
%                 parameters a and b
%
% Input: hyper - struct containing the hyper parameters of
%        stat - struct containing the sufficient statistics
%        a - shape parameter of the proposal parameter
%        b - scale parameter of the proposal parameter
%
% Output: LP -  Log posterior probability
%
lambda_i=a./b;

a_post=hyper.a_lambda+stat.N;
b_post=hyper.b_lambda./(stat.T_2*hyper.b_lambda+1);

LP=-gamlike([a_post,b_post],lambda_i);
