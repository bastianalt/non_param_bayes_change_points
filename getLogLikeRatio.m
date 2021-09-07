
function LLR=getLogLikeRatio(path,path_ast)
%% getLogLikeRatio Calculates the log-likelihood ratio of proposal and old path
%
% Input: path - the old path
%        path_ast - the proposal path
%        tau_Y - holding times of the point process
%        MCMC_param - MCMC paramter
%        param -  hyper parameter
%
% Output: LLR -  Likelihood ratio
%

LL=getLogLike(path);
LL_ast=getLogLike(path_ast);

LLR=LL_ast-LL;
end