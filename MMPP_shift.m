function path_new=MMPP_shift(path,Y,tau_Y,MCMC_param,param) %#ok<INUSD>
%% MMPP_shift Creates a new proposal by shifting a random jump
%
% Input: path - the last path
%        Y - Pointprocess
%        tau_Y - holding times of the point process
%        MCMC_param - MCMC paramter
%        param -  hyper parameter
%
% Output: path_new -  New metropolis hastings sample
%

%If there is no jump in the path, no shift can be made
if path.c==0
    path_new=path;
    return
end

%% Calculate proposal
%draw a random jump
jump_ind=randsample(length(path.t)-2,1)+1;
%Save the time of the jump
t_org=path.t(jump_ind);

%Shift the random jump by a draw from a truncated gaussian
sigma_t=MCMC_param.sigma_t;
pd = makedist('Normal','mu',t_org,'sigma',sigma_t);
%Lower bound
t_min=path.t(jump_ind-1);
%Upper bound
t_max=path.t(jump_ind+1);
pd_tr = truncate(pd,t_min,t_max);
%Shifted jump time
t_ast= random(pd_tr,1,1);

%% Save proposal
path_ast=path;
path_ast.t(jump_ind)=t_ast;

%Recalculate Statistics of the changed segments/ states
path_ast=updateSegStat(path_ast,Y,tau_Y,[jump_ind-1,jump_ind]);
path_ast=updateStateStat(path_ast,path_ast.k_i([jump_ind-1,jump_ind]));

%% LLR
LLR=getLogLikeRatio(path,path_ast);

%% Q ratio
LQR=log(normcdf((t_max-t_org)/sigma_t)-normcdf((t_min-t_org)/sigma_t))...
    -log(normcdf((t_max-t_ast)/sigma_t)-normcdf((t_min-t_ast)/sigma_t));
%% Prior Ratio
LPR=0;

%% acceptance prob
p_MH=min(1,exp(LQR+(LLR+LPR)*((1/MCMC_param.T)^MCMC_param.annealing)));

%% MH
if rand<p_MH
    path_new=path_ast;
else
    path_new=path;
end

end
