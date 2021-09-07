function path_new=MMPP_switch(path,Y,tau_Y,MCMC_param,param)
%% MMPP_switch Creates a new proposal by switching a random segment to a
%              different state
%
% Input: path - the last path
%        Y - Point process
%        tau_Y - holding times of the point process
%        MCMC_param - MCMC paramter
%        param -  hyper parameter
%
% Output: path_new -  New metropolis hastings sample
%

%save parameter
alpha=param.alpha;
q_n=MCMC_param.q_n;

%Save proposal
path_ast=path;

%draw a random segment
sw_idx=randsample(1:path.c+1,1);

%old state
i=path.k_i(sw_idx);

create_new_param=rand<q_n;
if create_new_param
    
    %Draw a new parameter set
    stat.N=path_ast.stat.n_y_seg(sw_idx);
    stat.T_1=path_ast.stat.sum_log_y_seg(sw_idx);
    stat.T_2=path_ast.stat.sum_y_seg(sw_idx);
    new_param=PostPropRnd(MCMC_param.hyper,stat);
    
    %Update the path parameter and statistics
    path_ast.params=[path.params;new_param];
    j=path.s+1;
    path_ast.k_i(sw_idx)=j;
    path_ast.s=path.s+1;
    
    path_ast=updateStateStat(path_ast,[i;j]);
    
else
    %Choose an existing parameter set proportinal to the posterior
    stat.N=path_ast.stat.n_y_seg(sw_idx);
    stat.T_1=path_ast.stat.sum_log_y_seg(sw_idx);
    stat.T_2=path_ast.stat.sum_y_seg(sw_idx);
    [ w_ast ] = postPropWeights(MCMC_param.hyper,stat,path.params(:,1),path.params(:,2));
    
    j=randsample(1:path.s,1,true,w_ast);
    
    %Update the path parameter and statistics
    path_ast.k_i(sw_idx)=j;
    
    path_ast=updateStateStat(path_ast,[i;j]);
    
    
end

% Check if the old state is still used in the proposal
if sum(path_ast.k_i==i)==0
    still_used=false;
else
    still_used=true;
end

%% LLR
LLR=getLogLikeRatio(path,path_ast);

%% Psi
num_i=sum(path.k_i==i);
num_ast_j=sum(path_ast.k_i==j);
if still_used  && ~create_new_param
    
    p_seg_ast=w_ast./sum(w_ast);
    
    
    stat.N=path.stat.n_y_seg(sw_idx);
    stat.T_1=path.stat.sum_log_y_seg(sw_idx);
    stat.T_2=path.stat.sum_y_seg(sw_idx);
    [ w ] = postPropWeights(MCMC_param.hyper,stat,path_ast.params(:,1),path_ast.params(:,2));
    
    p_seg=w./sum(w);
    
    LQR=log(p_seg(i))-log(p_seg_ast(j));
    LPR=log(num_ast_j-1)-log(num_i-1);
    
elseif still_used && create_new_param
    
    stat.N=path.stat.n_y_seg(sw_idx);
    stat.T_1=path.stat.sum_log_y_seg(sw_idx);
    stat.T_2=path.stat.sum_y_seg(sw_idx);
    [ w ] = postPropWeights(MCMC_param.hyper,stat,path_ast.params(:,1),path_ast.params(:,2));
    
    p_seg=w./sum(w);
    
    stat.N=path_ast.stat.n_y_i(j);
    stat.T_2=path_ast.stat.sum_y_i(j);
    LogPost_j=getLogPostProp(MCMC_param.hyper,stat,path_ast.params(j,1),path_ast.params(j,2));
    LogPri_j=getLogPri(param.hyper,path_ast.params(j,1),path_ast.params(j,2));
    
    LQR=log((1-q_n)*p_seg(i))-log(q_n)-LogPost_j;
    LPR=log(alpha)+LogPri_j-log(num_i-1);
    
elseif ~still_used && ~create_new_param
    
    p_seg_ast=w_ast./sum(w_ast);
    
    stat.N=path.stat.n_y_i(i);
    stat.T_2=path.stat.sum_y_i(i);
    LogPost_i=getLogPostProp(MCMC_param.hyper,stat,path.params(i,1),path.params(i,2));
    LogPri_i=getLogPri(param.hyper,path.params(i,1),path.params(i,2));
    
    LQR=log(q_n)+LogPost_i-log((1-q_n)*p_seg_ast(j));
    LPR=log(num_ast_j-1)-log(alpha)-LogPri_i;
    
elseif  ~still_used && create_new_param
    
    stat.N=path.stat.n_y_i(i);
    stat.T_2=path.stat.sum_y_i(i);
    LogPost_i=getLogPostProp(MCMC_param.hyper,stat,path.params(i,1),path.params(i,2));
    LogPri_i=getLogPri(param.hyper,path.params(i,1),path.params(i,2));
    
    stat.N=path_ast.stat.n_y_i(j);
    stat.T_2=path_ast.stat.sum_y_i(j);
    LogPost_j=getLogPostProp(MCMC_param.hyper,stat,path_ast.params(j,1),path_ast.params(j,2));
    LogPri_j=getLogPri(param.hyper,path_ast.params(j,1),path_ast.params(j,2));
    
    LQR=LogPost_i-LogPost_j;
    LPR=LogPri_j-LogPri_i;
end

%% acceptance prob
p_MH=min(1,exp(LQR+(LLR+LPR)*((1/MCMC_param.T)^MCMC_param.annealing)));

%% MH
if rand<p_MH
    % Delete old state if no used anymore
    if ~still_used
        path_ast=removeState(path_ast,i);
    end
    
    path_new=path_ast;
else
    path_new=path;
end


end
