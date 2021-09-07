function path_new=MMPP_remove(path,Y,tau_Y,MCMC_param,param)
%% MMPP_remove Creates a new proposal by removing a random jump
%
% Input: path - the last path
%        Y - Pointprocess
%        tau_Y - holding times of the point process
%        MCMC_param - MCMC paramter
%        param -  hyper parameter
%
% Output: path_new -  New metropolis hastings sample
%

%save parameter
f=param.f;
alpha=param.alpha;
q_a=MCMC_param.q_man(2);
q_r=MCMC_param.q_man(3);
q_n=MCMC_param.q_n;
c=path.c;
T=path.T;

%If there are no jumps --> nothing can be removed
if c==0
    path_new=path;
    return
end

%Draw a random jump
del_idx=randsample(1:path.c,1)+1;
%State before deletion
i_old=path.k_i(del_idx);

%save proposal
path_ast=path;
path_ast.t(del_idx)=[];
path_ast.k_i(del_idx)=[];

%Update segment statistics
path_ast.stat.sum_log_y_seg(del_idx)=[];
path_ast.stat.sum_y_seg(del_idx)=[];
path_ast.stat.n_y_seg(del_idx)=[];
path_ast=updateSegStat(path_ast,Y,tau_Y,del_idx-1);
path_ast.c=path.c-1;

%Number of times the old state is used in the proposal
num_i_ast=sum(path_ast.k_i==i_old);
if num_i_ast==0
    %old state is not used anymore
    rem_case=1;
    
    %Delete state and update statistics
    path_ast=removeState(path_ast,i_old);
    path_ast=updateStateStat(path_ast,path_ast.k_i(del_idx-1));
    
else
    %old state is still used
    rem_case=2;
    
    %update statistics
    path_ast=updateStateStat(path_ast,[path.k_i(del_idx-1);path.k_i(del_idx)]);
end


%% LLR
LLR=getLogLikeRatio(path,path_ast);


%% LQR LPR
switch rem_case
    case 1
        
        stat.N=path.stat.n_y_i(i_old);
        stat.T_2=path.stat.sum_y_i(i_old);
        LogPost_i=getLogPostProp(MCMC_param.hyper,stat,path.params(i_old,1),path.params(i_old,2));
        LogPri_i=getLogPri(param.hyper,path.params(i_old,1),path.params(i_old,2));
        
        LQR=log(q_a*q_n*c)+LogPost_i-log(q_r*T);
        LPR=log(alpha+c)-log(f*alpha)-LogPri_i;
    case 2
        
        
        stat.N=path.stat.n_y_seg(del_idx);
        stat.T_1=path.stat.sum_log_y_seg(del_idx);
        stat.T_2=path.stat.sum_y_seg(del_idx);
        [ w ] = postPropWeights(MCMC_param.hyper,stat,path_ast.params(:,1),path_ast.params(:,2)); 
        p_seg=w./sum(w);
        
        num_i=sum(path.k_i==i_old);
        LQR=log(q_a*(1-q_n)*c*p_seg(i_old))-log(q_r*T);
        LPR=log(alpha+c)-log(f*(num_i-1));
end


%% acceptance prob
p_MH=min(1,exp(LQR+(LLR+LPR)*((1/MCMC_param.T)^MCMC_param.annealing)));

%% MH
if rand<p_MH
    path_new=path_ast;
else
    path_new=path;
end

end