function path_new=MMPP_add(path,Y,tau_Y,MCMC_param,param)
%% MMPP_add Creates a new proposal by adding a random jump
%
% Input: path - the last path
%        Y - Pointprocess
%        tau_Y - holding times of the point process
%        MCMC_param - MCMC paramter
%        param -  hyper parameter
%
% Output: path_new -  New metropolis hastings sample
%


%Save the parameters
f=param.f;
alpha=param.alpha;
q_n=MCMC_param.q_n;
q_a=MCMC_param.q_man(2);
q_r=MCMC_param.q_man(3);
T=param.T;

%Draw a new jump time in the interval [0,T]
t_new=rand*T;

%Add the new jump time to the time vector
path_ast=path;
[path_ast.t,sort_idx]=sort([path.t;t_new]);
%Index of the new jump in the time vector
add_idx=find(sort_idx==length(path.t)+1);

%a segment is added
path_ast.c=path.c+1;

%Update the statistic vectors for the segments
path_ast.stat.n_y_seg=[path.stat.n_y_seg(1:add_idx-1);0;path.stat.n_y_seg(add_idx:end)];
path_ast.stat.sum_log_y_seg=[path.stat.sum_log_y_seg(1:add_idx-1);0;path.stat.sum_log_y_seg(add_idx:end)];
path_ast.stat.sum_y_seg=[path.stat.sum_y_seg(1:add_idx-1);0;path.stat.sum_y_seg(add_idx:end)];
path_ast=updateSegStat(path_ast,Y,tau_Y,[add_idx-1,add_idx]);


create_new_param=rand<q_n || path.s<3; %used to be 3
if  create_new_param
    
    %Draw a new parameter set
    stat.N=path_ast.stat.n_y_seg(add_idx);
    stat.T_1=path_ast.stat.sum_log_y_seg(add_idx);
    stat.T_2=path_ast.stat.sum_y_seg(add_idx);
    new_param=PostPropRnd(MCMC_param.hyper,stat);
    
    %Update the path parameter and statistics
    path_ast.params=[path.params;new_param];
    i=path.s+1;
    path_ast.k_i=[path.k_i(1:add_idx-1);i;path.k_i(add_idx:end)];
    path_ast=updateStateStat(path_ast,[path_ast.k_i(add_idx-1);i]);
    path_ast.s=path.s+1;
    
else
    %Choose an existing parameter set proportinal to the posterior
    stat.N=path_ast.stat.n_y_seg(add_idx);
    stat.T_1=path_ast.stat.sum_log_y_seg(add_idx);
    stat.T_2=path_ast.stat.sum_y_seg(add_idx);
    [ w_ast ] = postPropWeights(MCMC_param.hyper,stat,path.params(:,1),path.params(:,2));
    
    %Update the path parameter and statistics
    i=randsample(1:path.s,1,true,w_ast);
    path_ast.k_i=[path.k_i(1:add_idx-1);i;path.k_i(add_idx:end)];
    path_ast=updateStateStat(path_ast,[path_ast.k_i(add_idx-1);i]);
    
    
end


%% LLR
LLR=getLogLikeRatio(path,path_ast);


%% Q Ratio LPR
if create_new_param
    
    stat.N=path_ast.stat.n_y_i(path.s+1);
    stat.T_2=path_ast.stat.sum_y_i(path.s+1);
    LogPost_i=getLogPostProp(MCMC_param.hyper,stat,path_ast.params(path.s+1,1),path_ast.params(path.s+1,2));
    LogPri_i=getLogPri(param.hyper,path_ast.params(path.s+1,1),path_ast.params(path.s+1,2));
    
    LQR=log(q_r*T)-log(q_a*q_n*(path.c+1))-LogPost_i;
    LPR=log(f*alpha)+LogPri_i-log(alpha+path.c+1);
else
    p_seg_ast=w_ast./sum(w_ast);
    num_ast_i=sum(path_ast.k_i==i);
    
    LQR=log(q_r*T)-log(q_a*(1-q_n)*(path.c+1)*p_seg_ast(i));
    LPR=log(f*(num_ast_i-1))-log(alpha+path.c+1);
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