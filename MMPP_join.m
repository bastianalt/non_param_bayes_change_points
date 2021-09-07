function path_new= MMPP_join(path,Y,tau_Y,MCMC_param,param)
%% MMPP_join Creates a new proposal by joining two random states
%
% Input: path - the last path
%        Y - Pointprocess
%        tau_Y - holding times of the point process
%        MCMC_param - MCMC paramter
%        param -  hyper parameter
%
% Output: path_new -  New metropolis hastings sample
%

%If there is just one state nothing can be joined
if path.s==1
    path_new=path;
    return
end

%Save parameter
alpha=param.alpha;
q_j=MCMC_param.q_man(5);
q_d=MCMC_param.q_man(6);

%Save proposal
path_ast=path;

%Draw two different random states (without replacement)
i=randsample(1:path.s,2,false);
i1=i(1);
i2=i(2);

%Generate new parameter by sampling from the posterior
assigned2j=path.k_i==i1 | path.k_i==i2;

stat.N=sum(path_ast.stat.n_y_seg(assigned2j));
stat.T_1=sum(path_ast.stat.sum_log_y_seg(assigned2j));
stat.T_2=sum(path_ast.stat.sum_y_seg(assigned2j));
new_param=PostPropRnd(MCMC_param.hyper,stat);

%Generate proposal by assigning all segments to the new parameter
path_ast.k_i(assigned2j)=i1;
path_ast.params(i1,:)=new_param;
path_ast=removeState(path_ast,i2);

j=find(path_ast.params(:,1)==new_param(1));
j=j(1);
path_ast=updateStateStat(path_ast,j);

%% LLR
LLR=getLogLikeRatio(path,path_ast);

%% Psi
num_i1=sum(path.k_i==i1);
num_i2=sum(path.k_i==i2);
num_j_ast=sum(path_ast.k_i==j);

%Choose an existing parameter set proportinal to the posterior



w_seg_prime=zeros(sum(assigned2j),2);
stat.N=path_ast.stat.n_y_seg(assigned2j);
stat.T_1=path_ast.stat.sum_log_y_seg(assigned2j);
stat.T_2=path_ast.stat.sum_y_seg(assigned2j);
w_seg_prime(:,1)= postPropWeights(MCMC_param.hyper,stat,path.params(i1,1),path.params(i1,2));
w_seg_prime(:,2)= postPropWeights(MCMC_param.hyper,stat,path.params(i2,1),path.params(i2,2));

w_seg_prime(w_seg_prime==[0,0])=1;
p_seg_prime=w_seg_prime./sum(w_seg_prime,2);

omega_prime=[path.k_i(assigned2j)==i1,path.k_i(assigned2j)==i2];
log_p_par=sum(log(p_seg_prime(omega_prime)));

s_1=sum(sum(path_ast.k_i==1:path_ast.s)>1);



stat.N=path_ast.stat.n_y_i(j);
stat.T_2=path_ast.stat.sum_y_i(j);
LogPost_j=getLogPostProp(MCMC_param.hyper,stat,path_ast.params(j,1),path_ast.params(j,2));
LogPri_j=getLogPri(param.hyper,path_ast.params(j,1),path_ast.params(j,2));

stat.N=path.stat.n_y_i(i1);
stat.T_2=path.stat.sum_y_i(i1);
LogPost_i1=getLogPostProp(MCMC_param.hyper,stat,path.params(i1,1),path.params(i1,2));
LogPri_i1=getLogPri(param.hyper,path.params(i1,1),path.params(i1,2));

stat.N=path.stat.n_y_i(i2);
stat.T_2=path.stat.sum_y_i(i2);
LogPost_i2=getLogPostProp(MCMC_param.hyper,stat,path.params(i2,1),path.params(i2,2));
LogPri_i2=getLogPri(param.hyper,path.params(i2,1),path.params(i2,2));


LQR=log_p_par+log(q_d*path.s*(path.s-1))+LogPost_i1+LogPost_i2...
    -log(q_j *s_1)-LogPost_j;
LPR=log(factorial(num_j_ast-1))+LogPri_j-log(alpha*factorial(num_i1-1)*factorial(num_i2-1))-LogPri_i1-LogPri_i2;

%% acceptance prob
p_MH=min(1,exp(LQR+(LLR+LPR)*((1/MCMC_param.T)^MCMC_param.annealing)));

%% MH
if rand<p_MH
    path_new=path_ast;
else
    path_new=path;
end

end