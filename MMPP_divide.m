function path_new= MMPP_divide(path,Y,tau_Y,MCMC_param,param)
%% MMPP_divide Creates a new proposal by dividing one state into
%              two random states
%
% Input: path - the last path
%        Y - Pointprocess
%        tau_Y - holding times of the point process
%        MCMC_param - MCMC paramter
%        param -  hyper parameter
%
% Output: path_new -  New metropolis hastings sample
%

%Find the states which are assigned to at least two segments
valid_states=find(sum(path.k_i==1:path.s)>1);
%No division possible
if isempty(valid_states)
    path_new=path;
    return
end

%Save parameter
alpha=param.alpha;
q_j=MCMC_param.q_man(5);
q_d=MCMC_param.q_man(6);


if length(valid_states)==1
    i=valid_states;
else
    i=randsample(valid_states,1);
end

%Sample new parameters for the two new states
new_param=zeros(2,2);

%Draw a new parameter set
stat.N=path.stat.n_y_seg(i);
stat.T_1=path.stat.sum_log_y_seg(i);
stat.T_2=path.stat.sum_y_seg(i);
new_param(1,:)=PostPropRnd(MCMC_param.hyper,stat);
new_param(2,:)=PostPropRnd(MCMC_param.hyper,stat);

%Which segments have to be assigned
toBeAssigned=path.k_i==i;

%assign proportional to likelihood
w_seg_prime=zeros(sum(toBeAssigned),2);

%Choose an existing parameter set proportinal to the posterior
stat.N=path.stat.n_y_seg(toBeAssigned);
stat.T_1= path.stat.sum_log_y_seg(toBeAssigned);
stat.T_2=path.stat.sum_y_seg(toBeAssigned);
w_seg_prime(:,1)= postPropWeights(MCMC_param.hyper,stat,new_param(1,1),new_param(1,2));
w_seg_prime(:,2)= postPropWeights(MCMC_param.hyper,stat,new_param(2,1),new_param(2,2));

w_seg_prime(w_seg_prime==[0,0])=1;
w_seg_prime(isnan(w_seg_prime))=1;
%Number of segments to be assigned
num_i=sum(path.k_i==i);

% Assign all old state values to the new ones
omega=zeros(num_i,1);
for l=1:num_i
    omega(l)=randsample(1:2,1,true,w_seg_prime(l,:));
end

%Check that both states are used
if all(omega==1)
    omega(end)=2;
elseif all(omega==2)
    omega(end)=1;
end

%Save proposal
path_ast=path;
path_ast.k_i(toBeAssigned)=omega+path.s;
path_ast.params=[path.params;new_param];
path_ast=updateStateStat(path_ast,[path.s+1;path.s+2]);
path_ast.s=path.s+2;
path_ast=removeState(path_ast,i);

%New states
j1=path_ast.s-1;
j2=path_ast.s;

%% LLR
LLR=getLogLikeRatio(path,path_ast);

%% Psi
s_1=sum(valid_states);
num_j1_ast=sum(path_ast.k_i==j1);
num_j2_ast=sum(path_ast.k_i==j2);
num_i=sum(path.k_i==i);

p_seg_prime=w_seg_prime./sum(w_seg_prime,2);
omega_prime=[ones(num_i,1)==omega, 2*ones(num_i,1)==omega];
log_p_par=sum(log(p_seg_prime(omega_prime)));



stat.N=path.stat.n_y_i(i);
stat.T_2=path.stat.sum_y_i(i);
LogPost_i=getLogPostProp(MCMC_param.hyper,stat,path.params(i,1),path.params(i,2));
LogPri_i=getLogPri(param.hyper,path.params(i,1),path.params(i,2));

stat.N=path_ast.stat.n_y_i(j1);
stat.T_2=path_ast.stat.sum_y_i(j1);
LogPost_j1=getLogPostProp(MCMC_param.hyper,stat,path_ast.params(j1,1),path_ast.params(j1,2));
LogPri_j1=getLogPri(param.hyper,path_ast.params(j1,1),path_ast.params(j1,2));

stat.N=path_ast.stat.n_y_i(j2);
stat.T_2=path_ast.stat.sum_y_i(j2);
LogPost_j2=getLogPostProp(MCMC_param.hyper,stat,path_ast.params(j2,1),path_ast.params(j1,2));
LogPri_j2=getLogPri(param.hyper,path_ast.params(j2,1),path_ast.params(j2,2));

LQR=log(q_j*s_1)+LogPost_i...
    -log_p_par-log(q_d*path_ast.s*(path_ast.s-1))-LogPost_j1-LogPost_j2;
LPR=log(alpha*factorial(num_j1_ast-1)*factorial(num_j2_ast-1))+LogPri_j1+LogPri_j2-log(factorial(num_i-1))-LogPri_i;

%% acceptance prob
p_MH=min(1,exp(LQR+(LLR+LPR)*((1/MCMC_param.T)^MCMC_param.annealing)));

%% MH
if rand<p_MH
    path_new=path_ast;
else
    path_new=path;
end

end