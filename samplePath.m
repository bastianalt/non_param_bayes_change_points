function [path] = samplePath(param)
%SAMPLEDATA Creates a sample from the prior of a path lambda_0_T
%           from a MMPPs with infinite number of states
%
%
%   Input: f -rate parameter for the jumps
%          alpha - concentration parameter of the dirichlet process
%          T - length of the data samples
%
%   Output: lambda_0_T- struct
%           lambda_0_T.t - jump times
%           lambda_0_T.tau - soujourn times
%           lambda_0_T.lambdas - all rate values
%           lambda_0_T.lambda_i - vector of corresponding rates between the jumps
%           lambda_0_T.lambda_ind - index of the rate in the lambdas vector
%
% f=0.02;
% alpha=3;
% T=1000;

f=param.f;
alpha=param.alpha;

T=param.T;

%Initalize output
path=struct('t',[],'params',[],'k_i',[]);

%Sample the jumps
t=0;
tau=[];
while t<T
    newtau=exprnd(1/f);
    t=t+newtau;
    tau=[tau;newtau]; %#ok<AGROW>
end
%Delete last soujourn time because it is larger than T
tau=[tau(1:end-1);T-sum(tau(1:end-1))];

%Number of jumps
c=length(tau)-1;
%Time vector
t=[0;cumsum(tau)];

%Lambda values for the jumps
params=[];
k_i=zeros(c+1,1);

for i=0:c
    p_new=alpha/(alpha+i);
    
    if rand<p_new || i<2
        %sample from the prior
        new_param=priorRnd(param.hyper);
        params=[params;new_param]; %#ok<AGROW>
        k_i(i+1)=size(params,1);
    else     
        [~,k_i(i+1)]=datasample(1:c,1,1);
    end
end

% %delete obsolete jumps
% t(find(diff(k_i)==0)+2)=[];
% k_i(find(diff(k_i)==0)+1)=[];


path.t=t;
path.params=params;
% path.param_i=param_i;
path.k_i=k_i;
path.T=T;
path.c=size(k_i,1)-1;
path.s=size(params,1);

path.stat.sum_log_y_seg=zeros(path.c+1 ,1);
path.stat.sum_y_seg=zeros(path.c+1 ,1);
path.stat.n_y_seg=zeros(path.c+1 ,1);

path.stat.sum_log_y_i=zeros(path.s,1);
path.stat.sum_y_i=zeros(path.s,1);
path.stat.n_y_i=zeros(path.s,1);

end



