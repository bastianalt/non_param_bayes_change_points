%%%%%%%%%%%%%% Example %%%%%%%%%%%%%%

%%%%%%%%%%%%%% In silico Data %%%%%%%%%%%%%%
%% Generate GT Data
%pseudo counts controlling the shape parameter
param_GT.hyper.nu_a=100;
%pseudo sufficient statistics of the shape parameter
param_GT.hyper.chi_a=1.5e-2;

%pseudo counts controlling the scale parameter
param_GT.hyper.nu_b=5;
%pseudo sufficient statistics of the scale parameter
param_GT.hyper.chi_b=.1*.5;

%Rate parameter of the jump times
param_GT.f=0.02;
%conncentration parameter CRP
param_GT.alpha=3;
%Time length
param_GT.T=1000;

%% Ground_truth
GT = samplePath(param_GT);

%% Generate the observation data
Y=sampleObservations(GT);

%%%%%%%%%%%%%% Inference %%%%%%%%%%%%%%
%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%
%% MCMC Parameters
%Number of burnin samples
MCMC_param.burnin=1;
%Total number of samples
MCMC_param.N_MCMC=1e4+MCMC_param.burnin;
%Just use every "MCMC_param.discard"'s sample
MCMC_param.discard=1;

%Random Walk Parameters
MCMC_param.sigma_t=1e6; %Variance for shifting of jumps
MCMC_param.q_n=.1;  %Probability of creating a new value
MCMC_param.q_man=ones(6,1)/6; %The probabilities for the proposal actions (shift,add,remove,switch,join and divide)
% hyper parameter of the sampler of the propsal parameters
MCMC_param.hyper.a_lambda=1e-2;
MCMC_param.hyper.b_lambda=1e2;

%Annealing and the temperture grid
MCMC_param.annealing=true;
MCMC_param.T_grid=logspace(0,-10,MCMC_param.N_MCMC);


%% Hyper parameters
%Shape and Scale
param.hyper.nu_a=1e-2;
param.hyper.chi_a=1*1e2;
param.hyper.nu_b=1e-2;
param.hyper.chi_b=1e2;

%Jump rate
param.a_f=1e-2;
param.b_f=1e2;
%It is useful to start with a rate >0 to ensure a faster convergence
param.f=0.01;

%Concentration parameter
param.alpha=param_GT.alpha;

%Sequence length
param.T=max(Y);


%%%%%%%%%%%%%% MCMC Sampling %%%%%%%%%%%%%%
samples=inhomogenousGammaInference(MCMC_param,param,Y);

%%%%%%%%%%%%%% Plot posterior %%%%%%%%%%%%%%
%Resample Y for plots
Y_resample=resample(Y,1,cite(length(Y)/200));

%% Posterior mean and quantiles vs GT
dt=1;
t=0:dt:param.T-dt;
N_samples=size(samples,1);
a_sam=zeros(length(t),N_samples);
b_sam=zeros(length(t),N_samples);

for n_sam=1:N_samples
    for n=1:length(t)
        idx=find(t(n)>=samples{n_sam,1}.t,1,'last');
        a_sam(n,n_sam)=samples{n_sam,1}.params(samples{n_sam,1}.k_i(idx),1);
        b_sam(n,n_sam)=samples{n_sam,1}.params(samples{n_sam,1}.k_i(idx),2);
    end
end

figure,
subplot(2,1,1)
shadedErrorBar(t,mean(a_sam,2),[quantile(a_sam,.95,2)'-mean(a_sam,2)';mean(a_sam,2)'-quantile(a_sam,.05,2)']')
hold on
myPlotA(GT)
ylabel('a')
xlim([0,1000]);
ylim([0,3]);
hold on
for i=1:length(Y_resample)
    plot([Y_resample(i),Y_resample(i)],[0,.5],'bla')
end

subplot(2,1,2)
shadedErrorBar(t,mean(b_sam,2),[quantile(b_sam,.95,2)'-mean(b_sam,2)';mean(b_sam,2)'-quantile(b_sam,.05,2)']')
hold on
myPlotB(GT)
xlim([0,1000]);
xlabel('time')
ylabel('b')
hold on
for i=1:length(Y_resample)
    plot([Y_resample(i),Y_resample(i)],[.8,1],'bla')
end

