function samples=inhomogenousGammaInference(MCMC_param,param,Y)
%Calculate holding times
tau_Y=diff([0;Y]);

%initalize by a draw from the Prior
path=samplePath(param);

%Total number of samples to be saved
N_samples=floor((MCMC_param.N_MCMC-MCMC_param.burnin)/MCMC_param.discard);
samples=cell(N_samples,2);
n_sam=1;

%Calculate the Statistics for the inital path
path=updateSegStat(path,Y,tau_Y,1:path.c+1);
path=updateStateStat(path,1:path.s);

%% Metropolis Hastings!
for n_MCMC=1:MCMC_param.N_MCMC
    
    %echo
    if mod(n_MCMC,floor(N_samples/10))==0
        disp(['MCMC Sample: ' num2str(n_MCMC)]);
    end
    
    %current tempreature
    MCMC_param.T=MCMC_param.T_grid(n_MCMC);
    %current sample number
    MCMC_param.n_sam=n_sam;
    
    %% MH step
    switch randsample(6,1,true,MCMC_param.q_man)
        case 1
            path=MMPP_shift(path,Y,tau_Y,MCMC_param,param);
        case 2
            path=MMPP_add(path,Y,tau_Y,MCMC_param,param);
        case 3
            path=MMPP_remove(path,Y,tau_Y,MCMC_param,param);
        case 4
            path=MMPP_switch(path,Y,tau_Y,MCMC_param,param);
        case 5
            path= MMPP_join(path,Y,tau_Y,MCMC_param,param);
        case 6
            path= MMPP_divide(path,Y,tau_Y,MCMC_param,param);
    end
    
    
    %% Gibbs step
    param.f=gamrnd(param.a_f+path.c,param.b_f/(path.T*param.b_f+1));
    path.params=PostRnd(param.hyper,path.stat);

    %% Save MCMC samples
    if mod(n_MCMC,MCMC_param.discard)==0 && n_MCMC>MCMC_param.burnin
        samples{n_sam,1}=path;
        samples{n_sam,2}=param;
        n_sam=n_sam+1;
    end
    
    
end