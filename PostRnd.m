function [R] = PostRnd(hyper,stat)
%POSTRND Samples new values for the shape and scale parameters given a path
%        and the data
%
% Input: hyper - struct containing the hyper parameters 
%        stat - struct containing the sufficient statistics
%
% Output: R -  new set of sampled parameters [a,b]
%


%MH parameter
burnin=100;
sigma=.1;

%% Draw for each parameter set
s=size(stat.n_y_i,1);
R=zeros(s,2);

for i=1:s
    
    %Current set of statistics
    stat_i.N=stat.n_y_i(i);
    stat_i.T_1=stat.sum_log_y_i(i);
    stat_i.T_2=stat.sum_y_i(i);
    
    %% Initalize
    b_opt=@(a) -(stat_i.N*a*hyper.chi_b-hyper.nu_b*hyper.chi_b+hyper.chi_b)/2+...
        sqrt((stat_i.N*a*hyper.chi_b-hyper.nu_b*hyper.chi_b+hyper.chi_b)^2/4+hyper.chi_b*stat_i.T_2);
    
    %     fun=@(a) gammaln(a)*stat_i.N+stat_i.N*a*log(b_opt(a))-(a-1)*stat_i.T_1+stat_i.T_2/b_opt(a)...
    %         -(hyper.nu_a-1)*log(a)+a/hyper.chi_a-(hyper.nu_b-1)*log(b_opt(a))+b_opt(a)/hyper.chi_b;
    %     a=fminsearch(fun,1); 
    
    a=1;
    b=b_opt(a);  
    if b==0
        b=exprnd(.1);
    end

    %% Random walk!
    for t=1:burnin
        %% MH
        a_ast=exp(randn*sigma)*a;
        b_ast=exp(randn*sigma)*b;
        
        %aceptance probability
        p_MH=min(1,...
            exp( logPostWeights(hyper,stat_i,a_ast,b_ast) - logPostWeights(hyper,stat_i,a,b) )...
            );
        
        if rand<p_MH
            a=a_ast;
            b=b_ast;
        end

    end
    R(i,:)=[a,b];
    
end



