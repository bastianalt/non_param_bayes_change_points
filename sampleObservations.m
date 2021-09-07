function Y=sampleObservations(path)

%% generate s different poisson processes
Y_cell=cell(path.s,1);
for j=1:path.s
    Y_i=[];
    t_Y=0;
    while t_Y<path.T
        
        %Convert moment and variance to shape and scale
        a=path.params(j,1);
        b=path.params(j,2);
        %         a=2*m/(sqrt(m^2+4*sigma2)-m)+1;
        %         b=(sqrt(m^2+4*sigma2)-m)/2;
        
        h=gamrnd(a,b);
        
        t_Y=t_Y+h;
        Y_i=[Y_i;t_Y]; %#ok<AGROW>
    end
    Y_i=Y_i(1:end-1);
    Y_cell{j}=Y_i;
end

Y=[];
for n=1:path.c+1
    
    tmin=path.t(n);
    tmax=path.t(n+1);
    
    Y_i=Y_cell{path.k_i(n)};
    Y=[Y;Y_i(tmin<=Y_i & Y_i<tmax)]; %#ok<AGROW>
    
end