function myPlotVar(path)
%% Plots the variance of the path

 hold on
    t=path.t;
param_i=path.params(path.k_i,:);
%Plot
dt=.1;
c=size(param_i,1)-1;

for i=1:c+1
    t_i=t(i):dt:t(i+1);
    
    
    plot(t_i,repmat((param_i(i,1)*param_i(i,2)^2),length(t_i),1),'blu-');
   
    
    if i<=c

        plot([t(i+1),t(i+1)],[(param_i(i,1)*param_i(i,2)^2),(param_i(i+1,1)*param_i(i+1,2)^2)],'blu:')
    end
end
hold off
end