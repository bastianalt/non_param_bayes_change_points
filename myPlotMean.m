function myPlotMean(path)
%% Plots the mean of the path

 hold on
    t=path.t;
param_i=path.params(path.k_i,:);
%Plot
dt=.1;
c=size(param_i,1)-1;

for i=1:c+1
    t_i=t(i):dt:t(i+1);
    
    
    plot(t_i,repmat((1./(param_i(i,1)*param_i(i,2))),length(t_i),1),'red--');
   
    
    if i<=c

        plot([t(i+1),t(i+1)],[(1./(param_i(i,1)*param_i(i,2))),(1./(param_i(i+1,1)*param_i(i+1,2)))],'red:')
    end
end
hold off
end