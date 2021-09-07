function myPlot(path)
%% Plots the paths for shape a and scale b

 hold on
t=path(1).t;
param_i=path.params(path.k_i,:);
%Plot
dt=.1;
c=size(param_i,1)-1;
yyaxis left
for i=1:c+1
    t_i=t(i):dt:t(i+1);
    
    
    plot(t_i,repmat(param_i(i,1),length(t_i),1),'blu-');
    if i<=c

        plot([t(i+1),t(i+1)],[param_i(i,1),param_i(i+1,1)],'blu:')
    end
end
ylabel('Shape a')
yyaxis right
for i=1:c+1
    t_i=t(i):dt:t(i+1);
    
    plot(t_i,repmat(param_i(i,2),length(t_i),1),'red-');
    
    if i<=c    
        plot([t(i+1),t(i+1)],[param_i(i,2),param_i(i+1,2)],'red:')
        
    end
end
ylabel('Scale b')
hold off
end