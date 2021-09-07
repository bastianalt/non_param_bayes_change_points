function path=updateStateStat(path,update_i_ind)

%Which states have changed?
for j=1:length(update_i_ind)
    path.stat.sum_log_y_i(update_i_ind(j),1)=sum(path.stat.sum_log_y_seg(path.k_i==update_i_ind(j)));
    path.stat.sum_y_i(update_i_ind(j),1)=sum(path.stat.sum_y_seg(path.k_i==update_i_ind(j)));
    path.stat.n_y_i(update_i_ind(j),1)=sum(path.stat.n_y_seg(path.k_i==update_i_ind(j)));
end


end
