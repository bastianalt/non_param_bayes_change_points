function path=updateSegStat(path,Y,tau_Y,update_seg_ind)

for l=1:length(update_seg_ind)
    
    idx_seg=update_seg_ind(l);
    idx_n=path.t(idx_seg) <= Y & Y <= path.t(idx_seg+1);
    
    path.stat.sum_log_y_seg(idx_seg)= sum(log(tau_Y(idx_n)));
    path.stat.sum_y_seg(idx_seg)=sum(tau_Y(idx_n));
    path.stat.n_y_seg(idx_seg)=length(Y(idx_n));

end 

end
