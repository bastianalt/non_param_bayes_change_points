function path=removeState(path,i)
%% REMOVESTATE removes a state i from a path
%
% Input: path - path
%        i - state to be deleted
%
% Output: path - updated path
%

%delete parameter
path.params(i,:)=[];
%delete statistics
path.stat.sum_log_y_i(i)=[];
path.stat.sum_y_i(i)=[];
path.stat.n_y_i(i)=[];

path.k_i(path.k_i>i)=path.k_i(path.k_i>i)-1;
path.s=path.s-1;

end