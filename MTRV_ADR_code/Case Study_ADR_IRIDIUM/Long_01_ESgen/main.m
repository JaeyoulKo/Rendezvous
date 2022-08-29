%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Elementary Solution Generation                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Output
% x_mat: DB for all rendezvous trajectories
%        [index(debris i to j), departure time, transfer time, delV, r_drift orbit, i_drift orbit]
% x_min: minimum-cost rendezvous trajectories
%        [index(debris i to j), departure time, transfer time, delV, r_drift orbit, i_drift orbit,
%         index of initial row in x_mat, index of final row in x_mat]

clear all; close all; clc;
global Q count
count = 0;

load 'target_list.mat'
n = length(Q(:,1));
t_max = 360  % days

x_mat = []; x_min = zeros(n*1000+n,8); l_e = 0;
for i=0:n
    for j=0:n
        if(i~=j)
            disp([i,j])
            x_tmp = []; l_s = l_e+1;
            if(i==0)
                x_tmp = [j,-6,0,0,0,0];
            elseif(j==0)
                x_tmp = [1000*i+j,t_max,0,0,0,0];
            else
                xe = findlocal(i,j,t_max);
                x_tmp = [(1000*i+j)*ones(length(xe(:,1)),1),xe];
            end
            if(~isempty(x_tmp))
                l_e = l_e+length(x_tmp(:,1));
                x_mat = [x_mat; x_tmp];
                x_min(1000*i+j,:) = [x_tmp(1,:),l_s,l_e];
            end
        end
    end
end
% save('x_mat.mat','x_mat');
% save('x_min.mat','x_min');