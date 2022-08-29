%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Find Locally Optimal Departure/Transfer times b/w i and j             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xe = findlocal(i,j,tmax)
global targeti targetj t_max ttr_lb ttr_ub

%% Input parameters
ttr_lb = 0;     % lb for transfer time
ttr_ub = 30;    % ub for transfer time
del_td = 30;    % size of grid - departure time
del_ttr = 30;   % size of grid - transfer time

targeti = i; targetj = j; t_max = tmax;

%% Time grid formulation
td_grid = 0:del_td:t_max;
ttr_grid = ttr_lb:del_ttr:ttr_ub;
if(td_grid(length(td_grid))<t_max),    td_grid = [td_grid,t_max];    end
if(ttr_grid(length(ttr_grid))<ttr_ub), ttr_grid = [ttr_grid,ttr_ub]; end

x_grid = []; xe = [];
for ind1=1:length(td_grid)-1
    for ind2=1:length(ttr_grid)-1
        if(td_grid(ind1)+ttr_grid(ind2)<t_max)
            x_grid = [x_grid; td_grid(ind1), td_grid(ind1+1), ttr_grid(ind2), ttr_grid(ind2+1)];
        end
    end
end

%% Find optimal (td,ttr) for each grid
if(~isempty(x_grid))
    options = optimoptions('fmincon','Display','off');
    for ind1=1:length(x_grid(:,1))
        grid_lb = [x_grid(ind1,1),x_grid(ind1,3)];
        grid_ub = [x_grid(ind1,2),x_grid(ind1,4)];
        x_tmp1 = fmincon(@myfun,0.5*(grid_lb+grid_ub),ones(1,2),t_max,[],[],grid_lb,grid_ub,@mycon,options);
        [dv,r_dr,i_dr] = delV(targeti,targetj,x_tmp1(1),sum(x_tmp1));
        xe = [xe; x_tmp1,dv,r_dr,i_dr];
    end
end
if(isempty(xe))
    x_tmp1 = [t_max-ttr_ub,ttr_ub];
    [dv,r_dr,i_dr] = delV(targeti,targetj,x_tmp1(1),sum(x_tmp1));
    xe = [x_tmp1,dv,r_dr,i_dr];
    x_tmp2 = [0,ttr_ub];
    [dv,r_dr,i_dr] = delV(targeti,targetj,x_tmp2(1),sum(x_tmp2));
    xe = [xe; x_tmp2,dv,r_dr,i_dr];
end

[a,ai] = sort(xe(:,3),1,'ascend');
x_tmp1 = xe(ai,:);
l = length(x_tmp1(:,1));
check = ones(l,1);
for k=2:l
    for kk=1:k-1
        if((x_tmp1(k,1)<=x_tmp1(kk,1) && x_tmp1(k,1)+x_tmp1(k,2)>=x_tmp1(kk,1)+x_tmp1(kk,2)) || x_tmp1(k,3)>1)
            check(k) = 0;
            break;
        end
    end
end
ispareto = bitor(eq(check,1),eq(check,1));
xe = x_tmp1(ispareto,:);
% save('xe.mat','xe');
end

function J = myfun(x) % x=[td;ttr]
global targeti targetj
J = delV(targeti,targetj,x(1),sum(x));
J = J+0.01*x(2);
end

function [c,ceq] = mycon(x)
global t_max
c(1)=sum(x)-t_max;
ceq=[];
end