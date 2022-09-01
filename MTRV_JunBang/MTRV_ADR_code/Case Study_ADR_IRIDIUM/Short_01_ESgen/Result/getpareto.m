clear all; close all; clc;

n = 100;
tmax = 7*24*3600;
Jmax = 1;

x_mat = load('x_local0.txt');
x_min = zeros(n*1000+n,6);
for i0=1:n
    x_min(i0,:) = [i0,0,0,0,i0,i0];
end

for i0=1:n
    x_mat = [x_mat;i0*1000,tmax,0,0];
    i_s = length(x_mat(:,1));
    x_min(i0*1000,:) = [i0*1000,tmax,0,0,i_s,i_s];

    fname = strcat("x_local",int2str(i0),".txt");
    x_local = load (fname);
    ind = x_local(1,1); ind_check = 1; ind_start = 1; i_s = i_s+1;
    while(ind_check)
        disp(ind)
        x_tmp = []; ind_check = 0;
        for i=ind_start:length(x_local(:,1))
            if(x_local(i,1)==ind)
                x_tmp = [x_tmp; x_local(i,:)];
            elseif(x_local(i,1)>ind)
                ind_check = 1;
                ind_end = i-1;
                break;
            end
        end
        if(~isempty(x_tmp))
            [a,ai] = sort(x_tmp(:,4),1,'ascend');
            x_tmp = x_tmp(ai,:);

            l = length(x_tmp(:,1));
            check = ones(l,1);
            for k=2:l
                for j=1:k-1
                    if((x_tmp(k,2)<x_tmp(j,2) && x_tmp(k,2)+x_tmp(k,3)>x_tmp(j,2)+x_tmp(j,3)) || x_tmp(k,4)>Jmax || (abs(x_tmp(k,2)-x_tmp(j,2))<3*3600 && abs(x_tmp(k,3)-x_tmp(j,3))<3600))
                        check(k) = 0;
                        break;
                    end
                end
            end
            ispareto = bitor(eq(check,1),eq(check,1));
            x_pareto = x_tmp(ispareto,:);

            x_mat = [x_mat; x_pareto];
            i_e = length(x_mat(:,1));
            x_min(x_pareto(1,1),:) = [x_pareto(1,:),i_s, i_e];
            i_s = i_e+1;

            ind = x_local(i,1);
            ind_start = ind_end+1;
        end
    end
end

save('x_mat.mat','x_mat');
save('x_min.mat','x_min');
