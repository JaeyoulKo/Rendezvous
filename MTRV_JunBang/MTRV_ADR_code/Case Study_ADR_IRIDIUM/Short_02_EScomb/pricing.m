function [col_pr,status] = pricing(p_mat,dv_max,p_cutoff,ncol_max,nsol_max)
global x_mat x_min col

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Solve Pricing Problem                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Identify columns that can improve the current solution
% 
% col_pr: set of columns (profit > p_cutoff & cost < dv_max)
% status = 0: no additional column 
%          1: find optimal solution of the pricing problem
%          2: exceed nsol_max
%          3: exceed time limit
% p_mat: profit matrix
% dv_max: maximum delV
% p_cutoff: lower bound for the profit
% ncol_max: maximum number of columns
% nsol_max: maximum number of feasible solutions explored

tser = 7;  % service time

col_pr = {}; status = 0;    
n = length(p_mat);
Q = [0:n]; p_mat = [0;p_mat];
x_tsp = []; n_arc = zeros(n+1,n+1); p_obj = [];
for i=1:n+1
    for j=1:n+1
        if(i~=j)
            ind = Q(i)*1000+Q(j);
            x_tmp = x_mat(x_min(ind,length(x_min(1,:))-1):x_min(ind,length(x_min(1,:))),2:4);
            n_arc(i,j) = length(x_tmp(:,1));
            x_tsp = [x_tsp;(i*1000+j-1001)*ones(n_arc(i,j),1),x_tmp];
            p_obj = [p_obj;p_mat(j)*ones(n_arc(i,j),1)];
        end
    end
end
l = length(x_tsp(:,1)); M = 1e6;

A1 = zeros(n,l); A2 = zeros(n+1,l); A3 = zeros(1,l); A4 = zeros(1,l); A5 = zeros(n,l);
ind = 0;
for i=0:n
    for j=0:n
        if(i~=j)
            for p=1:n_arc(i+1,j+1)
                ind = ind+1;
                if(i~=0),   A1(i,ind) = 1; end
                A2(i+1,ind) = -1;
                if(i==0),   A3(1,ind) = 1; end
                A4(1,ind) = x_tsp(ind,4);
                if(i~=0),   A5(i,ind) = M-x_tsp(ind,2); end
            end
            ind2 = getind(j,i,1,n_arc)-1;
            for p=1:n_arc(j+1,i+1)
                ind2 = ind2+1;
                A2(i+1,ind2) = 1;
                if(i~=0),   A5(i,ind2) = x_tsp(ind2,2)+x_tsp(ind2,3); end
            end
        end
    end
end
sense_mat = char([60*ones(n,1);61*ones(n+2,1);60*ones(n+1,1)]); % 61:'=', 60:'<'

try
    clear model;
    model.A = sparse([A1;A2;A3;A4;A5]);
    model.obj = p_obj;
    model.rhs = [ones(n,1);zeros(n+1,1);1;dv_max;(M-tser)*ones(n,1)];
    model.sense = sense_mat;
    model.vtype = 'B';
    model.modelsense = 'max';

    clear params;
    params.outputflag = 0;
    params.cutoff = p_cutoff+1e-4;
    params.improvestartnodes = 10;
    params.improvestarttime = 10;
    params.poolsearchmode = 2;
    params.poolsolutions = nsol_max;
    params.solutionlimit = nsol_max;
    params.timelimit = 300;

    result = gurobi(model, params);

catch gurobiError
    fprintf('Error reported\n');
end

if(length(result.status)==7 && sum(result.status=='OPTIMAL')==7)
    status = 1;
elseif(length(result.status)==14 && sum(result.status=='SOLUTION_LIMIT')==14)
    status = 2;
elseif(length(result.status)==10 && sum(result.status=='TIME_LIMIT')==10)
    status = 3;
end
if(status==0 || ~isfield(result,'pool')), return; end

for j=1:length(result.pool)
    s_tmp = [];
    for i=1:l
        if(abs(result.pool(j).xn(i)-1)<1e-4)
    %         disp([x_tsp(i,1),x_tsp(i,2:4)])
            ind = x_tsp(i,1);
            ii = floor(ind/1000);
            jj = ind-1000*ii;
            s_tmp = [s_tmp; ii+1,jj+1];
        end
    end
    s_tmp2 = s_tmp(1,2);
    for i=1:length(s_tmp(:,1))-2
        ind = 0;
        for j=1:length(s_tmp(:,1))
            if(s_tmp(j,1)==s_tmp2(i))
                ind = j;
                break;
            end
        end
        s_tmp2 = [s_tmp2,s_tmp(ind,2)];
    end
    s = Q(1,s_tmp2);
    
    col_tmp = sort(s);
    if(~isincluded(col_tmp,col_pr) && ~isincluded(col_tmp,col))
        col_pr = [col_pr;col_tmp];
        if(length(col_pr)>=ncol_max)
            break;
        end
    end
end
end

function ind = getind(i,j,p,n_arc)
    n = length(n_arc(:,1))-1;
    if(i>n || j>n || p>n_arc(i+1,j+1))
        disp('input error')
    else
        ind = 0;
        if(i>0)
            ind = sum(sum(n_arc(1:i,:)));
        end
        if(j>0)
            ind = ind+sum(n_arc(i+1,1:j));
        end
        if(p>0)
            ind = ind+p;
        end
    end
end