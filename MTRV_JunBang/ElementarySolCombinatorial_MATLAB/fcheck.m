function isfeasible = fcheck(col,dv_max)
global x_mat x_min

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Check the Feasibility of the Column                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% isfeasible = 0: min-cost to visit all debris in a subset > max delV
%              1: min-cost to visit all debris in a subset < max delV
% col: column (debris sequence)
% dv_max: maximum delV

tser = 5;    % service time

n = length(col);  
col = [0,col];

%% Solve TSP with time constraints
x_tsp = []; n_arc = zeros(n+1,n+1);
for i=1:n+1
    for j=1:n+1
        if(i~=j)
            ind = col(i)*1000+col(j);
            x_tmp = x_mat(x_min(ind,length(x_min(1,:))-1):x_min(ind,length(x_min(1,:))),2:4);
            n_arc(i,j) = length(x_tmp(:,1));
            x_tsp = [x_tsp;(i*1000+j-1001)*ones(length(x_tmp(:,1)),1),x_tmp];
        end
    end
end
l = length(x_tsp(:,1)); c_obj = x_tsp(:,4);

A1 = zeros(n+1,l); A2 = zeros(n+1,l); A3 = zeros(n,l);
ind = 0;
for i=0:n
    for j=0:n
        if(i~=j)
            for p=1:n_arc(i+1,j+1)
                ind = ind+1;
                A1(j+1,ind) = 1;
                A2(i+1,ind) = 1;
                if(i~=0),   A3(i,ind) = -x_tsp(ind,2); end
            end
            ind2 = getind(j,i,1,n_arc)-1;
            for p=1:n_arc(j+1,i+1)
                ind2 = ind2+1;
                if(i~=0),   A3(i,ind2) = x_tsp(ind2,2)+x_tsp(ind2,3); end
            end
        end
    end
end
sense_mat = char([61*ones(2*n+2,1); 60*ones(n,1)]);   % 61:'=', 60:'<'

try
    clear model;
    model.A = sparse([A1;A2;A3]);
    model.obj = c_obj;
    model.rhs = [ones(2*n+2,1);-1*tser*ones(n,1)];
    model.sense = sense_mat;
    model.vtype = 'B';
    model.modelsense = 'min';

    clear params;
    params.outputflag = 0;
    params.cutoff = dv_max;
    params.bestobjstop = dv_max;
    params.improvestartnodes = 10;
    params.improvestarttime = 1;

    result = gurobi(model, params);

catch gurobiError
    fprintf('Error reported\n');
end

if(isfield(result,'objval') && result.objval<dv_max)
    isfeasible = 1;
else
    isfeasible = 0;
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