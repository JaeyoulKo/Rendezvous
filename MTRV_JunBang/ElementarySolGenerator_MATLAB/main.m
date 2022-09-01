clear all; close all; clc;
global x_min x_mat col

n = 123;            % number of debris candidates
K = 7;              % number of spacecraft
dv_max = 3.454;     % maximum delV

ngen_max = 200;     % max number of generation
ncol_max = 100;     % max number of column for each genaration
neval_max = 2000;  % max number of tsp evaluation for each genration

load 'x_mat.mat'    % load DB of elementary solutions
load 'x_min.mat'    % load DB of min-cost elementary solutions
p_mat = ones(n,1);  % profit value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Column Generation                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ngen = 0; neval_total = 0;
col_eval = {}; debset_eval = {};

%% Initial column
col = cell(n,1);
for i=1:n
    col(i) = {i};
end
fprintf('No.initial column: %d\n',length(col));

%% Start column generation procedure
while(1)
    %% Solve dual problem
    l = length(col);
    p_obj = zeros(l,1); A1 = zeros(n,l); A2 = ones(1,l);
    for r=1:l
        s_tmp = col{r};
        p_obj(r) = sum(p_mat(s_tmp));
        A1(s_tmp,r) = 1;
    end
    AA1 = [A1',A2'];
    sense_mat = char(62*ones(l,1));    % 62: '>'
    try
        clear model;
        model.A = sparse(AA1);
        model.obj = [ones(n,1);K];
        model.rhs = p_obj;
        model.sense = sense_mat;
        model.vtype = 'C';
        model.modelsense = 'min';
        model.lb = zeros(n+1,1);
        
        clear params;
        params.outputflag = 0;
        
        result_dual = gurobi(model, params);
        
    catch gurobiError
        fprintf('Error reported\n');
    end
    if(ngen>0)
        fprintf('No.generated column #%d: %d / Obj = %2.4f -> %2.4f\n',ngen,length(col_gen),obj_old,result_dual.objval);
        if(isempty(col_gen) || abs(result_dual.objval-sum(p_mat))<1e-4), break; end  % termination condition: no more column has been generated
    end
    obj_old = result_dual.objval;
    if(ngen>=ngen_max), break; end

    %% Generate columns
    ngen = ngen+1;  neval = 0;
    gencheck = 0;   % 0: no more column
                    % 1: exceed ncol_max
                    % 2: exceed eval_max
    col_gen = {};
    p_new = p_mat-result_dual.x(1:n);
    [aa,ai] = sort(p_new,1,'descend');
    debset = ai(aa>1e-6)';  % subset of debris with positive adjusted reward (descending order)
    if(length(debset)>15), debset = debset(1:15); end
    if(~isincluded(sort(debset),debset_eval))
        debset_eval = [debset_eval;sort(debset)];
        col_tmp_i = 1;
    else
        col_tmp_i = [];
    end

    while(~isempty(col_tmp_i))  % explore columns in lexicographic order
        col_tmp = sort(debset(col_tmp_i));
        [iseval,ind] = isincluded(col_tmp,col_eval);
        if(isincluded(col_tmp,col))
            col_tmp_i = getnextcol(col_tmp_i,length(debset),0);
        elseif(iseval)
            if(col_eval{ind,2}==1)
                col_tmp_i = getnextcol(col_tmp_i,length(debset),0);
            else
                col_tmp_i = getnextcol(col_tmp_i,length(debset),1);
            end
        else
            isfeasible = fcheck(col_tmp,dv_max);
            if(isfeasible==1)
                if(sum(p_new(col_tmp))>result_dual.x(n+1))
                    col_gen = [col_gen;col_tmp];
                end
                col_eval = [col_eval;{},col_tmp,1];
                col_tmp_i = getnextcol(col_tmp_i,length(debset),0);
            else
                col_eval = [col_eval;{},col_tmp,0];
                col_tmp_i = getnextcol(col_tmp_i,length(debset),1);
            end
            neval_total = neval_total+1;
        end
        neval = neval+1;

        if(length(col_gen)>=ncol_max)
            gencheck = 1;
            break;
        elseif(neval>neval_max)
            gencheck = 2;
            break;
        end
    end        

    if(isempty(col_gen))
        [col_pr,status] = pricing(p_new,dv_max,result_dual.x(n+1),ncol_max,2*ncol_max);
        if(~isempty(col_pr))
            col_gen = col_pr;
            for i=1:length(col_pr)
                col_eval = [col_eval;{},col_pr{i},1];
            end
        end
    end
    col = [col;col_gen];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Solve MPLR with current column set                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l = length(col);
p_obj = zeros(l,1); A1 = zeros(n,l); A2 = ones(1,l);
for r=1:l
    s_tmp = col{r};
    p_obj(r) = sum(p_mat(s_tmp));
    A1(s_tmp,r) = 1;
end
sense_mat = char(60*ones(n+1,1));   % 60: '<'
try
    clear model;
    model.A = sparse([A1;A2]);
    model.obj = p_obj;
    model.rhs = [ones(n,1);K];
    model.sense = sense_mat;
    model.vtype = 'C';
    model.modelsense = 'max';

    clear params;
    params.outputflag = 0;

    result_primal = gurobi(model, params);
catch gurobiError
    fprintf('Error reported\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Solve MP                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l = length(col);
p_obj = zeros(l,1); A1 = zeros(n,l); A2 = ones(1,l);
for r=1:l
    s_tmp = col{r};
    p_obj(r) = sum(p_mat(s_tmp));
    A1(s_tmp,r) = 1;
end
sense_mat = char(60*ones(n+1,1));   % 60: '<'
try
    clear model;
    model.A = sparse([A1;A2]);
    model.obj = p_obj;
    model.rhs = [ones(n,1);K];
    model.sense = sense_mat;
    model.vtype = 'B';
    model.modelsense = 'max';

    clear params;
    params.outputflag = 0;

    result_final = gurobi(model, params);

catch gurobiError
    fprintf('Error reported\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Post-process                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if((result_primal.objval-result_final.objval)/result_final.objval>0.05)
    while(1)
    final_route = {}; aa = [];
    for i=1:l
        if(result_final.x(i)>1e-2)
            final_route = [final_route; {},col{i}];
            aa = [aa; length(col{i})];
        end
    end
    [aa,ai] = sort(aa,1,'ascend');
    final_route = final_route(ai,:);

    deb_visit = []; deb_nonvisit = []; col_gen = {};
    for i=1:length(final_route(:,1))
        deb_visit = [deb_visit,final_route{i}];
    end
    for i=1:n
        if(~ismember(i,deb_visit)), deb_nonvisit = [deb_nonvisit,i]; end
    end
    fprintf('Current obj = %2.4f / No.nonvisited debris: %d\n',result_final.objval,length(deb_nonvisit));
    if(isempty(deb_nonvisit) || abs(result_final.objval-obj_old)<1e-4), break; end

    for i=1:length(final_route(:,1))
        if(length(final_route{i})<15)
            p_new = zeros(n,1);
            p_new([final_route{i},deb_nonvisit]) = 1;
            [col_pr,status] = pricing(p_new,dv_max,length(final_route{i}),5,5);
            if(~isempty(col_pr))
                col_gen = [col_gen;col_pr];
                for j=1:length(col_pr(:,1))
                    col_eval = [col_eval;{},col_pr(j),1];
                end
                break;
            end
        end
    end
    if(isempty(col_gen))
        break;
    else
        col = [col;col_gen];
        obj_old = result_final.objval;
    end
    
    l = length(col);
    p_obj = zeros(l,1); A1 = zeros(n,l); A2 = ones(1,l);
    for r=1:l
        s_tmp = col{r};
        p_obj(r) = sum(p_mat(s_tmp));
        A1(s_tmp,r) = 1;
    end
    sense_mat = char(60*ones(n+1,1));   % 60: '<'
    try
        clear model;
        model.A = sparse([A1;A2]);
        model.obj = p_obj;
        model.rhs = [ones(n,1);K];
        model.sense = sense_mat;
        model.vtype = 'B';
        model.modelsense = 'max';

        clear params;
        params.outputflag = 0;

        result_final = gurobi(model, params);

    catch gurobiError
        fprintf('Error reported\n');
    end
    end
end

final_route = {};
for i=1:l
    if(result_final.x(i)>1e-2)
        [s_tmp,c_tmp] = tsp(col{i});
        final_route = [final_route; {},s_tmp,sum(p_mat(s_tmp)),c_tmp];
    end
end
fprintf('\n');
for i=1:length(final_route(:,1))
    fprintf('Route %d: %s\n',i,mat2str(final_route{i,1}));
end
fprintf('\nFinal obj: %f\n',result_final.objval);
fprintf('Optimality gap: %f\n',(result_primal.objval-result_final.objval)/result_final.objval*100);
