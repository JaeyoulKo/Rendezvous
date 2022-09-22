function [nSol, Sol] = GetElementarySolution(startOrbit,targetOrbit, tmax , ttfUb , ttfLb)
% Local Optima (DelV Min) For Single Target Lambert Rendezvous
% startOrbit : 6 orbit elements [km & rad]
% targetOrbit : 6 orbit elements [km & rad]
% td : departing time [sec]
% ttf : transfer time [sec]
% Sol : [t_d;t_tr; J(=delV)]
nSol=0;
Sol = zeros(3,3000);
mu = 398600.442; %km^3/sec^2

% Step 1. Partition the whol solution space into subspaces,
% size shoud be smaller than (del_td,A by del_ttr,A)
% basic setting : (0.9deltd_a by 0.9delttr_a)

n1 = sqrt(mu/startOrbit(1)^3);
n2 = sqrt(mu/targetOrbit(1)^3);
deltdA = pi/(n1);
delttfA = 2*pi/n2;
uGridTd = deltdA*0.9;
uGridTtf = delttfA*0.9;
nGridTd = ceil(tmax / uGridTd);
nGridTtf = ceil((ttfUb-ttfLb)/uGridTtf);
tdmin = (0:nGridTd-1)*uGridTd;
tdmax = (1:nGridTd)*uGridTd;
tdmax(end) = tmax;
ttfmin = (0:nGridTtf-1)*uGridTtf+ttfLb;
ttfmax = (1:nGridTtf)*uGridTtf+ttfLb;
ttfmax(end) = ttfUb;

%%Step 2. For each grid implement the gradient-based algorithm


% for idx = 0:nGridTd*nGridTtf-1
%     idxTd = fix(idx/nGridTtf)+1;
%     idxTtf = rem(idx,nGridTtf)+1;
%     disp(idxTd);
%     disp(idxTtf);
% end
idx = 1;
idxTd = 1;
idxTtf = 1;
while (idx<=nGridTd*nGridTtf)
    [num, solution] = GetElementarySolutionsInGrid(startOrbit,targetOrbit,[tdmin(idxTd),tdmax(idxTd),ttfmin(idxTtf),ttfmax(idxTtf)]);
    Sol(:,nSol+1:nSol+num) = solution;
    nSol=nSol+num;
    idx=idx+1;
    if (idxTtf == nGridTtf)
        idxTtf=1;
        idxTd=idxTd+1;
    else
        idxTtf=idxTtf+1;
    end
end
Sol = Sol(:,1:nSol);
end

function [nSol, Sol] = GetElementarySolutionsInGrid(startOrbit,targetOrbit,grid)

nSol=0;
Sol = zeros(3,20);
minGridSize = 1000;
x0eps = 0.001;

% grid = [td_min, td_max, ttf_min, ttf_max]
conA = [1 0 ; -1 0 ; 0 1 ; 0 -1];
conb = [grid(2) ; -grid(1) ; grid(4) ; -grid(3)];
% options = optimoptions(@fmincon,'Algorithm','sqp','Display','off');
options = optimoptions(@fmincon,'Algorithm','sqp','Display','notify-detailed','MaxFunctionEvaluations',5000,'MaxIteration',5000,'StepTolerance',10^-15);

[xLL , valLL] = fmincon(@(t)CalcDelVForLamRVFromOrbitsTdTtf(startOrbit,targetOrbit,t) ...
    ,[grid(1)+x0eps grid(3)+x0eps],conA,conb,[],[],[],[],[],options);
[xLR , valLR] = fmincon(@(t)CalcDelVForLamRVFromOrbitsTdTtf(startOrbit,targetOrbit,t) ...
    ,[grid(2)-x0eps grid(3)+x0eps],conA,conb,[],[],[],[],[],options);
[xUL , valUL] = fmincon(@(t)CalcDelVForLamRVFromOrbitsTdTtf(startOrbit,targetOrbit,t) ...
    ,[grid(1)+x0eps grid(4)-x0eps],conA,conb,[],[],[],[],[],options);
[xUR , valUR] = fmincon(@(t)CalcDelVForLamRVFromOrbitsTdTtf(startOrbit,targetOrbit,t) ...
    ,[grid(2)-x0eps grid(4)-x0eps],conA,conb,[],[],[],[],[],options);

x=[xLL;xLR;xUL;xUR];
J=[valLL;valLR;valUL;valUR];
%%Step 3. If all the solutions are located on the boundary stop exploring
if isAllBoundary(x,grid)
    Sol = nan;
    return;
elseif isAllConvergeSinglePoint(x,J)    %step 4
    nSol = 1;
    [val, idx] = min(J);
    Sol = [x(idx,:)' ; val];
elseif (grid(2)-grid(1)<minGridSize && grid(4)-grid(3)<minGridSize)
    Sol = nan;
    return;
else        %step 5
    newGridLL = [grid(1), (grid(1)+grid(2))/2 , grid(3), (grid(3)+grid(4))/2];
    newGridLR = [(grid(1)+grid(2))/2 , grid(2), grid(3), (grid(3)+grid(4))/2];
    newGridUL = [grid(1), (grid(1)+grid(2))/2 , (grid(3)+grid(4))/2 , grid(4)];
    newGridUR = [(grid(1)+grid(2))/2 , grid(2), (grid(3)+grid(4))/2 , grid(4)];
    [nSolLL, solLL] = GetElementarySolutionsInGrid(startOrbit,targetOrbit,newGridLL);
    [nSolLR, solLR] = GetElementarySolutionsInGrid(startOrbit,targetOrbit,newGridLR);
    [nSolUL, solUL] = GetElementarySolutionsInGrid(startOrbit,targetOrbit,newGridUL);
    [nSolUR, solUR] = GetElementarySolutionsInGrid(startOrbit,targetOrbit,newGridUR);
    if (nSolLL~=0)
        Sol(:,nSol+1:nSol+nSolLL) = solLL;
        nSol=nSol+nSolLL;
    end
    if (nSolLR~=0)
        Sol(:,nSol+1:nSol+nSolLR) = solLR;
        nSol=nSol+nSolLR;
    end
    if (nSolUL~=0)
        Sol(:,nSol+1:nSol+nSolUL) = solUL;
        nSol=nSol+nSolUL;
    end
    if (nSolUR~=0)
        Sol(:,nSol+1:nSol+nSolUR) = solUR;
        nSol=nSol+nSolUR;
    end
    Sol = Sol(:,1:nSol);

end
end

    function flag = isBoundary(X,grid)
        boundaryTol = 0.01;
        flag = (abs(X(1,1)-grid(1))<boundaryTol) || (abs(X(1,1)-grid(2))<boundaryTol) || (abs(X(1,2)-grid(3))<boundaryTol) || (abs(X(1,2)-grid(4))<boundaryTol);
    end
    function flag = isAllBoundary(X,grid)
        flag = isBoundary(X(1,:),grid)&&isBoundary(X(2,:),grid)&&isBoundary(X(3,:),grid)&&isBoundary(X(4,:),grid);
    end

    function flag = isAllConvergeSinglePoint(X,J)
        sameTimeTol = 200;
        sameJTol = 0.01;
        minX = min(X);
        maxX = max(X);
        flag = (max(J)-min(J)<sameJTol) && (maxX(1)-minX(1) < sameTimeTol) && (maxX(2)-minX(2) < sameTimeTol);
    end

    function flag = isAddSol(X,grid) % C code와 유사하게 구현. 필요시, is AllConvergeSinglePoint를 이 함수로 대치해 사용.
        minGridSize=1800;
        flagBnd = (~isBoundary(X(1,:),grid))&&(~isBoundary(X(2,:),grid))&&(~isBoundary(X(3,:),grid))&&(~isBoundary(X(4,:),grid));
        flagGrid = (grid(2)-grid(1)<minGridSize || grid(4)-grid(3)<minGridSize);
        flag = flagBnd||flagGrid;
    end

