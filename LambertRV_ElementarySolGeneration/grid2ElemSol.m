function [nSol, Sol] = grid2ElemSol(startOrbit,targetOrbit,grid)

nSol=0;
Sol = zeros(3,20);

x0eps = 0.001;

% grid = [td_min, td_max, ttf_min, ttf_max]
conA = [1 0 ; -1 0 ; 0 1 ; 0 -1];
conb = [grid(2) ; -grid(1) ; grid(4) ; -grid(3)];
% options = optimoptions(@fmincon,'Algorithm','sqp','Display','off');
options = optimoptions(@fmincon,'Display','off');
[xLL , valLL] = fmincon(@(t)DelV_LambRend_tdttf(startOrbit,targetOrbit,t) ...
    ,[grid(1)+x0eps grid(3)+x0eps],conA,conb,[],[],[],[],[],options);
[xLR , valLR] = fmincon(@(t)DelV_LambRend_tdttf(startOrbit,targetOrbit,t) ...
    ,[grid(2)-x0eps grid(3)+x0eps],conA,conb,[],[],[],[],[],options);
[xUL , valUL] = fmincon(@(t)DelV_LambRend_tdttf(startOrbit,targetOrbit,t) ...
    ,[grid(1)+x0eps grid(4)-x0eps],conA,conb,[],[],[],[],[],options);
[xUR , valUR] = fmincon(@(t)DelV_LambRend_tdttf(startOrbit,targetOrbit,t) ...
    ,[grid(2)-x0eps grid(4)-x0eps],conA,conb,[],[],[],[],[],options);

x=[xLL;xLR;xUL;xUR];
%%Step 3. If all the solutions are located on the boundary stop exploring
if isAllBoundary(x,grid)
    Sol = nan;
    return;
elseif isAddSol(x,grid)    %step 4
    nSol = 1;
    J=[valLL;valLR;valUL;valUR];
    [val, idx] = min(J);
    Sol = [x(idx,:)' ; val];
else        %step 5
    newGridLL = [grid(1), (grid(1)+grid(2))/2 , grid(3), (grid(3)+grid(4))/2];
    newGridLR = [(grid(1)+grid(2))/2 , grid(2), grid(3), (grid(3)+grid(4))/2];
    newGridUL = [grid(1), (grid(1)+grid(2))/2 , (grid(3)+grid(4))/2 , grid(4)];
    newGridUR = [(grid(1)+grid(2))/2 , grid(2), (grid(3)+grid(4))/2 , grid(4)];
    [nSolLL, solLL] = grid2ElemSol(startOrbit,targetOrbit,newGridLL);
    [nSolLR, solLR] = grid2ElemSol(startOrbit,targetOrbit,newGridLR);
    [nSolUL, solUL] = grid2ElemSol(startOrbit,targetOrbit,newGridUL);
    [nSolUR, solUR] = grid2ElemSol(startOrbit,targetOrbit,newGridUR);
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

    function delV = DelV_LambRend_tdttf(startOrbit,targetOrbit, t)
        Mu = 398600.442; %km^3/sec^2
        td=t(1);
        ttf=t(2);
        % addpath('./../Lambert_Gooding/')
        [Xs,Xsdot] = getStateAtT(startOrbit(1),startOrbit(2),startOrbit(3),startOrbit(4),startOrbit(5),startOrbit(6),td);
        [Xt,Xtdot] = getStateAtT(targetOrbit(1),targetOrbit(2),targetOrbit(3),targetOrbit(4),targetOrbit(5),targetOrbit(6),td+ttf);
        delV = DelV_LambRend(Mu, Xs, Xsdot, Xt, Xtdot, ttf);
    end

    function flag = isBoundary(X,grid)
        boundaryTol = 0.01;
        flag = (abs(X(1,1)-grid(1))<boundaryTol) || (abs(X(1,1)-grid(2))<boundaryTol) || (abs(X(1,2)-grid(3))<boundaryTol) || (abs(X(1,2)-grid(4))<boundaryTol);
    end
    function flag = isAllBoundary(X,grid)
        flag = isBoundary(X(1,:),grid)&&isBoundary(X(2,:),grid)&&isBoundary(X(3,:),grid)&&isBoundary(X(4,:),grid);
    end

    function flag = isAllConvergeSinglePoint(X,J)
        sameTimeTol = 20;
        sameJTol = 0.0001;
        minX = min(X);
        maxX = max(X);
        flag = (max(J)-min(J)<sameJTol) && (maxX(1)-minX(1) < sameTimeTol) && (maxX(2)-minX(2) < sameTimeTol);
    end

    function flag = isAddSol(X,grid)
        minGridSize=1800;
        flagBnd = (~isBoundary(X(1,:),grid))&&(~isBoundary(X(2,:),grid))&&(~isBoundary(X(3,:),grid))&&(~isBoundary(X(4,:),grid));
        flagGrid = (grid(2)-grid(1)<minGridSize || grid(4)-grid(3)<minGridSize);
        flag = flagBnd||flagGrid;
    end