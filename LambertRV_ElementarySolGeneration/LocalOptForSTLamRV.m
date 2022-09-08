function [nSol, Sol] = LocalOptForSTLamRV(startOrbit,targetOrbit, tmax , ttfUb , ttfLb)
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
    [num, solution] = Grid2ElemSol(startOrbit,targetOrbit,[tdmin(idxTd),tdmax(idxTd),ttfmin(idxTtf),ttfmax(idxTtf)]);
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

