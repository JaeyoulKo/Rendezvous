function delV = CalcDelVForLamRVFromOrbitsTdTtf(startOrbit,targetOrbit, t)
    Mu = 398600.442; %km^3/sec^2
    td=t(1);
    ttf=t(2);
    % addpath('./../Lambert_Gooding/')
    [Xs,Xsdot] = Keplerian2Geocentric(startOrbit(1),startOrbit(2),...
        startOrbit(3),startOrbit(4),startOrbit(5),...
        NuAtT(startOrbit(1),startOrbit(2),startOrbit(6),td));
    [Xt,Xtdot] = Keplerian2Geocentric(targetOrbit(1),targetOrbit(2),...
        targetOrbit(3),targetOrbit(4),targetOrbit(5),...
        NuAtT(targetOrbit(1),targetOrbit(2),targetOrbit(6),td+ttf));
    delV = CalcDelVForLamRV(Mu, Xs, Xsdot, Xt, Xtdot, ttf);
end

function delV = CalcDelVForLamRV(Mu, r1vec, v1vec, r2vec, v2vec, delT)

[~, delV1Vec, delV2Vec] = CalcDelVVectorForLamRV(Mu, r1vec, v1vec, r2vec, v2vec, delT);
delV = min(vecnorm(delV1Vec)+vecnorm(delV2Vec));

end

function [nSol, delV1Vec, delV2Vec] = CalcDelVVectorForLamRV(Mu, r1vec, v1vec, r2vec, v2vec, delT)
% addpath("./../Lambert_Gooding/");
[nSol, v1Lam, v2Lam] = Lambert(Mu, r1vec, r2vec, delT);

delV1Vec = v1Lam-v1vec;
delV2Vec = -v2Lam+v2vec;

end

