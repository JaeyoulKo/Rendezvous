function delV = DelV_LambRend_tdttf(startOrbit,targetOrbit, t)
    Mu = 398600.442; %km^3/sec^2
    td=t(1);
    ttf=t(2);
    % addpath('./../Lambert_Gooding/')
    [Xs,Xsdot] = getStateAtT(startOrbit(1),startOrbit(2),startOrbit(3),startOrbit(4),startOrbit(5),startOrbit(6),td);
    [Xt,Xtdot] = getStateAtT(targetOrbit(1),targetOrbit(2),targetOrbit(3),targetOrbit(4),targetOrbit(5),targetOrbit(6),td+ttf);
    delV = DelV_LambRend(Mu, Xs, Xsdot, Xt, Xtdot, ttf);
end