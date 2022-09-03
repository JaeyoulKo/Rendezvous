function [X,Xdot] = getStateAtT(a,e,i,o,w,nu0,t)
%Aerospace Toolbox가 있다면, keplerian2ijk 를 사용. (단위 다름)
nu = TrueAnomalyAtdeltaT(a,e,i,o,w,nu0,t);
[X,Xdot]=OrbitElem2RV(a,e,i,o,w,nu);
end

