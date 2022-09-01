function [X,Xdot] = getStateAtT(a,e,i,o,w,nu0,t)
nu = TrueAnomalyAtdeltaT(a,e,i,o,w,nu0,t);
[X,Xdot]=OrbitElem2RV(a,e,i,o,w,nu);
end

