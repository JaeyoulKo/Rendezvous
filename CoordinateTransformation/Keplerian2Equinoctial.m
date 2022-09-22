function [p,f,g,h,k,L] = Keplerian2Equinoctial(a,e,i,o,w,nu)
% Keplerian Orbit Element 를 Equinoctial Orbit Element로 변환
p = a*(1-e.^2);
f = e*cos(w+o);
g = e*sin(w+o);
h = tan(i/2)*cos(o);
k = tan(i/2)*sin(o);
L = o+w+nu;
end

