function r = Equinoctial2Distance(p,f,g,h,k,L)
%EQUINOCTIAL Orbit Element 를 입력받아 Central Body와의 거리 r 을 return
e = sqrt(f^2+g^2);
o = atan2(k,h);
w = atan2(g,f)-o; %?
nu = L-o-w;
r = p/(1+e*cos(nu));
end

