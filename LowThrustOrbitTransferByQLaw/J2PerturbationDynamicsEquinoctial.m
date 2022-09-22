function fPerturb = J2PerturbationDynamicsEquinoctial(state,parameter)
mu = parameter(1);
Re = parameter(2);
J2 = parameter(3);

p = state(1);
f = state(2);
g = state(3);
h = state(4);
k = state(5);
L = state(6);

% addpath("./../CoordinateTransformation/")
r = Equinoctial2Distance(p,f,g,h,k,L);
fr = -3*mu*J2*Re^2/(2*r^4) * ( 1 - 12*(h*sin(L) - k*cos(L))^2 / (1+h^2+k^2)^2 );
ft = -12*mu*J2*Re^2/(r^4) * ( (h*sin(L) - k*cos(L))*(h*cos(L)+k*sin(L)) / (1+h^2+k^2)^2 );
fh = -6*mu*J2*Re^2/(r^4) * ( (1-h^2-k^2)*(h*sin(L) - k*cos(L)) / (1+h^2+k^2)^2 );

fPerturb = [ft;fr;fh];
end

