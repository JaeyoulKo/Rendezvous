function stateDot = TwoBodyVOPDynamicsEquinoctial(state, controlInput, parameter)
% TwoBodyVOPDynamicsEquinoctial
mu = parameter(1); Tmax = parameter(2); c = parameter(3);
p = state(1); f = state(2); g = state(3); h = state(4); k = state(5); L = state(6);
% Mass = state(7)
% constrolInput = [f_t ; f_r ; f_n] = Fcos(beta)cos(alpha) ; Fcos(beta)sin(alpha) ; Fsin(beta)

stateDot = zeros(7,1);
A=zeros(6,3);

w   = 1 + f*cos(L) + g*sin(L);
sSqr  = 1 + h^2 + k^2; 
A(1,1) = 2 * p/w * sqrt(p/mu);
A(1,2) = 0;
A(1,3) = 0;
A(2,1) = sqrt(p/mu)/w * ( (w+1)*cos(L) + f );
A(2,2) = sqrt(p/mu) * sin(L);
A(2,3) = - sqrt(p/mu) * g/w * ( h*sin(L) - k*cos(L) );
A(3,1) = sqrt(p/mu)/w * ( (w+1)*sin(L) + g );
A(3,2) = - sqrt(p/mu) * cos(L);
A(3,3) = sqrt(p/mu) * f/w * ( h*sin(L) - k*cos(L) );
A(4,1) = 0;
A(4,2) = 0;
A(4,3) = sqrt(p/mu) * sSqr * cos(L) / (2*w);
A(5,1) = 0;
A(5,2) = 0;
A(5,3) = sqrt(p/mu) * sSqr * sin(L) / (2*w);
A(6,1) = 0;
A(6,2) = 0;
A(6,3) = sqrt(p/mu) * ( h*sin(L) - k*cos(L) ) / w;

%% control modified
stateDot(1:6,1) = A * (controlInput) + [0;0;0;0;0;sqrt(mu*p)*(w/p)^2];
stateDot(7,1)   = - Tmax/c;
stateDot = real(stateDot);
end