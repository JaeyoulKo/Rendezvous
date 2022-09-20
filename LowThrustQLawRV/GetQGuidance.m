function [QDot, ControlInput] = GetQGuidance(chaserState, targetState, F, mu, W)
% Q law Guidance의 output 값인 control Input 값과 QDot 을 계산하는 함수
% TODO : QDot 관련 Unit Test & Debuggin 필요
% WeightParameter = [Wa; Wf; Wg; Wh; Wk; WL] / norm([Wa; Wf; Wg; Wh; Wk; WL]);
% chaser & target State = [p;f;g;h;k;L;Mass(not use)]

% QGrad = ComputeGradientOfQ(chaserState, targetState, F, mu, W);
QGrad = ComputeGradientOfQ(chaserState, targetState, 1, mu, W);

D = CalcD(chaserState,QGrad,mu);
% QDot_TEMP = D;

alpha = atan2(-D(2),-D(1));
beta  = atan( -D(3)/sqrt(D(1)^2+D(2)^2) ); % Direction of Thrust Force
ControlInput = D/norm(D) * F; % Magnitude of Thrust Force
QDot = D(1)*cos(beta)*cos(alpha) + D(2)*cos(beta)*sin(alpha) + D(3)*sin(beta); % TODO: Value Check, Debugging

end

function gradientVector = ComputeGradientOfQ(chaserState, targetState, F, mu, WeightParameter)
% Q값의 Numerical Gradient를 구하는 함수.
% chaser & target State = [p;f;g;h;k;L;Mass(not use)]
% gradientVector = [dQ/da ; dQ/df ; dQ/dg ; ... ]

gradientVector = zeros(5,1);
% gradientVector = zeros(6,1);
dElem = 1e-9;
f = chaserState(2);
g = chaserState(3);
delState = zeros(7,1);
delState(1) = dElem*(1-f^2-g^2);
gradientVector(1) = ( CalcQ(chaserState+delState,targetState,F,mu,WeightParameter) ...
    - CalcQ(chaserState-delState,targetState,F,mu,WeightParameter) ) / (2*dElem);  %need unit test and double check
delState(1)=0;
for orbitIdx=2:6 %6
    delState(orbitIdx)=dElem;
    gradientVector(orbitIdx) = ( CalcQ(chaserState+delState,targetState,F,mu,WeightParameter) ...
    - CalcQ(chaserState-delState,targetState,F,mu,WeightParameter) ) / (2*dElem);
    delState(orbitIdx)=0;
end
end

function D = CalcD(state, QDotPartial, mu)
% TwoBodyVOPDynamicsEquinoctial
p = state(1); f = state(2); g = state(3); h = state(4); k = state(5); L = state(6);
% Mass = state(7)
% constrolInput = [f_t ; f_r ; f_n] = Fcos(beta)cos(alpha) ; Fcos(beta)sin(alpha) ; Fsin(beta)

A=zeros(6,3);

a = p*(1-f^2-g^2);
e = sqrt(f^2+g^2);
nu = L-atan2(g,f);
w   = 1 + f*cos(L) + g*sin(L);
r = p/(1+e*cos(nu));
sSqr  = 1 + h^2 + k^2; 
%% A = round orbit element / round F_t, F_r, F_n
A(1,1) = (2 * a^2)/sqrt( mu * p ) * p/r;
A(1,2) = (2 * a^2)/sqrt( mu * p ) * e*sin(nu);
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

D = -A' * QDotPartial;
end