function Q = CalcQ(chaserState, targetState, F, mu, WeightParameter)
% Q law Guidance에 사용되는 Q 값을 계산하는 함수.
% WeightParameter = [Wa; Wf; Wg; Wh; Wk] / norm([Wa; Wf; Wg; Wh; Wk]);
% chaser & target State = [p;f;g;h;k;L;Mass(not use)]
p=chaserState(1);
f=chaserState(2);
g=chaserState(3);
sSqr  = 1 + chaserState(4)^2 + chaserState(5)^2;

%% convert p to a
chaserState(1) = chaserState(1) / (1 - chaserState(2)^2 - chaserState(3)^2);
targetState(1) = targetState(1) / (1 - targetState(2)^2 - targetState(3)^2);

%now State = [a;f;g;h;k;L]
a=chaserState(1);

%% maximum rates of change of the orbit elements
aDotxx = 2 * F * a * sqrt(a/mu) * sqrt( ( 1+sqrt(f^2+g^2) ) / ( 1-sqrt(f^2+g^2) ) );
fDotxx = 2 * F * sqrt(p/mu);
gDotxx = 2 * F * sqrt(p/mu);
hDotxx = 0.5 * F * sqrt(p/mu) * sSqr / ( sqrt(1-g^2)+f );
kDotxx = 0.5 * F * sqrt(p/mu) * sSqr / ( sqrt(1-f^2)+g );
% LDotxx = F * sqrt(p/mu) * sqrt(sSqr-1);

maxRateOfChange=[aDotxx;fDotxx;gDotxx;hDotxx;kDotxx];
% maxRateOfChange=[aDotxx;fDotxx;gDotxx;hDotxx;kDotxx;LDotxx];

%% multiply Scaling Factor S to WeightParameter
% m=3; n=4; r=2;
% WeightParameter(1) = WeightParameter(1) * (1+(abs(chaserState(1)-targetState(1))/targetState(1)/m)^n)^(1/r);

% r_pmin = 6600;
% k=10;
% P = exp(k*(1-p/(1+e)/r_pmin));

%% Calculate Q value
% Q =sum ((1+Wp*P).*WeightParameter.*((chaserState(1:6)-targetState(1:6))./maxRateOfChange).^2);
diffState = chaserState(1:5)-targetState(1:5);
% diffState = chaserState(1:6)-targetState(1:6);
% diffState(6) = rem(diffState(6),2*pi);
% if diffState(6)>pi
%     diffState(6) = diffState(6)-2*pi;
% end
Q =sum (WeightParameter.*((diffState)./maxRateOfChange).^2);

end

