function dt = TimeStepForQLaw(chaserState,Q,QDot,mu)
% Q Law Guidauance Input을 Update할 Time Step을 결정
% chaser State = [p;f;g;h;k;L;Mass(not use)]
% dL : hyperparameter : 튜닝요소
dL = 10 * pi/180;
% addpath("./../CoordinateTransformation/")
r = Equinoctial2Distance(chaserState(1),chaserState(2),chaserState(3),chaserState(4),chaserState(5),chaserState(6));
dt = dL * sqrt(r^3/mu);
if Q + QDot*dt < 0
    dt = - Q / QDot;
end
end

