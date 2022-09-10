clear all
clc
addpath("./../CoordinateTransformation/");
tic
%% gravitational constant and spacecraft 
%--------------- Earth ------------------%
mu        = 3.986004418*10^5;       % km^3/s^2
Re        = 6378.1;                 % km
ThrustMax = 1;                      % N
Mmax      = 300;                    % kg
Isp       = 3100;                   % sec
g0        = 9.80665;                % m/s^2
c         = Isp*g0;
F         = ThrustMax / Mmax;
J2        = 0;
%--------------- Sun ------------------%
% mu        = 1.32712440018*10^11;    % km^3/s^2
% Re        = 149597870.7;            % km
% Thrustmax = 0.5;                    % N
% Mmax      = 2000;                   % kg
% Isp       = 2000;                   % sec
% g0        = 9.80665;                % m/s^2
% c         = Isp*g0;
% F         = Tmax / Mmax;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AU2km   = Re;
TU2sec  = sqrt(AU2km^3 / mu);
TU2days = TU2sec / (60*60*24);
MUkg    = Mmax;

mu   = mu   / (AU2km^3/TU2sec^2);
Re   = Re   / AU2km;
ThrustMax = ThrustMax / (MUkg * AU2km * 1000 / TU2sec^2); %?
Mmax = Mmax / MUkg; %?
Isp  = Isp  / TU2sec;
g0   = g0   / (AU2km*1000 / TU2sec^2);
c    = Isp * g0;
F    = ThrustMax / Mmax;

%% Weight parameters
Wa = 1;
Wf = 1;
Wg = 1;
Wh = 1;
Wk = 1;

W = [Wa; Wf; Wg; Wh; Wk] / norm([Wa; Wf; Wg; Wh; Wk]);

%% chaser initial state 
aChaser = 7000/AU2km;
% aChaser = 7000;
eChaser = 0.01;
iChaser = 0.05*pi/180;
oChaser = 0;
wChaser = 0;
nuChaser  = 0;
initialMass = Mmax;
[pChaser,fChaser,gChaser,hChaser,kChaser,LChaser] = Keplerian2Equinoctial(aChaser,eChaser,iChaser,wChaser,oChaser,nuChaser);
chaserState = [pChaser; fChaser; gChaser; hChaser; kChaser; LChaser; initialMass];

aTarget = 42000/AU2km;
% aTarget = 42000;
eTarget = 0.01;
iTarget = 0.05*pi/180;
oTarget = 0;
wTarget = 0;
nuTarget  = 0;
[pTarget,fTarget,gTarget,hTarget,kTarget,LTarget] = Keplerian2Equinoctial(aTarget,eTarget,iTarget,wTarget,oTarget,nuTarget);

targetState = [pTarget; fTarget; gTarget; hTarget; kTarget; LTarget; nan]; %? initial Mass .... ?

%% initialize
timeStepUnitSize=3000;
time = zeros(timeStepUnitSize,1);
QHist = zeros(timeStepUnitSize,1);
chaserTrajectory = zeros(7,timeStepUnitSize);
targetTrajectory = zeros(7,timeStepUnitSize);

chaserTrajectory(:,1) = chaserState;
targetTrajectory(:,1) = targetState;

Q = CalcQ(chaserState, targetState, F, mu, W);
QHist(1) = Q;

%% Q-Guidance Simulation Start
i = 1;
tol = 100;
while 1
    
    F = ThrustMax / chaserState(7);
    i = i+1;
    if (rem(i,timeStepUnitSize)==0)
        time = [time ; zeros(timeStepUnitSize,1)];
        QHist = [QHist ; zeros(timeStepUnitSize,1)];
        chaserTrajectory = [chaserTrajectory, zeros(7,timeStepUnitSize)];
        targetTrajectory = [targetTrajectory, zeros(7,timeStepUnitSize)];
    end
    [QDot, controlInput] = GetQGuidance(chaserState, targetState, F, mu, W);
    dt = TimeStepForQLaw(chaserState,Q,QDot,mu);
    time(i) = time(i-1) + dt;    
    %% chaser trajectory
    inputWithPerturb=controlInput+J2PerturbationDynamicsEquinoctial(chaserState, [mu;Re;J2]);
%     odeOption = odeset('MaxStep',1);
%     [time , st] = ode45(@(time,state)TwoBodyVOPDynamicsEquinoctial(state,inputWithPerturb,[mu;ThrustMax;c]),[0 dt],chaserState,odeOption);
    chaserState = RK4IntegralGaussEq(chaserState,inputWithPerturb,dt,[mu;ThrustMax;c]);
    chaserTrajectory(:,i) = chaserState;

    %% target trajectory
    targetPerturbation = J2PerturbationDynamicsEquinoctial(targetState, [mu;Re;J2]);
    targetState = RK4IntegralGaussEq(targetState,targetPerturbation,dt,[mu;ThrustMax;c]);
    targetTrajectory(:,i) = targetState;

    Q = CalcQ(chaserState, targetState, F, mu, W);
    QHist(i) = Q;
    if Q < tol
        time = time(1:i);
        QHist = QHist(1:i);
        chaserTrajectory = chaserTrajectory(:,1:i);
        targetTrajectory = targetTrajectory(:,1:i);
        break
    end
    
end
toc
%% data for chaser orbit plot
pChaser = ones(1,100) * pChaser;
fChaser = ones(1,100) * fChaser;
gChaser = ones(1,100) * gChaser;
hChaser = ones(1,100) * hChaser;
kChaser = ones(1,100) * kChaser;
LChaser = linspace(LChaser,LChaser+2*pi,100);

[aChaser,eChaser,iChaser,oChaser,wChaser,nuChaser] = Equinoctial2KeplerianMIMO(pChaser,fChaser,gChaser,hChaser,kChaser,LChaser);
[rChaserOrbit,~] = Keplerian2GeocentricMIMO(pChaser,eChaser,iChaser,oChaser,wChaser,nuChaser);

%% data for target orbit plot
pTarget = ones(1,100) * pTarget;
fTarget = ones(1,100) * fTarget;
gTarget = ones(1,100) * gTarget;
hTarget = ones(1,100) * hTarget;
kTarget = ones(1,100) * kTarget;
LTarget = linspace(LTarget,LTarget+2*pi,100);

[aTarget,eTarget,iTarget,oTarget,wTarget,nuTarget] = Equinoctial2KeplerianMIMO(pTarget,fTarget,gTarget,hTarget,kTarget,LTarget);
[rTargetOrbit,~] = Keplerian2GeocentricMIMO(pTarget,eTarget,iTarget,oTarget,wTarget,nuTarget);

%% Low Thrust RV trajectory
[a,e,i,o,w,nu] = Equinoctial2KeplerianMIMO(chaserTrajectory(1,:), chaserTrajectory(2,:), chaserTrajectory(3,:), chaserTrajectory(4,:), chaserTrajectory(5,:) ,chaserTrajectory(6,:));

[r,~] = Keplerian2GeocentricMIMO(a,e,i,o,w,nu);

%% plot
figure
subplot(2,3,1)
hold on
plot3(rChaserOrbit(1,:),rChaserOrbit(2,:),rChaserOrbit(3,:),'b','LineWidth',3)
plot3(rTargetOrbit(1,:),rTargetOrbit(2,:),rTargetOrbit(3,:),'r','LineWidth',3)
plot3(r(1,:),r(2,:),r(3,:),'k')
plot3(r(1,1),r(2,1),r(3,1),'bd','MarkerSize',10,'LineWidth',3)
plot3(r(1,end),r(2,end),r(3,end),'rs','MarkerSize',10,'LineWidth',3)
grid on
xlabel('x')
ylabel('y')
zlabel('z')
legend('initial orbit','final orbit','transfer orbit')

pp1 = chaserTrajectory(1,:);
ff1 = chaserTrajectory(2,:);
gg1 = chaserTrajectory(3,:);
hh1 = chaserTrajectory(4,:);
kk1 = chaserTrajectory(5,:);
LL1 = chaserTrajectory(6,:);

[aa1,ecc1,inc1,Ome1,ome1,nu1] = lowThrustMee2Coe(pp1, ff1, gg1, hh1, kk1, LL1);

pp2 = targetTrajectory(1,:);
ff2 = targetTrajectory(2,:);
gg2 = targetTrajectory(3,:);
hh2 = targetTrajectory(4,:);
kk2 = targetTrajectory(5,:);
LL2 = targetTrajectory(6,:);

[aa2,ecc2,inc2,Ome2,ome2,nu2] = lowThrustMee2Coe(pp2, ff2, gg2, hh2, kk2, LL2);

subplot(2,3,2)
plot(time*TU2days,aa1*AU2km)
hold on
plot(time*TU2days,aa2*AU2km)
xlabel('time (days)')
ylabel('semi-major axis (km)')

subplot(2,3,3)
plot(time*TU2days,ecc1)
hold on
plot(time*TU2days,ecc2)
xlabel('time (days)')
ylabel('eccentricity (-)')

subplot(2,3,4)
plot(time*TU2days,inc1*180/pi)
hold on
plot(time*TU2days,inc2*180/pi)
xlabel('time (days)')
ylabel('inclination (deg)')

subplot(2,3,5)
plot(time*TU2days,Ome1*180/pi)
hold on
plot(time*TU2days,Ome2*180/pi)
xlabel('time (days)')
ylabel('RAAN (deg)')

subplot(2,3,6)
plot(time*TU2days,ome1*180/pi)
hold on
plot(time*TU2days,ome2*180/pi)
xlabel('time (days)')
ylabel('argument of perigee (deg)')