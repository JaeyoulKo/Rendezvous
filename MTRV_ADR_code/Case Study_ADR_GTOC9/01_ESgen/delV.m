%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Rendezvous cost based on three-step transfer using RAAN drift         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [J_opt,r_dr_opt,i_dr_opt] = delV(i,j,td,ta)
global t_d t_a
global mu re d2r J2term Q
global ai_d ei_d ii_d omi_0 omi_d 
global aj_a ej_a ij_a omj_0 omj_a
global count
count = count+1;

J_opt = 100; r_dr_opt = NaN; i_dr_opt = NaN;
if(td>=ta), return; end

%% Input parameters
minalt = 222;    % minimum altitude
maxalt = 1000;   % maximum altitude

t_d = td; t_a = ta;
mu = 398600.4418;
re = 6378.137;
d2r = pi/180;
J2term = -1.5*1.08261e-3*sqrt(mu)*re^2*24*3600/d2r;  % deg/day
ai_d = Q(i,1); ei_d = Q(i,2); ii_d = Q(i,3); omi_0 = Q(i,4);
aj_a = Q(j,1); ej_a = Q(j,2); ij_a = Q(j,3); omj_0 = Q(j,4);
omdoti = J2term*cos(ii_d*d2r)/ai_d^3.5/(1-ei_d^2)^2;    omi_d = mod(omi_0+omdoti*t_d,360);
omdotj = J2term*cos(ij_a*d2r)/aj_a^3.5/(1-ej_a^2)^2;    omj_a = mod(omj_0+omdotj*t_a,360);
dom = mod(omj_a-omi_d,360);
checkf = 1;
if(dom/(ta-td)>-J2term*(re+minalt)^-3.5)
    dom = dom-360;
    if(dom/(ta-td)<J2term*(re+minalt)^-3.5)
        checkf = 0;
    end
end
if(~checkf), return; end

%% Find optimal drift orbit (r_dr,i_dr)
options = optimoptions('fmincon','Algorithm','sqp','Display','off','TolX',1e-10);
tmp1 = dom*(re+minalt)^3.5/J2term/(ta-td);
tmp2 = dom*(re+maxalt)^3.5/J2term/(ta-td);
if(tmp1>1), tmp1=1; end;  if(tmp1<-1), tmp1=-1; end
if(tmp2>1), tmp2=1; end;  if(tmp2<-1), tmp2=-1; end
tmp1 = acos(tmp1)/d2r;
tmp2 = acos(tmp2)/d2r;
if(tmp1>90)
    [x1,J1] = fmincon(@myfun,[re+0.5*(minalt+maxalt);tmp1],[],[],[],[],[re+minalt;tmp1],[re+maxalt;tmp2],@mycon,options);
else
    [x1,J1] = fmincon(@myfun,[re+0.5*(minalt+maxalt);tmp1],[],[],[],[],[re+minalt;tmp2],[re+maxalt;tmp1],@mycon,options);
end
J_opt = J1; r_dr_opt = x1(1); i_dr_opt = x1(2);
end

function J = myfun(x) % x=[r_dr;i_dr]
global mu d2r J2term
global t_d t_a
global ai_d ei_d ii_d omi_d aj_a ej_a ij_a omj_a
r_dr = x(1); i_dr = x(2);
omdotdr = J2term*cos(i_dr*d2r)/r_dr^3.5;
omdr_a = omi_d+omdotdr*(t_a-t_d);

%% Hohmann to drift orbit
r = ai_d*(1+ei_d); deli = abs(i_dr-ii_d);
vtmp1 = sqrt(2*mu/r-mu/ai_d); vtmp2 = sqrt(2*mu/r-2*mu/(r+r_dr)); 
vtmp3 = sqrt(2*mu/r_dr-2*mu/(r+r_dr)); vtmp4 = sqrt(mu/r_dr);
if(r_dr>r)
    tmp1 = abs(vtmp2-vtmp1);
    tmp2 = abs(vtmp4*cos(deli*d2r)-vtmp3);
    tmp3 = abs(vtmp4*sin(deli*d2r));
    delv1 = tmp1+sqrt(tmp2*tmp2+tmp3*tmp3);
else
    tmp1 = abs(vtmp2*cos(deli*d2r)-vtmp1);
    tmp2 = abs(vtmp2*sin(deli*d2r));
    tmp3 = abs(vtmp4-vtmp3);
    delv1 = sqrt(tmp1*tmp1+tmp2*tmp2)+tmp3;
end

%% Hohmann to debris orbit
r = aj_a*(1+ej_a); % deli = abs(i_dr-ij_a);
deli = acos(cos(i_dr*d2r)*cos(ij_a*d2r)+sin(i_dr*d2r)*sin(ij_a*d2r)*cos(omdr_a*d2r-omj_a*d2r))/d2r;
vtmp1 = sqrt(mu/r_dr); vtmp2 = sqrt(2*mu/r_dr-2*mu/(r+r_dr)); 
vtmp3 = sqrt(2*mu/r-2*mu/(r+r_dr)); vtmp4 = sqrt(2*mu/r-mu/aj_a);
if(r_dr<r)
    tmp1 = abs(vtmp2-vtmp1);
    tmp2 = abs(vtmp4*cos(deli*d2r)-vtmp3);
    tmp3 = abs(vtmp4*sin(deli*d2r));
    delv2 = tmp1+sqrt(tmp2*tmp2+tmp3*tmp3);
else
    tmp1 = abs(vtmp2*cos(deli*d2r)-vtmp1);
    tmp2 = abs(vtmp2*sin(deli*d2r));
    tmp3 = abs(vtmp4-vtmp3);
    delv2 = sqrt(tmp1*tmp1+tmp2*tmp2)+tmp3;
end

J = delv1+delv2;
end

function [c,ceq] = mycon(x)
global d2r J2term
global t_d t_a
global omi_d omj_a 
r_dr = x(1); i_dr = x(2);

c=[];
omdotdr = J2term*cos(i_dr*d2r)/r_dr^3.5;
omdr_a = omi_d+omdotdr*(t_a-t_d);
tmp1 = mod(omj_a-omdr_a,360);
if(tmp1>180), tmp1=tmp1-360; end
ceq=tmp1;
end