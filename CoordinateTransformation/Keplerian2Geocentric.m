function [X,Xdot] = Keplerian2Geocentric(a,e,i,o,w,nu)
% Keplerian Orbit Element 를 geocentric equatorial coordinate system 으로 변환


% use keplerian2ijk

% a : semimajor axis
% e : eccentricity
% i[rad] : inclination
% o[rad] : Longitude of the ascending node
% w[rad] : Argument of periapsis
% nu[rad] : True anomaly

mu = 398600.442;
p = a*(1 - e^2);
r = p / (1 + e*cos(nu));
h = sqrt(mu*p);

tmp1 = cos(w + nu);
tmp2 = sin(w + nu);

X=zeros(3,1);
Xdot=zeros(3,1);

X(1) = r*(cos(o)*tmp1 - sin(o)*tmp2*cos(i));
X(2) = r*(sin(o)*tmp1 + cos(o)*tmp2*cos(i));
X(3) = r*tmp2*sin(i);
Xdot(1) = X(1) * h*e*sin(nu) / r / p - h*(cos(o)*tmp2 + sin(o)*tmp1*cos(i)) / r;
Xdot(2) = X(2) * h*e*sin(nu) / r / p - h*(sin(o)*tmp2 - cos(o)*tmp1*cos(i)) / r;
Xdot(3) = X(3) * h*e*sin(nu) / r / p + h*tmp1*sin(i) / r;
end

