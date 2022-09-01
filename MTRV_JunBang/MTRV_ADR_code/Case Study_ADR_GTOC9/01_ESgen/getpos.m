function [r,v] = getpos(q,t)
global mu
if(isempty(mu)), mu = 398600.4418;    end
d2r = pi/180;

a = q(1);
e = q(2);
i = q(3)*d2r;
omega = q(4)*d2r;
w = q(5)*d2r;
M0 = q(6)*d2r;

M = mod(M0+sqrt(mu/a^3)*t,2*pi);
if(e==0)
    nu = M;
else
    if(M>pi)
        E = M-e;
    else
        E = M+e;
    end
    while(1)
        E_new = E+(M-E+e*sin(E))/(1-e*cos(E));
        error = abs(E_new-E);
        if(error<1e-8)
            break;
        else
            E = E_new;
        end
    end
    sinnu = sin(E)*sqrt(1-e^2)/(1-e*cos(E));
    cosnu = (cos(E)-e)/(1-e*cos(E));
    nu = atan2(sinnu,cosnu);
end
p = a*(1-e^2);
r = p/(1+e*cos(nu));
h = sqrt(mu*p);

x = r*(cos(omega)*cos(w+nu)-sin(omega)*sin(w+nu)*cos(i));
y = r*(sin(omega)*cos(w+nu)+cos(omega)*sin(w+nu)*cos(i));
z = r*sin(w+nu)*sin(i);
xdot = x*h*e*sin(nu)/r/p-h*(cos(omega)*sin(w+nu)+sin(omega)*cos(w+nu)*cos(i))/r;
ydot = y*h*e*sin(nu)/r/p-h*(sin(omega)*sin(w+nu)-cos(omega)*cos(w+nu)*cos(i))/r;
zdot = z*h*e*sin(nu)/r/p+h*cos(w+nu)*sin(i)/r;

r = [x, y, z];
v = [xdot, ydot, zdot];

end