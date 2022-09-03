function nu = TrueAnomalyAtdeltaT(a,e,i,o,w,nu0,t)

mu = 398600.442;
M0 = nu0 - 2 * e*sin(nu0) + (3 * e^2 / 4 + e^4 / 8)*sin(2 * nu0) - e^3*sin(3 * nu0) / 3 + 5*e^4*sin(4 * nu0) / 32;
% M0 = nu0;
M = M0 + sqrt(mu / a^3)*t;
% Mean anomaly - (Measured in Radian) - Mean anomaly is the fraction of an elliptical orbit's period that has elapsed since the orbiting body passed periapsis.
if (M < 0)
	M = M + 2 * pi;
end
if (M > 2 * pi)
	M = mod(M, 2 * pi);
end
if (e == 0)
	nu = M;
else
    if (M > pi)
	    E = M - e;
    else
	    E = M + e;
    end
    while(1)
	    E_new = E + (M - E + e*sin(E)) / (1 - e*cos(E));
	    err = abs(E_new - E);
        if (err < 1e-8)
	        break;
        else
	        E = E_new;
        end
    end
	sinNu = sin(E)*sqrt(1 - e^2) / (1 - e*cos(E));
	cosNu = (cos(E) - e) / (1 - e*cos(E));
	nu = atan2(sinNu, cosNu);
end

