function [N, Vr11, Vt11, Vr12, Vt12, Vr21, Vt21, Vr22, Vt22] = VLAMB(Mu, r1, r2, theta, delT)

Vr11 = NaN; Vt11 = NaN; Vr12 = NaN; Vt12 = NaN;
Vr21 = NaN; Vt21 = NaN; Vr22 = NaN; Vt22 = NaN;
% return NaN for N<1 case

m = floor(theta/(2*pi)); %number of revolution
theta2 = theta/2 - m*pi; % normalized value of theta
cSqr = (r1-r2)^2 + 4*r1*r2*sin(theta2)^2; % c^2
c = sqrt(cSqr); % chord length
s = (r1 + r2 + c) /2; % 둘레
GMS = sqrt(Mu*s/2);
qSqrFM1 = c/s; %1-q^2
q = sqrt(r1*r2) * cos(theta2) /s;
if(c~=0)
    RHO = (r1-r2)/c;
    SIG = 4*r1*r2*sin(theta2)^2/cSqr;
else
    RHO=0;
    SIG=1;
end
T = 4*GMS*delT/s^2;
[N,X1,X2] = XLAMB(m,q,qSqrFM1,T);
% PROCEED FOR SINGLE SOLUTION, OR A PAIR
for I=1:N
    if (I==1)
        X=X1;
    else
        X=X2;
    end
    [~, QZMINX, QZPLX, ZPLQX]=TLAMB(m,q,qSqrFM1,X,-1);
    Vt2 = GMS*ZPLQX*sqrt(SIG);
    Vr1 = GMS*(QZMINX - QZPLX*RHO)/r1;
    Vt1 = Vt2/r1;
    Vr2 = -GMS*(QZMINX + QZPLX*RHO) /r2;
    Vt2 = Vt2/r2;
    if (I==1)
        Vr11 = Vr1;
        Vt11 = Vt1;
        Vr12 = Vr2;
        Vt12 = Vt2;
    else
        Vr21 = Vr1;
        Vt21 = Vt1;
        Vr22 = Vr2;
        Vt22 = Vt2;
    end
end

end

function [N,x0,x00] = XLAMB(m,q,qSqrFM1,TIn)
x0=NaN; % N=0 인 경우 X0 값으로 NaN 반환
x00=NaN; % N=1 인 경우 x00 값으로 NaN 반환
toler = 3*10^-7;
C0 = 1.7;
C1 = 0.5;
C2 = 0.03;
C3 = 0.15;
C41 = 1;
C42 = 0.24;
theta2 = atan2(qSqrFM1, 2*q)/pi;
if (m==0)
    %SINGLE REVOLUTION STARTER FROM T
    N=1;
    [T0,~,~,~]=TLAMB(m,q,qSqrFM1,0,0);
    TDiff = TIn - T0;
    if (TDiff <= 0)
        x0 = T0*TDiff / (-4+TIn);
    else
        x0 = -TDiff/(TDiff+4);
        W = x0+C0*sqrt(2*(1-theta2));
        if (W<0)
            x0 = x0 - sqrt ((-W)^(1/8))*(x0+sqrt(TDiff/(TDiff + 1.5*T0)));
        end
        W = 4 / (4+TDiff);
        x0 = x0 * (1+x0*(C1*W - C2*x0*sqrt(W)));
    end
else
    %multirevs
    xm = 1/(1.5*(m+0.5)*pi);
    if (theta2 < 0.5)
        xm = (2*theta2)^(1/8)*xm;
    end
    if (theta2 > 0.5 )
        xm = (2 * (2-2*theta2)^(1/8))*xm;
    end
    for i=1:12
        [TMin,dT,d2T,d3T] = TLAMB(m,q,qSqrFM1,xm,3);
        if (d2T == 0)
            break;
        end
        xmOld = xm;
        xm = xm - dT*d2T/(d2T*d2T-dT*d3T/2);
        xTest = abs(xmOld/xm - 1);
        if (xTest<=toler)
            break;
        end
        if (i==12)
            N=-1;
            return;
        end
    end
    TDiffM = TIn - TMin;
    if (TDiffM<0)
        N=0;
        return;
    elseif(TDiffM==0)
            x0 = xm;
            N = 1;
            return;
    else
        N = 3;
        if (d2T==0)
            d2T = 6*m*pi;
        end
        x0 = sqrt(TDiffM/(d2T/2 + TDiffM/(1-xm)^2));
        W = xm + x0;
        W = W*4/(4+TDiffM) + (1-W)^2;
        x0 = x0*(1-(1+m+C41*(theta2-0.5))/(1+C3*m)*x0*(C1*W + C2*x0*sqrt(W))) + xm;
        d2T2 = d2T/2;
        if (x0>=1)
            N=1;
            %go to 3
            [T0,dT,d2T,d3T] = TLAMB(m,q,qSqrFM1,0,0);
            TDiff0 = T0 - TMin;
            TDiff = TIn - T0;
            if (TDiff<=0)
                x0 = xm - sqrt(TDiffM/(d2T2-TDiffM*(d2T2/TDiff0-1/xm^2)));
            else
                x0= -TDiff/(TDiff+4);
                W = x0 + C0* sqrt(2*(1*theta2));
                if (W<0)
                    x0 = x0 - sqrt((-W)^(1/8))*(x0+sqrt(TDiff/(TDiff+1.5*T0)));
                end
                W = 4/(4+TDiff);
                x0 = x0*(1+(1+m+C42*(theta2-0.5))/(1+C3*m)*x0*(C1*W-C2*x0*sqrt(W)));
                if (x0<=-1)
                    N=N-1;
                    if (N==1)
                        x0 = x00;
                    end
                end
            end
        end
    end
end
%starting condition 완료
while(1)
    for i=1:3
        [T,dT,d2T,~] = TLAMB(m,q,qSqrFM1,x0,2);
        T = TIn-T;
        if (dT~=0)
            x0 = x0 + T*dT/ (dT*dT+T*d2T/2);
        end
    end
    if (N~=3)
        return; % Exit if only one solution, normally when M=0
    end
    N=2;
    x00=x0;
    [T0,~,~,~] = TLAMB(m,q,qSqrFM1,0,0);
    TDiff0 = T0 - TMin;
    TDiff = TIn - T0;
    if (TDiff<=0)
        x0 = xm - sqrt(TDiffM/(d2T2-TDiffM*(d2T2/TDiff0-1/xm^2)));
    else
        x0= -TDiff/(TDiff+4);
        W = x0 + C0* sqrt(2*(1*theta2));
        if (W<0)
            x0 = x0 - sqrt((-W)^(1/8))*(x0+sqrt(TDiff/(TDiff+1.5*T0)));
        end
        W = 4/(4+TDiff);
        x0 = x0*(1+(1+m+C42*(theta2-0.5))/(1+C3*m)*x0*(C1*W-C2*x0*sqrt(W)));
        if (x0<=-1)
            N=N-1;
            if (N==1)
                x0 = x00;
            end
        end
    end
end
end


function [T,dT,d2T,d3T] = TLAMB(m,q,qSqrFM1,x,nDot)
T=NaN; % N=-1 인 경우, T= NaN 으로 출력
LM1 = (nDot==-1);
L1 = (nDot>=1);
L2 = (nDot>=2);
L3 = (nDot==3);
qSqr= q*q;
xSqr = x*x;
U = (1-x)*(1+x); % -E
if (~LM1)
    %needed if series, and otherwise usful when z =0
    dT = 0;
    d2T = 0;
    d3T = 0;
end
if (LM1 || (m>0) || (x<0) || (abs(U)>0.4))
    %direct computation (not series)
    y = sqrt(abs(U));
    z = sqrt(qSqrFM1 + qSqr*xSqr);
    qx = q*x;
    if (qx<=0)
        A = z - qx;
        B = q*z - x;
    end
    if (qx<0 && LM1)
        AA = qSqrFM1 / A;
        BB = qSqrFM1 * (qSqr*U - xSqr)/B;
    end
    if (((qx==0) && LM1) || (qx>0))
        AA = z + qx;
        BB = q*z + x;
    end
    if (qx>0)
        A = qSqrFM1/AA;
        B = qSqrFM1*(qSqr*U - xSqr) /BB;
    end
    if (~LM1)
        if (qx*U >= 0)
            G = x*z + q*U;
        else
            G = (xSqr - qSqr*U) / (x*z - q*U);
        end
        F = A*y;
        if (x<=1)
            T = m*pi + atan2(F,G);
        else
            if(F>0.4)
                T = log(F+G);
            else
                FG1 = F/(G+1);
                TERM = 2*FG1;
                FG1Sqr = FG1*FG1;
                T = TERM;
                TWOI1 = 1;
                TWOI1 = TWOI1 + 2;
                TERM = TERM * FG1Sqr;
                TOld = T;
                T = T + TERM/TWOI1;
                while(T~=TOld)
                    TWOI1 = TWOI1 + 2;
                    TERM = TERM * FG1Sqr;
                    TOld = T;
                    T = T + TERM/TWOI1;
                end
            end
        end
        T = 2 * (T/y + B)/U;
        if (L1 && (z~=0))
            dT = (3*x*T - 4*(A + qx*qSqrFM1) / z) /U;
            if (L2)
                d2T= ( 3*T + 5*x*dT + 4*(q/z)^3*qSqrFM1 ) /U;
            end
            if (L3)
                d3T = (8*dT+7*x*d2T-12*(q/z)^5*x*qSqrFM1)/U;
            end
        end
    else
        dT=B;
        d2T = BB;
        d3T = AA;
    end
else
    %compute by series
    u0=1;
    if (L1)
        u1 = 1;
    end
    if (L2)
        u2 = 1;
    end
    if (L3)
        u3 = 1;
    end
    TERM = 4;
    TQ = q*qSqrFM1;
    I = 0;
    if (q<5/10)
        TQSUM = 1 - q*qSqr;
    end
    if (q>=5/10)
        TQSUM = (1/(1+q) + q)*qSqrFM1;
    end
    TTMOLD = TERM/3;
    T = TTMOLD*TQSUM;
    % Start of Loop
    while (1)
        I = I + 1;
        P = I;
        u0 = u0*U;
        if (L1 && (I>1))
            u1 = u1*U;
        end
        if (L2 && (I>2))
            u2 = u2*U;
        end
        if (L3 && (I>3))
            u3 = u3*U;
        end
        TERM = TERM * (P-0.5)/P;
        TQ = TQ*qSqr;
        TQSUM = TQSUM + TQ;
        TOld = T;
        TTERM = TERM/(2*P + 3);
        TQTERM = TTERM * TQSUM;
        T = T - u0*((1.5*P + 0.25)*TQTERM / (P*P - 0.25 - TTMOLD*TQ));
        TTMOLD = TTERM;
        TQTERM = TQTERM*P;
        if (L1)
            dT = dT+ TQTERM*u1;
        end
        if (L2)
            d2T = d2T+ TQTERM*u2*(P-1);
        end
        if (L3)
            d3T = d3T+ TQTERM*u3*(P-1)*(P-2);
        end
        if ~((I<nDot) || (T~=TOld))
            break;
        end
    end
    if (L3)
        d3T = 8*x*(1.5*d2T - xSqr * d3T);
    end
    if (L2)
        d2T = 2*(2*xSqr*d2T - dT);
    end
    if (L1)
        dT = -2D0*x*dT;
    end
    T = T/xSqr;
end
end
