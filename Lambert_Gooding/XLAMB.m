function [N,X,XPL] = XLAMB(M,Q,QSQFM1,TIN)

TOL = 3*10^-7;
C0 = 1.7;
C1 = 0.5;
C2 = 0.03;
C3 = 0.15;
C41 = 1;
C42 = 0.24;
THR2 = atan2(QSQFM1, 2*Q)/pi;
if (M==0)
    %SINGLE REVOLUTION STARTER FROM T
    N=1;
    [T0,~,~,~]=TLAMB(M,Q,QSQFM1,0,0);
    TDIFF = TIN - T0;
    if (TDIFF <= 0)
        X = T0*TDIFF / (-4+TIN);
    else
        X = -TDIFF/(TDIFF+4);
        W = X+C0*sqrt(2*(1-THR2));
        if (W<0)
            X = X - sqrt (D8RT(-W))*(X+sqrt(TDIFF/(TDIFF + 1.5*T0)));
        end
        W = 4 / (4+TDIFF);
        X = X * (1+X*(C1*W - C2*X*sqrt(W)));
    end
else
    %multirevs
    XM = 1/(1.5*(M+0.5)*PI);
    if (THR2 < 0.5)
        XM = D8RT(2*THR2)*XM;
    end
    if (THR > 0.5 )
        XM = (2 * D8RT(2-2*THR2))*XM;
    end
    for I=1:12
        [TMIN,DT,D2T,D3T] = TLAMB(M,Q,QSQFM1,XM,3);
        if (D2T == 0)
            break;
        end
        XMOLD = XM;
        XM = XM - DT*D2T/(D2T*D2T-DT*D3T/2);
        XTEST = abs(XMOLD/XM - 1);
        if (XTEST<=TOL)
            break;
        end
        if (I==12)
            N=-1;
            return;
        end
    end
    TDIFFM = TIN - TMIN;
    if (TDIFFM<0)
        N=0;
        return;
    elseif(TDIFFM==0)
            X = XM;
            N = 1;
            return;
    else
        N = 3;
        if (D2T==0)
            D2T = 6*M*pi;
        end
        X = sqrt(TDIFFM/(D2T/2 + TDIFFM/(1-XM)^2));
        W = XM + X;
        W = W*4/(4+TDIFFM) + (1-W)^2;
        X = X*(1-(1+M+C41*(THR2-0.5))/(1+C3*M)*X*(c1*W + C2*X*sqrt(W))) + XM;
        D2T2 = D2T/2;
        if (X>=1)
            N=1;
        end
    end
end
%starting condition 완료
if (X>=1)
    %논문의 GOTO 3 구현..
    [T0,~,~,~] = TLAMB(M,Q,QSQFM1,0,0);
    TDIFF0 = T0 - TMIN;
    TDIFF = TIN - T0;
    if (TDIFF<=0)
        X = XM - sqrt(TDIFFM/(D2T2-TDIFFM*(D2T2/TDIFF0-1/XM^2)));
    else
        X= -TDIFF/(TDIFF+4);
        W = X + C0* sqrt(2*(1*THR2));
        if (W<0)
            X = X - sqrt(D8RT(-W))*(X+sqrt(TDIFF/(TDIFF+1.5*T0)));
        end
        W = 4/(4+TDIFF);
        X = X*(1+(1+M+C42*(THR2-0.5))/(1+C3*M)*X*(C1*W-C2*X*sqrt(W)));
        if (X<=-1)
            N=N-1;
            if (N==1)
                X = XPL;
            end
        end
    end
end

for i=1:3
    [T,DT,D2T,~] = TLAMB(M,Q,QSQFM1,X,2);
    T = TIN - T;
    if (DT~=0)
        X = X + T*DT/(DT*DT + T*D2T/2);
    end
end
if(N~=3)
    return;
end
N=2;
XPL = X;
[T0,~,~,~] = TLAMB(M,Q,QSQFM1,0,0);
TDIFF0 = T0 - TMIN;
TDIFF = TIN - T0;
if (TDIFF<=0)
    X = XM - sqrt(TDIFFM/(D2T2-TDIFFM*(D2T2/TDIFF0-1/XM^2)));
else
    X= -TDIFF/(TDIFF+4);
    W = X + C0* sqrt(2*(1*THR2));
    if (W<0)
        X = X - sqrt(D8RT(-W))*(X+sqrt(TDIFF/(TDIFF+1.5*T0)));
    end
    W = 4/(4+TDIFF);
    X = X*(1+(1+M+C42*(THR2-0.5))/(1+C3*M)*X*(C1*W-C2*X*sqrt(W)));
    if (X<=-1)
        N=N-1;
        if (N==1)
            X = XPL;
        end
    end
end

end


function d8rt = D8RT(x) 
d8rt= sqrt(sqrt(sqrt(x)));
end