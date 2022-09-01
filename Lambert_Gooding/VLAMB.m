function [N, VR11, VT11, VR12, VT12, VR21, VT21, VR22, VT22] = VLAMB(GM, R1, R2, TH, TDELT)

VR11 = NaN; VT11 = NaN; VR12 = NaN; VT12 = NaN;
VR21 = NaN; VT21 = NaN; VR22 = NaN; VT22 = NaN;
% return NaN for N<1 case

TWOPI=2*pi;
M = floor(TH/TWOPI);
THR2 = TH/2 - M*pi;
DR = R1-R2;
R1R2 = R1*R2;
R1R2TH = 4*R1R2*sin(THR2)^2;
CSQ = DR^2 + R1R2TH;
C = sqrt(CSQ);
S = (R1 + R2 + C) /2;
GMS = sqrt(GM*S/2);
QSQFM1 = C/S;
Q = sqrt(R1R2) * cos(THR2) /S;
if(C~=0)
    RHO = DR/C;
    SIG = R1R2TH/CSQ;
else
    RHO=0;
    SIG=1;
end
T = 4*GMS*TDELT/S^2;
[N,X1,X2] = XLAMB(M,Q,QSQFM1,T);
% PROCEED FOR SINGLE SOLUTION, OR A PAIR
for I=1:N
    if (I==1)
        X=X1;
    else
        X=X2;
    end
    [~, QZMINX, QZPLX, ZPLQX]=TLAMB(M,Q,QSQFM1,X,-1);
    VT2 = GMS*ZPLQX*sqrt(SIG);
    VR1 = GMS*(QZMINX - QZPLX*RHO)/R1;
    VT1 = VT2/R1;
    VR2 = -GMS*(QZMINX + QZPLX*RHO) /R2;
    VT2 = VT2/R2;
    if (I==1)
        VR11 = VR1;
        VT11 = VT1;
        VR12 = VR2;
        VT12 = VT2;
    else
        VR21 = VR1;
        VT21 = VT1;
        VR22 = VR2;
        VT22 = VT2;
    end
end

end

function [N,X,XPL] = XLAMB(M,Q,QSQFM1,TIN)
X=NaN; % N=0 인 경우 X 값으로 NaN 반환
XPL=NaN; % N=1 인 경우 XPL 값으로 NaN 반환
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
    XM = 1/(1.5*(M+0.5)*pi);
    if (THR2 < 0.5)
        XM = D8RT(2*THR2)*XM;
    end
    if (THR2 > 0.5 )
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
        X = X*(1-(1+M+C41*(THR2-0.5))/(1+C3*M)*X*(C1*W + C2*X*sqrt(W))) + XM;
        D2T2 = D2T/2;
        if (X>=1)
            N=1;
            %go to 3
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
    end
end
%starting condition 완료
while(1)
    for I=1:3
        [T,DT,D2T,D3T] = TLAMB(M,Q,QSQFM1,X,2);
        T = TIN-T;
        if (DT~=0)
            X = X + T*DT/ (DT*DT+T*D2T/2);
        end
    end
    if (N~=3)
        return; % Exit if only one solution, normally when M=0
    end
    N=2;
    XPL=X;
    %논문의 3 구현..
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

function [T,DT,D2T,D3T] = TLAMB(M,Q,QSQFM1,X,N)
T=NaN; % N=-1 인 경우, T= NaN 으로 출력
SW = 0.4;
LM1 = (N==-1);
L1 = (N>=1);
L2 = (N>=2);
L3 = (N==3);
QSQ = Q*Q;
XSQ = X*X;
U = (1-X)*(1+X); % -E
if (~LM1)
    %needed if series, and otherwise usful when z =0
    DT = 0;
    D2T = 0;
    D3T = 0;
end
if (LM1 || (M>0) || (X<0) || (abs(U)>SW))
    %direct computation (not series)
    Y = sqrt(abs(U));
    Z = sqrt(QSQFM1 + QSQ*XSQ);
    QX = Q*X;
    if (QX<=0)
        A = Z - QX;
        B = Q*Z - X;
    end
    if (QX<0 && LM1)
        AA = QSQFM1 / A;
        BB = QSQFM1 * (QSQ*U - XSQ)/B;
    end
    if (((QX==0) && LM1) || (QX>0))
        AA = Z + QX;
        BB = Q*Z + X;
    end
    if (QX>0)
        A = QSQFM1/AA;
        B = QSQFM1*(QSQ*U - XSQ) /BB;
    end
    if (~LM1)
        if (QX*U >= 0)
            G = X*Z + Q*U;
        else
            G = (XSQ - QSQ*U) / (X*Z - Q*U);
        end
        F = A*Y;
        if (X<=1)
            T = M*pi + atan2(F,G);
        else
            if(F>SW)
                T = log(F+G);
            else
                FG1 = F/(G+1);
                TERM = 2*FG1;
                FG1SQ = FG1*FG1;
                T = TERM;
                TWOI1 = 1;
                TWOI1 = TWOI1 + 2;
                TERM = TERM * FG1SQ;
                TOLD = T;
                T = T + TERM/TWOI1;
                while(T~=TOLD)
                    TWOI1 = TWOI1 + 2;
                    TERM = TERM * FG1SQ;
                    TOLD = T;
                    T = T + TERM/TWOI1;
                end
            end
        end
        T = 2 * (T/Y + B)/U;
        if (L1 && (Z~=0))
            QZ = Q/Z;
            QZ2 = QZ*QZ;
            QZ = QZ*QZ2;
            DT = (3*X*T - 4*(A + QX*QSQFM1) / Z) /U;
            if (L2)
                D2T= ( 3*T + 5*X*DT + 4*QZ*QSQFM1 ) /U;
            end
            if (L3)
                D3T = (8*DT+7*X*D2T-12*QZ*QZ2*X*QSQFM1)/U;
            end
        end
    else
        DT=B;
        D2T = BB;
        D3T = AA;
    end
else
    %compute by series
    U0I=1;
    if (L1)
        U1I = 1;
    end
    if (L2)
        U2I = 1;
    end
    if (L3)
        U3I = 1;
    end
    TERM = 4;
    TQ = Q*QSQFM1;
    I = 0;
    if (Q<5/10)
        TQSUM = 1 - Q*QSQ;
    end
    if (Q>=5/10)
        TQSUM = (1/(1+Q) + Q)*QSQFM1;
    end
    TTMOLD = TERM/3;
    T = TTMOLD*TQSUM;
    % Start of Loop
    I = I + 1;
    P = I;
    U0I = U0I*U;
    if (L1 && (I>1))
        U1I = U1I*U;
    end
    if (L2 && (I>2))
        U2I = U2I*U;
    end
    if (L3 && (I>3))
        U3I = U3I*U;
    end
    TERM = TERM * (P-0.5)/P;
    TQ = TQ*QSQ;
    TQSUM = TQSUM + TQ;
    TOLD = T;
    TTERM = TERM/(2*P + 3);
    TQTERM = TTERM * TQSUM;
    T = T - U0I*((1.5*P + 0.25)*TQTERM / (P*P - 0.25) - TTMOLD*TQ);
    TTMOLD = TTERM;
    TQTERM = TQTERM*P;
    if (L1)
        DT = DT+ TQTERM*U1I;
    end
    if (L2)
        D2T = D2T+ TQTERM*U2I*(P-1);
    end
    if (L3)
        D3T = D3T+ TQTERM*U3I*(P-1)*(P-2);
    end
    while ((I<N) || (T~=TOLD))
        I = I + 1;
        P = I;
        U0I = U0I*U;
        if (L1 && (I>1))
            U1I = U1I*U;
        end
        if (L2 && (I>2))
            U2I = U2I*U;
        end
        if (L3 && (I>3))
            U3I = U3I*U;
        end
        TERM = TERM * (P-0.5)/P;
        TQ = TQ*QSQ;
        TQSUM = TQSUM + TQ;
        TOLD = T;
        TTERM = TERM/(2*P + 3);
        TQTERM = TTERM * TQSUM;
        T = T - U0I((1.5*P + 0.25)*TQTERM / (P*P - 0.25 - TTMOLD*TQ));
        TTMOLD = TTERM;
        TQTERM = TQTERM*P;
        if (L1)
            DT = DT+ TQTERM*U1I;
        end
        if (L2)
            D2T = D2T+ TQTERM*U2I*(P-1);
        end
        if (L3)
            D3T = D3T+ TQTERM*U3I*(P-1)*(P-2);
        end
    end
    if (L3)
        D3T = 8*X*(1.5*D2T - XSQ * D3T);
    end
    if (L2)
        D2T = 2*(2*XSQ*D2T - DT);
    end
    if (L1)
        DT = -2D0*X*DT;
    end
    T = T/XSQ;
end
end
