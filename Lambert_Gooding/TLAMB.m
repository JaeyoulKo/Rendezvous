function [T,DT,D2T,D3T] = TLAMB(M,Q,QSQFM1,X,N)
SW = 0.4;
LM1 = (N==-1);
L1 = (N>=1);
L2 = (N>=2);
L3 = (N==3);
QSQ = Q*Q;
XSQ = X*X;
U = (1-X)*(1+X);
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
