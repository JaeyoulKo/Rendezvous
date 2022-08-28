function [N, VR11, VT11, VR12, VT12, VR21, VT21, VR22, VT22] = VLAMB(GM, R1, R2, TH, TDELT)

TWOPI=2*pi;
M = TH/TWOPI;
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
[N,X1,X2] = XLMAB(M,Q,QSQFM1,T);
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
