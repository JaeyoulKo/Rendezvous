function [nSol, v1vec, v2vec] = Lambert(Mu, r1vec, r2vec, delT)
nSol=0;
r1 = norm(r1vec);
r2 = norm(r2vec);
ur1 = r1vec/r1;
ur2 = r2vec/r2;
ut1 = cross(cross(r1vec,r2vec),ur1);
ut1 = ut1/norm(ut1);
ut2 = cross(cross(r1vec,r2vec),ur2);
ut2 = ut2/norm(ut2);
theta = acos(dot(ur1,ur2));

preAlocSize=100;
v1vec = zeros(3,preAlocSize); % return vector preallocation
v2vec = zeros(3,preAlocSize);

numIter=0;

while(1)
    numIter = numIter + 1;
    [N, Vr11, Vt11, Vr12, Vt12, Vr21, Vt21, Vr22, Vt22] = ...
        VLAMB(Mu,r1,r2,theta,delT);
    if (N==0)
        theta=theta-(numIter-1)*2*pi;
        break;
    end
    if (N==1)
        nSol=nSol+1;
        v1 = Vr11*ur1 + Vt11 * ut1;
        v2 = Vr12*ur2 + Vt12 * ut2;
        v1vec(1:3,nSol)=v1;
        v2vec(1:3,nSol)=v2;
    end
    if (N==2)
        nSol=nSol+1;
        v1 = Vr11*ur1 + Vt11 * ut1;
        v2 = Vr12*ur2 + Vt12 * ut2;
        v1vec(1:3,nSol)=v1;
        v2vec(1:3,nSol)=v2;

        nSol=nSol+1;
        v1 = Vr21*ur1 + Vt21 * ut1;
        v2 = Vr22*ur2 + Vt22 * ut2;
        v1vec(1:3,nSol)=v1;
        v2vec(1:3,nSol)=v2;
    end
    theta=theta+2*pi;
end
theta = 2*pi-theta; %retrograde direction
while(1)
    numIter = numIter + 1;
    [N, Vr11, Vt11, Vr12, Vt12, Vr21, Vt21, Vr22, Vt22] = ...
        VLAMB(Mu,r1,r2,theta,delT);
    if (N==0)
        break;
    end
    if (N==1)
        nSol=nSol+1;
        v1 = Vr11*ur1 + Vt11 * (-ut1);
        v2 = Vr12*ur2 + Vt12 * (-ut2);
        v1vec(1:3,nSol)=v1;
        v2vec(1:3,nSol)=v2;
    end
    if (N==2)
        nSol=nSol+1;
        v1 = Vr11*ur1 + Vt11 * (-ut1);
        v2 = Vr12*ur2 + Vt12 * (-ut2);
        v1vec(1:3,nSol)=v1;
        v2vec(1:3,nSol)=v2;

        nSol=nSol+1;
        v1 = Vr21*ur1 + Vt21 * (-ut1);
        v2 = Vr22*ur2 + Vt22 * (-ut2);
        v1vec(1:3,nSol)=v1;
        v2vec(1:3,nSol)=v2;
    end
    theta=theta+2*pi;
    v1vec = v1vec(:,1:nSol);
    v2vec = v2vec(:,1:nSol);
end
end