clear all
clc
addpath('./../Lambert_Gooding/')
addpath('./../CoordinateTransformation/')

%Draw Chen et al. JGCD & Jun Bang 1st case Contour Plot 
% -> Jun Bang 논문 Case Study 1 부분의 Parameter는 typo 인 것 같다.
startOrbit = [6371+2000,0.002,60*pi/180,30*pi/180,0*pi/180,0*pi/180];
targetOrbit = [6371+36000,0.0002,55*pi/180,35*pi/180,-20*pi/180,30*pi/180];
Z = zeros(500,500);
i=1;
for td = linspace(20000,30000,500)
    j=1;
    for ttf = linspace(40000,70000,500)
        Z(j,i)=DelV_LambRend_tdttf(startOrbit,targetOrbit,[td;ttf]);
        j=j+1;
    end
    i=i+1;
end
contour(linspace(20000,30000,500),linspace(40000,70000,500),Z,50);

% %% Draw Jun Bang two-phase 논문 설명 contour plot
% startOrbit = [6678,0,0,0,0,0];
% targetOrbit = [13878,0.02,30*pi/180,20*pi/180,0,90*pi/180];
% Z = zeros(6000,6000);
% i=1;
% for td = linspace(300,20000,6000)
%     j=1;
%     for ttf = linspace(3000,20000,6000)
%         Z(j,i)=DelV_LambRend_tdttf(startOrbit,targetOrbit,[td;ttf]);
%         j=j+1;
%     end
%     i=i+1;
% end
% contour(linspace(300,20000,6000),linspace(3000,20000,6000),Z,50);

% startOrbit = [7081.884743,0.0012296,86.3839*pi/180,204.5888*pi/180,239.1039*pi/180,179.7762*pi/180];
% targetOrbit = [7158.977181,0.0016367,86.3876*pi/180,217.3307*pi/180,92.2436*pi/180,330.3728*pi/180];
% [a,e,i,o,w,nu] = ijk2keplerian([7241 4181 0],[-1.728 2.993 5.985]);
% startOrbit = [a,e,i*pi/180,o*pi/180,w*pi/180,nu*pi/180];
% [a,e,i,o,w,nu] = ijk2keplerian([31760 27391 6027],[-1.430 1.114 2.475]);
% targetOrbit = [a,e,i*pi/180,o*pi/180,w*pi/180,nu*pi/180];

% figure()
