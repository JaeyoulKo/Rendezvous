clear all; close all; clc;
finp = fopen('IRIDIUM33_TLE_171201.txt','r');
mu = 398600.442;
Q = [];

tline = fgetl(finp);
while(tline~=-1)
    if(tline(12)=='D')
        tline = fgetl(finp);
        a = str2num(tline(54:59));
        b = str2num(tline(60:61));
        Bstar = a*1e-5*10^b;
        
        tline = fgetl(finp);
        a = str2num(tline(1:26));
        b = 1e-7*str2num(tline(27:33));
        c = str2num(tline(34:69));

        ID = a(2);
        i = a(3);
        OMEGA = a(4);
        e = b(1);
        w = c(1);
        M = c(2);
        N = c(3);

        N = N*2*pi/24/3600;
        a = (mu/N/N)^(1/3);

        Q = [Q; ID,a,e,i,OMEGA,w,M,N,Bstar];
    else
        tline = fgetl(finp);
        tline = fgetl(finp);    
    end
    tline = fgetl(finp);    
end