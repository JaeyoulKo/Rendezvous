% [N, VR11, VT11, VR12, VT12, VR21, VT21, VR22, VT22] = ...
%     VLAMB(398600.442,7090.592562,7156.532950,111.439811*pi/180,1800.001000)
% % Solution 1 Case
% [N, VR11, VT11, VR12, VT12, VR21, VT21, VR22, VT22] = ...
%     VLAMB(398600.442,7090.592562,7154.349354,33.858423*pi/180,5408.247437)
% [N, VR11, VT11, VR12, VT12, VR21, VT21, VR22, VT22] = ...
%     VLAMB(398600.442,7090.592562,7154.349354,33.858423*pi/180+2*pi,5408.247437)
% % 3 Solution Case (1+2)
% [N, VR11, VT11, VR12, VT12, VR21, VT21, VR22, VT22] = ...
%     VLAMB(398600.442,7090.592562,7154.349354,33.858423*pi/180+4*pi,5408.247437)
% % 0 Solution due to large m(revolution) Case
% 
% [nSol, v1vec, v2vec] = Lambert(398600.442,...
%     [-3172.952007;-1872.974947;6058.122093], [5655.161661;4217.272457;1204.044631], 1800.001000)
% % 1 Solution Case
% [nSol, delV1vec, DelV2vec] = DelV_LambRend(398600.442,...
%     [-3172.952007;-1872.974947;6058.122093], [5.934684;2.447099;3.864907], ...
%     [5655.161661;4217.272457;1204.044631], [0.728885;1.138719;-7.340878], 1800.001000)
% 
% %  7090.592562 	 7156.532950 	 111.439811 	  
fileID = fopen("lambert_delV_Results.txt",'r');
tic
while (~feof(fileID))
    fgets(fileID);
    temp=fscanf(fileID,'%f,%f,%f');
    r1=temp(1:3);
    r2=temp(4:6);
    v1=temp(7:9);
    v2=temp(10:12);
    tf=temp(13);
    delV=DelV_LambRend(398600.442,r1,v1,r2,v2,tf);
    fgets(fileID);
    temp=fgets(fileID);
    delVComp=temp(8:end-1);
    delVComp=str2double(delVComp);
    if(abs(delV - delVComp)>0.0001)
        printf("error occurs");
    end
end
toc
fgets(fileID);
