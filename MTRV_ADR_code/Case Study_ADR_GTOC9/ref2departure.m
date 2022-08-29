load 'ref_element.mat'

t_d = 23467; % departure epoch
mu = 398600.4418;
re = 6378.137;
d2r = pi/180;
J2term = -1.5*1.08261e-3*sqrt(mu)*re^2*24*3600/d2r; % deg/day
Q = zeros(123,6);

for i=1:123
    t_ref = ref_element(i,1);  
    a_ref = ref_element(i,2); 
    e_ref = ref_element(i,3); 
    i_ref = ref_element(i,4); 
    om_ref = ref_element(i,5); 
    w_ref = ref_element(i,6); 
    M_ref = ref_element(i,7);

    omdot = J2term*cos(i_ref*d2r)/a_ref^3.5/(1-e_ref^2)^2;
    wdot = -0.5*J2term*(5*cos(i_ref*d2r)^2-1)/a_ref^3.5/(1-e_ref^2)^2;
    om_d = mod(om_ref+omdot*(t_d-t_ref),360);
    w_d = mod(w_ref+wdot*(t_d-t_ref),360);
    M_d = mod(M_ref+sqrt(mu/a_ref^3)*(t_d-t_ref)*86400/d2r,360);

    Q(i,:) = [a_ref,e_ref,i_ref,om_d,w_d,M_d];
end

save('target_list.mat','Q')