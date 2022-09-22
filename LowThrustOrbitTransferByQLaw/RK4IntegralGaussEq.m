function newState = RK4IntegralGaussEq(newState,input,dt, Parameters)
mu=Parameters(1);
ThrustMax=Parameters(2);
c=Parameters(3);
k1 = TwoBodyVOPDynamicsEquinoctial(newState, input, [mu;ThrustMax;c]);
k2 = TwoBodyVOPDynamicsEquinoctial(newState + k1*0.5 * dt, input, [mu;ThrustMax;c]);
k3 = TwoBodyVOPDynamicsEquinoctial(newState + k2*0.5 * dt, input, [mu;ThrustMax;c]);
k4 = TwoBodyVOPDynamicsEquinoctial(newState + k3 * dt, input, [mu;ThrustMax;c]);
newState = newState + dt*(k1 + 2*(k2 + k3) + k4)/6;
newState = real(newState);
end