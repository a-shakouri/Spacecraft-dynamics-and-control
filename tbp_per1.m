function [X,T] = tbp_per1 (x, tf, delta)

%Simulation of perturbed two-body problem (TBP) under J2 (zonal sphere harmonics)

%Inputs:
%         x: State vector [km;km/s]
%         tf: Final simulation time [s]
%         delta: Simulation time step [s]

%Outputs:
%         X: State history matrix [km;km/s]
%         T: Time values vector [s]

i=1; %Index number

X=zeros(6,floor(tf/delta)+1); %Initiation of state history matrix [km;km/s]
T=zeros(1,floor(tf/delta)+1); %Initiation of time vector [s]

for t=0:delta:tf
    
    X(:,i)=x; %Saving history matrix [km;km/s]
    
    T(i)=t; %Saving time values [s]
    
    x=rk4(x,delta,1); %4th order Runge-Kutta integration and updating [km;km/s]
    
    i=i+1; %Index updating

end

end

