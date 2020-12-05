function [X,T] = tbp_per2 (x, tf, delta)

%Simulation of perturbed two-body problem (TBP) under
%J2 (zonal harmonics) using approximated secular motion

%Inputs:
%         x: State vector [km;km/s]
%         tf: Final simulation time [s]
%         delta: Simulation time step [s]

%Outputs:
%         X: State history matrix [km;km/s]
%         T: Time values vector [s]

mu=3.986004418e+5; %Standard gravitational parameter [km^3/s^2]
J2=1.082635854e-3; %J2 zonal coefficient [dimensionless]
Re=6378.1363; %Earth radius [km]

oe=rv2oe(x(1:3),x(4:6)); %Obtaining the initial orbital elements

a0=oe(1); %Semimajor axis [km]
e0=oe(2); %Eccentricity [dimensionless]
inc0=oe(3); %Inclination [rad]
W0=oe(4); %RAAN [rad]
w0=oe(5); %Argument of periapsis [rad]
theta0=oe(6); %True anomaly [rad]
E0=2*atan2((1-e0^2)^0.5*sin(theta0)/(1+e0*cos(theta0)),...
           (e0+cos(theta0))/(1+e0*cos(theta0))); %Initial eccentric anomaly [rad]
M0=E0-e0*sin(E0); %Initial mean anomaly [rad]

i=1; %Index number

X=zeros(6,floor(tf/delta)+1); %Initiation of state history matrix [km;km/s]
T=zeros(1,floor(tf/delta)+1); %Initiation of time vector [s]

E=0; %An initial guess for Kepler probelm used only for the first loop [rad]
epsilon=1e-6; %Stopping criterion for solving the Kepler problem
N=1000; %Maximum iteration number for Kepler probelm

for t=0:delta:tf
    
    [r,v]=oe2rv(oe); %Orbital element to r and v [km,km/s]
    
    x=[r;v]; %State vector [km;km/s]
    
    X(:,i)=x; %Saving history matrix [km;km/s]
    
    T(i)=t; %Saving time values [s]
    
    a=a0; %Semimajor axis does not change [km]
    
    e=e0; %Eccentricity does not change [dimensionless]
    
    n=(mu/a^3)^0.5; %Mean anomaly (or mean angular motion) [1/s]
    
    p=a*(1-e^2); %Semi-latus rectum [km]
    
    inc=inc0; %Inclination does not change [rad]
    
    dW=-1.5*J2*Re^2*n/p^2*cos(inc); %Time rate of RAAN [rad/s]
    W=W0+dW*t; %Calculating RAAN at each time [rad]
    
    dw=-1.5*J2*Re^2*n/p^2*(2.5*sin(inc)^2-2); %Time rate of argument of periapsis [rad/s]
    w=w0+dw*t; %Calculating argument of periapsis at each time [rad]
    
    dM=n*(1-1.5*J2*Re^2/p^2*(1-e^2)^0.5*(1.5*sin(inc)^2-1)); %Time rate of Mean motion [rad/s]
    M=M0+dM*t; %Calculating mean motion at each time [rad]
    
    [theta,E]=Kepler(M,e,epsilon,N,E); %Solving Kepler problem [rad,rad]
    
    oe=[a;e;inc;W;w;theta]; %Vector of orbital elements
    
    i=i+1; %Index updating
    
end

end