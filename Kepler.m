function [theta,E] = Kepler (M,e,epsilon,N,E0)

%Solving Kepler problem using Newton-Raphson method 

%Inputs:
%         M: Mean anomaly [rad]
%         e: Eccentricity [dimensionless]--assumes e=0 if not specified
%         epsilon: Stopping error bound--assumes epsilon=1e-6 if not specified
%         N: Maximum number of iterations--assumes N=1e+2 if not specified
%         E0: An initial guess for the eccentric anomaly [rad]--assumes E0=0 if not specified

%Outputs:
%         E: Eccentric anomaly [rad]
%         theta: True anomaly [rad]

switch nargin
    case 1 %If only the first argument is specified
        e=0;
        epsilon=1e-6;
        N=1e+2;
        E0=0;
    case 2 %If only the first two arguments are specified
        epsilon=1e-6;
        N=1e+2;
        E0=0;
    case 3 %If only the first three arguments are specified
        N=1e+2;
        E0=0;
    case 4 %If only the first four arguments are specified
        E0=0;
end

E=E0; %Initiation of the algorithm

for i=1:N
    
    E=E-(M-E+e*sin(E))/(e*cos(E)-1); %Updating according to the Newton-Raphson method  
    
    if abs(M-E+e*sin(E))<epsilon
        
        break %Breaking if the error bound is satisfied as a stopping criterion
        
    end
    
end

sintheta=sin(E)*(1-e^2)^0.5/(1-e*cos(E)); %Sine of true anomaly
costheta=(cos(E)-e)/(1-e*cos(E)); %Cosine of true anomaly

theta=atan2(sintheta,costheta); %Obtaining true anomaly by tangent inverse

end

