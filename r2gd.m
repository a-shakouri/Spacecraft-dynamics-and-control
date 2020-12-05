function [lambda,phigd] = r2gd (r,epsilon,N)

%Calculation of geodetic latitude and logitude from ECEF position vector

%Inputs:
%         r: Position vector [km]
%         epsilon: Stopping error bound--assumes epsilon=1e-6 if not specified
%         N: Maximum number of iterations--assumes N=1e+2 if not specified

%Outputs:
%         lambda: Longitude [rad]
%         phigd: Geodetic latitude [rad]

aE=6378.1370; %Earth equatorial radius [km]
bE=6356.7523; %Earth polar radius [km]

fE=(aE-bE)/aE; %Earth flattening parameter [dimensionless]

eE=(2*fE-fE^2)^0.5; %Eccentricity of the Earth ellisoid [dimensionless]

lambda=atan2(r(2),r(1)); %Longitude [rad]

rxy=norm(r(1:2)); %Norm of the first two elements of position vector [km]

phigd=atan2(r(3),rxy); %Initial guess for geodetic latitude [rad]

switch nargin
    case 1 %If only the first argument is specified
        epsilon=1e-6;
        N=1e+2;
    case 2 %If only the first two arguments are specified
        N=1e+2;
end

for i=1:N
    
    phigd0=phigd; %Saving the previous value
    
    Rphi=aE/(1-eE^2*sin(phigd)^2); %Magnitude of Earth ellipsoid normal vector intersecting by polar axis [km]
    
    phigd=atan2(r(3)+Rphi*eE^2*sin(phigd),rxy); %Updating phigd by a fixed-point interation method
    
    if abs(phigd0-phigd)<epsilon
        
        break %Breaking if the error bound is satisfied as a stopping criterion
        
    end
    
end

end

