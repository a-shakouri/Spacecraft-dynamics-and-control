function [dx] = dyn (x, p)

%Two-body problem (TBP) unactuated state-space model

%Inputs:
%          x: State vector [km;km/s]
%          p: Switching between perturbed and unperturbed models:--assumes p=0 if not specified
%             p=0: Unperturbed TBP
%             p=1: Perturbed TBP under J2

%Outputs:
%         dx: Time derivative of state vector x [km/s;km/s^2]

mu=3.986004418e+5; %Standard gravitational parameter [km^3/s^2]
J2=1.082635854e-3; %J2 zonal coefficient [dimensionless]
Re=6378.1363; %Earth radius [km]

dx=zeros(6,1); %Initiation of state derivative vector [km/s;km/s^2]

r=x(1:3); %Position vector [km]

dx(1:3)=x(4:6); %Velocity vector [km/s]

switch nargin
    case 1 %If only the first argument is specified
        p=0;
end

if p==0
    
    dx(4:6)=-mu/norm(r)^3*r; %Unperturbed acceleration [km/s^2]
    
elseif p==1
    
    rj2=[(1+1.5*J2*(Re/norm(r))^2*(1-5*(r(3)/norm(r))^2))*r(1)
         (1+1.5*J2*(Re/norm(r))^2*(1-5*(r(3)/norm(r))^2))*r(2)
         (1+1.5*J2*(Re/norm(r))^2*(3-5*(r(3)/norm(r))^2))*r(3)]; %Modified position vector [km]
    
    dx(4:6)=-mu/norm(r)^3*rj2; %Perturbed acceleration [km/s^2]
    
else
    
    error('Parameter "p" should be 0 (for unperturbed case) or 1 (for perturbed case).')
    
end

end
