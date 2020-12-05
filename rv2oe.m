function [oe] = rv2oe(r,v)

%Obtaining the classical orbital elements from inertial position and velocity vectors

%Inputs:
%         r: Position vector [km]
%         v: Velocity vector [km/s]

%Outputs:
%         oe: Vector of orbital elements
%             oe(1): Semimajor axis [km]
%             oe(2): Eccentricity [dimensionless]
%             oe(3): Inclination [rad]
%             oe(4): RAAN [rad]
%             oe(5): Argument of periapsis [rad]
%             oe(6): True anomaly [rad]

mu=3.986004418e+5; %Standard gravitational parameter [km^3/s^2]

h=cross(r,v); %Specific angular momentum [km^2/s]
inc=acos(h(3)/norm(h)); %Inclination [rad]
 
n=[-h(2);h(1);0]; %Vector perpendicular to h and Earth's north pole

if n(2)>=0
    
    W=acos(n(1)/norm(n)); %RAAN [rad]
    
elseif n(2)<0
    
    W=2*pi-acos(n(1)/norm(n)); %RAAN [rad]
    
end

E=(1/mu)*(cross(v,h))-r/norm(r); %Eccentricity vector [dimensionless]

if E(3)>=0
    
    w=acos((E'*n)/(norm(n)*norm(E))); %Argument of periapsis [rad]

elseif E(3)<0
    
    w=2*pi-acos((E'*n)/(norm(n)*norm(E))); %Argument of periapsis [rad]

end

Vr=v'*r/norm(r); %Derivative of position magnitude [km/s]

if Vr>=0
    
    theta=acos((E'*r)/(norm(r)*norm(E))); %True anomaly [rad]
    
elseif Vr<0
    
    theta=2*pi-acos((E'*r)/(norm(r)*norm(E))); %True anomaly [rad]

end

a=(norm(h)^2/mu)/(1-norm(E)^2); %Semimajor axis [km]
e=norm(E); %Eccentricity [dimensionless]

oe=[a,e,inc,W,w,theta]; %Vector of orbital elements

end



