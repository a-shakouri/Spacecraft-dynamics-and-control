function [r,v] = oe2rv(oe)

%Obtaining the inertial position and velocity vectors from classical orbital elements

%Inputs:
%         oe: Vector of orbital elements:
%             oe(1): Semimajor axis [km]
%             oe(2): Eccentricity [dimensionless]
%             oe(3): Inclination [rad]
%             oe(4): RAAN [rad]
%             oe(5): Argument of periapsis [rad]
%             oe(6): True anomaly [rad]

%Outputs:
%         r: Position vector [km]
%         v: Velocity vector [km/s]

a=oe(1); %Semimajor axis [km]
e=oe(2); %Eccentricity [dimensionless]
inc=oe(3); %Inclination [rad]
W=oe(4); %RAAN [rad]
w=oe(5); %Argument of periapsis [rad]
theta=oe(6); %True anomaly [rad]

mu=3.986004418e+5; %Gravitational parameter [km^3/s^2]

gamma=atan2(e*sin(theta),1+e*cos(theta));

rn=a*(1-e^2)/(1+e*cos(theta)); %Position magnitude [km]

vn=(2*mu/rn-mu/a)^0.5; %Velocity magnitude [km/s]

r=rn*[cos(theta+w)*cos(W)-sin(theta+w)*cos(inc)*sin(W)
      cos(theta+w)*sin(W)+sin(theta+w)*cos(inc)*cos(W)
      sin(theta+w)*sin(inc)]; %Position vector in ECI [km]
  
v=vn*[-sin(theta+w-gamma)*cos(W)-cos(theta+w-gamma)*cos(inc)*sin(W)
      -sin(theta+w-gamma)*sin(W)+cos(theta+w-gamma)*cos(inc)*cos(W)
       cos(theta+w-gamma)*sin(inc)]; %Velocity vector ECI [km/s]

end

