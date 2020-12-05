function xnew=rk4(x,delta,p)

%Runge-Kutta integration

% Input(s):
%           x:     State vector at t
%           delta: Integration time step--assumes delta=1e-3 if not specified
%           p:     Vector of additional parameters--assumes not specified if not specified
%
% Output(s):
%           xnew:  State vector at t+delta.

switch nargin
    case 1 %If only the first argument is specified
        delta=1e-3;
        xd1=dyn(x,p);
        xd2=dyn(x+delta/2*xd1);
        xd3=dyn(x+delta/2*xd2);
        xd4=dyn(x+delta*xd3);
        xnew=x+delta/6*(xd1+2*xd2+2*xd3+xd4);
    case 2 %If only the first two arguments are specified
        xd1=dyn(x,p);
        xd2=dyn(x+delta/2*xd1);
        xd3=dyn(x+delta/2*xd2);
        xd4=dyn(x+delta*xd3);
        xnew=x+delta/6*(xd1+2*xd2+2*xd3+xd4);
    otherwise
        xd1=dyn(x,p);
        xd2=dyn(x+delta/2*xd1,p);
        xd3=dyn(x+delta/2*xd2,p);
        xd4=dyn(x+delta*xd3,p);
        xnew=x+delta/6*(xd1+2*xd2+2*xd3+xd4);
end

end