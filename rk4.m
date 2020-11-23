function xnew=rk4(x,delta,p)

% Input(s):
%           x:     State vector at $t$.
%           delta: Integration time step.
%           p:     Vector of additional parameters.
%
% Output(s):
%           xnew:  State vector at $t+delta$.



xd1=dyn(x,p);
xd2=dyn(x+delta/2*xd1,p);
xd3=dyn(x+delta/2*xd2,p);
xd4=dyn(x+delta*xd3,p);
xnew=x+delta/6*(xd1+2*xd2+2*xd3+xd4);

end