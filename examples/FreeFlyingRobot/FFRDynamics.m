function output = FFRDynamics(t,x,u,auxdata)

alpha = auxdata.alpha;
beta  = auxdata.beta;

y1    = x(:,1);
y2    = x(:,2);
y3    = x(:,3);
y4    = x(:,4);
y5    = x(:,5);
y6    = x(:,6);
u1    = u(:,1);
u2    = u(:,2);
u3    = u(:,3);
u4    = u(:,4);
T1    = u1-u2;
T2    = u3-u4;

y1dot   = y4;
y2dot   = y5;
y3dot   = y6;
y4dot   = (T1+T2).*cos(y3);
y5dot   = (T1+T2).*sin(y3);
y6dot   = alpha*T1-beta*T2;
output  = [y1dot, y2dot, y3dot, y4dot, y5dot, y6dot];
end
