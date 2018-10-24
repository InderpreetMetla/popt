function output = leeBioReactorDynamics(t,x,u,auxdata)
s = x;

x1 = s(:,1);
x2 = s(:,2);
x3 = s(:,3);
x4 = s(:,4);
x5 = s(:,5);
x6 = s(:,6);
x7 = s(:,7);
u1 = s(:,8);
u2 = s(:,9);

k1 = 0.09*x5./(0.034+x5);
g1 = (x3./(14.35+x3.*(1+x3./111.5))).*(x6+(0.22*x7)./(0.22+x5));
Rfp = (0.233*x3./(14.35+x3.*(1+x3./111.5))).*((0.0005+x5)./(0.022+x5));

x1dot = u1+u2;
x2dot = g1.*x2-(u1+u2).*x2./x1;
x3dot = 100*u1./x1-(u1+u2).*x3./x1-g1.*x2/0.51;
x4dot = Rfp.*x2-(u1+u2).*x4./x1;
x5dot = 4*u2./x1-(u1+u2).*x5./x1;
x6dot = -k1.*x6;
x7dot = k1.*(1-x7);

x8dot = u(:,1);
x9dot = u(:,2);

output = [x1dot x2dot x3dot x4dot x5dot x6dot x7dot x8dot x9dot];
end