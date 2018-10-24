function output = leeBioReactorPathObj(t,x,u,auxdata)
N = length(t);
rho = 0.1/N;
u1dot = u(:,1); u2dot = u(:,2); 
output = rho.*((u1dot).^2+(u2dot).^2);
end