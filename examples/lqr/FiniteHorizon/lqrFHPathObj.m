function output = lqrFHPathObj(t,x,u,auxdata)
output = 0.5*(x(:,1)+u(:,1).^2);
end