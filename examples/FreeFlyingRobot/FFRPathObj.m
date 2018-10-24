function output = FFRPathObj(t,x,u,auxdata)
gamma = auxdata.gamma;
u1    = u(:,1);
u2    = u(:,2);
u3    = u(:,3);
u4    = u(:,4);

output = (u1+u2+u3+u4).*gamma;
end


