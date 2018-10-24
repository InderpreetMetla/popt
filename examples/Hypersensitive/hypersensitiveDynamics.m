function output = hypersensitiveDynamics(t,x,u,auxdata)
output = -x.^3+u;
end