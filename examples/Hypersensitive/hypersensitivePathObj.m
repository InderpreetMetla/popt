function output = hypersensitivePathObj(t,x,u,auxdata)
output = 0.5*(x.^2+u.^2);
end