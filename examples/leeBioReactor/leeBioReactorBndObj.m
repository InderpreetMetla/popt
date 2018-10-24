function output = leeBioReactorBndObj(t0,tf,x0,xf,auxdata)
x1f = xf(1); x4f = xf(4);
output = -x1f.*x4f;
end