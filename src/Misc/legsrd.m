function [varargout]=legsrd(n)

% x=legsrd(n) returns n Legendre-Gauss-Radau points with x(1)=-1.
% [x,w]= legsrd(n) returns n Legendre-Gauss-Radau points and weights
% Eigenmethod is used for computing nodes.
% Last modified on August 30, 2011

% indices
j=0:n-2;
% Main diagonal
A=diag(1./((2*j+1).*(2*j+3)));   
j=1:n-2;
% Create Jacobi matrix
A=A+diag(sqrt(j.*(j+1))./(2*j+1),1) ...
    +diag(sqrt(j.*(j+1))./(2*j+1),-1);
% Compute eigenvalues     
x= sort(eig(sparse(A)));              
x=[-1;x];
varargout{1}=x;
if nargout==1
    return; 
end
y=lepoly(n-1,x); 
% Return the weights
varargout{2}= (1-x)./(n^2*y.^2); 
end


 






