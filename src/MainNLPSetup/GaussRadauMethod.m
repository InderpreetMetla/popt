function Gauss = GaussRadauMethod(n,breakpts,nInts)
%******************************************************************%
% This function computes the Gauss Radau Points, Gauss Radau Weights, 
% and the Gauss Radau Differentiation/Integration Matrices to use 
% in the Gauss Radau Pseudospectral Method.
%
% The Points and Weights are calculated by legsrd.m which is written
% by Professor Wang Li-Lian and available from the following ref.:
% "Shen, J., Wang, L. and Tang, T. (2011). 'Spectral Methods: 
% Algorithms, Analysis and Applications'. Springer Series in 
% Computational Mathematics." The MATLAB codes presented in the book
% can be found at <http://www.ntu.edu.sg/home/lilian/book.htm>.
%
% The Differentation Matrix is found using poldif.m which is obtained
% from the ref.:
% "Weideman, J.A.C. and Reddy, S.C. (2000). 'A MATLAB Differentiation 
% Matrix Suite'. ACM Transactions on Mathematical Software, 
% Vol. 26, No. 4". Available at: 
% <http://dip.sun.ac.za/~weideman/research/differ.html>. 
% 
% Inputs:
% - n        :     Number of nodes desired by user in a phase
% - breakpts :     Location of interval break points in a phase
% - nInts    :     Number of intervals in a phase
%
% Outputs:
%  - Gauss:     struct with following fields
%    - Gauss.Points     : Collocation Points
%    - Gauss.Weights    : Quadrature Weights as column vector
%    - Gauss.DM         : GRPM Differentiation Matrix
%    - Gauss.IM         : GRPM Integration Matrix - used for 
%                         Error calcs.
%    - Gauss.DM1        : Sparsity of GPM Differential Matrix
%    - Gauss.MM         : Multiplier Matrix for the Error Calcs
%                         This is just used to vectorize the calc
%                         of the Y0hat term in the Absolute Error.
%                         See computeError.m for more. 
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%******************************************************************%

rowIdxD = 0;colIdxD = 0;rowIdxI = 0;colIdxI = 0;
GP = []; GW = [];
for iInt = 1:nInts
    N = n(iInt);
    a = breakpts(iInt); b = breakpts(iInt+1);
    % Use legsrd to compute points and weights
    [Gauss.Points, Gauss.Weights] = legsrd(N);
    % Scale the Points and Weights according to a and b
    Points  = 0.5*(b - a)*(Gauss.Points + 1) + a;
    Weights = 0.5*(b - a)*Gauss.Weights;
    % Use poldiff to calc derivative of lagrange polynomials
    DM = poldif([Points;b],1);
    DM = DM(1:end-1,:);
    % Integration matrix is just the inverse of the differential matrix
    % from the 2nd column onwards
    IM = inv(DM(:,2:end));
    MultiplierMatrix = zeros(size(DM));
    MultiplierMatrix(:,1) = 1;
    [nRowD,nColD] = size(DM);
    [nRowI,nColI] = size(IM);
    % Store the matrices for each interval
    DI.DiffMatrix(1+rowIdxD:rowIdxD+nRowD,1+colIdxD:...
        colIdxD+nColD) = DM;
    DI.IntMatrix(1+rowIdxI:rowIdxI+nRowI,1+colIdxI:...
        colIdxI+nColI) = IM;
    DI.X0Multiplier(1+rowIdxD:rowIdxD+nRowD,1+colIdxD:...
        colIdxD+nColD) = MultiplierMatrix;
    rowIdxD = rowIdxD + nRowD;
    colIdxD = colIdxD + nColD - 1;
    rowIdxI = rowIdxI + nRowI;
    colIdxI = colIdxI + nColI;
    % Append all points and weights
    GP = [GP;Points];
    GW = [GW;Weights];
end
Gauss.Points  = GP;
Gauss.Weights = GW.';
% Remove first row as not included in collocation
Gauss.DM  = sparse(DI.DiffMatrix);
Gauss.IM  = sparse(DI.IntMatrix);
Gauss.DM1 = spones(DI.DiffMatrix);
Gauss.MM  = sparse(DI.X0Multiplier);
end


