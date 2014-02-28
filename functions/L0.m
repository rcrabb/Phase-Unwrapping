%L0 - L0 Norm
% Compute the L0 norm
%
% [X2] = L0( X, SIGMA, TRUNCATE )
% 
% Input:
%   X - vector of which the norm is determined
%   SIGMA - divides X
%   TRUNCATE - maximum value 
function [xOut] = L0(x,sigma,truncate)

if (nargin < 2 || truncate == 0)
    xOut = x==0;
else
    xOut = abs(x) > truncate;
end

end