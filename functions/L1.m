%L1 - L1 Norm
% Compute the L1 norm
%
% [X2] = L1( X, SIGMA, TRUNCATE )
% 
% Input:
%   X - vector of which the norm is determined
%   SIGMA - divides X
%   TRUNCATE - maximum value 
function [xOut] = L1(x,sigma,truncate)

if (nargin < 2 || truncate == 0)
    xOut = abs(x/sigma);
else
    xOut = min(abs(x/sigma),truncate/sigma);
end

end