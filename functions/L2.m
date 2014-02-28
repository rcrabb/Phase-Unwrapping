%L2 - L2 Norm
% Compute the L2 norm
%
% [X2] = L2( X, SIGMA, TRUNCATE )
% 
% Input:
%   X - vector of which the norm is determined
%   SIGMA - divides X
%   TRUNCATE - maximum value 
function [x2] = L2(x,sigma,truncate)


if (nargin < 2 || truncate == 0)
     x2 = x.^2/(2*(sigma.^2));
else
    x2 = min(x.^2/(2*(sigma.^2)),truncate^2*(.5*sigma^-2));
end

end