function [AB R fgmask] = combinesynth(AB1,AB2,R1,R2,fgmask1,fgmask2)

if ~isequal(size(AB1),size(AB2),size(R1),size(R2))
    disp('Shit ain''t the same size');
    return
end

if nargin < 6
    fgmask2 = ones(size(AB2));
end
if nargin < 5
    fgmask1 = ones(size(AB1));
end

[M N] = size(AB1);

for m = 1:M
    for n = 1:N
        if (~fgmask1(m,n) && ~fgmask2(m,n))
            fgmask(m,n) = 0;
            AB(m,n) = max(AB1(m,n),AB2(m,n));
            R(m,n) = min(R1(m,n),R2(m,n));
        elseif ~fgmask1(m,n)
            fgmask(m,n) = fgmask2(m,n);
            AB(m,n) = AB2(m,n);
            R(m,n) = R2(m,n);
        elseif ~fgmask2(m,n)
            fgmask(m,n) = fgmask1(m,n);
            AB(m,n) = AB1(m,n);
            R(m,n) = R1(m,n);
        elseif R1(m,n) <= R2(m,n)
            fgmask(m,n) = fgmask1(m,n);
            AB(m,n) = AB1(m,n);
            R(m,n) = R1(m,n);
        else
            fgmask(m,n) = fgmask2(m,n);
            AB(m,n) = AB2(m,n);
            R(m,n) = R2(m,n);
        end
    end
end
            