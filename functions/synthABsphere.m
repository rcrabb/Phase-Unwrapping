function [AB R fgmask] = synthABsphere(BrightnessConst,r,center,A,albedofg, bgdist, albedobg)

center = center + .00001;

if nargin < 7
    albedobg = .6;
end
if nargin < 6
    bgdist = center(3) + r*2;
end

if nargin < 5
    albedofg = .8;
end
if nargin < 4
    A = gethlut();
end
A(A==0) = -.00001;

[M N] = size(BrightnessConst);

[ABback] = synthAB(bgdist*ones([M N]),albedobg*ones([M N]),BrightnessConst,A,'hack');
Rback = bgdist./A(:,:,3);

% Solve quadratic equation to find intersection of sphere and pixel vector a

a = sum(A.^2,3);
b = -2*(dot(A,repmat(reshape(center,[1 1 3]),[M N 1]),3));
c = (sum(center.^2)-r^2)*ones(size(a));

d = (-b - real(sqrt(b.^2 - 4*a.*c))) ./ (2*a);
mask = b.^2 >= 4*a.*c;

xyz = repmat(d,[1 1 3]).*A;

normals = (xyz - repmat(reshape(center,[1 1 3]),[M N 1]))/r;

adotn = abs(dot(A,normals,3));

AB = BrightnessConst .* albedofg .* adotn ./ d.^2;
AB(~mask) = ABback(~mask);
R = d;
R(~mask) = Rback(~mask);
fgmask = mask;
return
%%
r = .25;
center = [-.2,.3,1.500];
bgdist = 2;
%%
r = 500;
center = [0 0 1500];
bgdist = 2;