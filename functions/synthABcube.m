function [AB fgmask] = synthABcube(BrightnessConst,sidelength,center,A,albedofg, bgdist, albedobg)

center = center + .00001;
s = sidelength/2;

if nargin < 7
    albedobg = .6;
end
if nargin < 6
    bgdist = center(3) + sidelength*3;
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

planedist = (center(3)-s)./A(:,:,3);
mask = planedist.*A(:,:,1) >= center(1) - s & planedist.*A(:,:,1) <= center(1) + s &  planedist.*A(:,:,2) >= center(2) - s & planedist.*A(:,:,2) <= center(2) + s;
planedist(~mask) = Inf;
d(:,:,1) = planedist;

planedist =(center(2)-s)./A(:,:,2);
mask = planedist.*A(:,:,1) >= center(1) - s & planedist.*A(:,:,1) <= center(1) + s &  planedist.*A(:,:,3) >= center(3) - s & planedist.*A(:,:,3) <= center(3) + s;
planedist(~mask) = Inf;
d(:,:,2) = planedist;

planedist = (center(2)+s)./A(:,:,2);
mask = planedist.*A(:,:,1) >= center(1) - s & planedist.*A(:,:,1) <= center(1) + s &  planedist.*A(:,:,3) >= center(3) - s & planedist.*A(:,:,3) <= center(3) + s;
planedist(~mask) = Inf;
d(:,:,3) = planedist;

planedist =(center(1)-s)./A(:,:,1);
mask = planedist.*A(:,:,2) >= center(2) - s & planedist.*A(:,:,2) <= center(2) + s & planedist.*A(:,:,3) >= center(3) - s & planedist.*A(:,:,3) <= center(3) + s;
planedist(~mask) = Inf;
d(:,:,4) = planedist;

planedist = (center(1)+s)./A(:,:,1);
mask = planedist.*A(:,:,2) >= center(2) - s & planedist.*A(:,:,2) <= center(2) + s & planedist.*A(:,:,3) >= center(3) - s & planedist.*A(:,:,3) <= center(3) + s;
planedist(~mask) = Inf;
d(:,:,5) = planedist;
%d(d<0) = Inf;
[d dind] = min(d,[],3);


adotn(:,:,1) = A(:,:,3);
adotn(:,:,2) = A(:,:,2);
adotn(:,:,3) = A(:,:,2);
adotn(:,:,4) = A(:,:,1);
adotn(:,:,5) = A(:,:,1);
for m = 1:M
    for n = 1:N
        adotn(m,n,1) = adotn(m,n,dind(m,n));
    end
end
adotn = abs(adotn(:,:,1));
% adotn = [A(:,:,3) A(:,:,2) A(:,:,2) A(:,:,1) A(:,:,1)];
% normals = repmat(reshape([0 0 1;0 1 0; 0 1 0; 1 0 0;1 0 0], [1 1 5 3]),[M N 1 1]);
% adotn = abs(dot(A,normals,3));

mask = d.*A(:,:,1) >= center(1) - s & d.*A(:,:,1) <= center(1) + s;
mask = mask & d.*A(:,:,2) >= center(2) - s & d.*A(:,:,2) <= center(2) + s;
mask = mask & d.*A(:,:,3) >= center(3) - s & d.*A(:,:,3) <= center(3) + s;

AB = BrightnessConst .* albedofg .* adotn ./ d.^2;
AB(~mask) = ABback(~mask);
fgmask = mask;
return
%%
sidelength = .25;
center = [-.2,.3,1.500];
bgdist = 2;
%%
r = .25
center = [0 0 1.500];
bgdist = 2;