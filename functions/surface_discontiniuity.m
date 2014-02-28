%% Examine scene to find a correlation between dAB and dZ
function [Bhist Rhist] = surface_discontiniuity()
cvprdata = 'C:\data\cvprdata\';
dirlist = getdir(cvprdata);

Bhist = [];
Rhist = [];

for testnum = 1:length(dirlist)
dataname = dirlist{testnum};
[R AB] = pgm2mat([cvprdata dataname]);

R = median(R,3);
AB = median(AB,3);

% find the deltas
dvAB = AB(1:end-1,:)-AB(2:end,:);
dhAB = AB(:,1:end-1)-AB(:,2:end);
dvR = R(:,1:end-1)-R(:,2:end);
dhR = R(1:end-1,:)-R(2:end,:);

% concatenate and remove too dim of pixels
vB = AB(1:end-1,:); vB(:,:,2) = AB(2:end,:);
vB = min(vB,[],3);
vR = R(1:end-1,:); vR(:,:,2) = R(2:end,:);
vR = min(vR,[],3);
hB = AB(:,1:end-1); hB(:,:,2) = AB(:,2:end);
hB = min(hB,[],3);
hR = R(:,1:end-1); hR(:,:,2) = R(:,2:end);
hR = min(hR,[],3);

dAB = [dvAB(:); dhAB(:)];
dR = [dvR(:); dhR(:)];
AB = [vB(:); hB(:)];
R = [vR(:); hR(:)];

dAB = dAB(AB>50);
dR = dR(AB>50);

Bhist = [Bhist; dAB];
Rhist = [Rhist; dR];

end