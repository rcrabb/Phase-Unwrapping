function [energy emap] = BPEnergy(phaseimg,sigma,measure,truncate,is8connect,probXs,nStates)
%(phaseimg,sigma,is8connect, sigmaDiag)

[Ny, Nx] = size(phaseimg);
vertshift = phaseimg(2:end,:)-phaseimg(1:end-1,:);
horshift = phaseimg(:,2:end)-phaseimg(:,1:end-1);
diagshiftUpLeft = phaseimg(2:end,2:end)-phaseimg(1:end-1,1:end-1);
diagshiftDownLeft = phaseimg(1:end-1,2:end)-phaseimg(2:end,1:end-1);

if (nargin < 5)
    is8connect = 0;
end

nConnect = 4+is8connect*4;

if (nargin < 7)
    probXs = 1;
end
    probXd  = 1-probXs;

diagFactor = 1;
df = 1/diagFactor;


% Using smooth surface assumption
if (nargin < 7)

% Various arguments
if (nargin < 2)
    sigma = (var([vertshift(:)' horshift(:)']))^.5;
    sigmaDiag = sigma;%(var([diagshiftUpLeft(:)' diagshiftDownLeft(:)']))^.5;
else
    sigmaDiag = sigma*diagFactor;
end

if (nargin < 3)
    measure = @L2;
end

if (nargin < 4)
    truncate = 0;
end

if (nargin < 2)
    sigma = (var([vertshift(:)' horshift(:)']))^.5;
    sigmaDiag = (var([diagshiftUpLeft(:)' diagshiftDownLeft(:)']))^.5;
else
    sigmaDiag = sigma*1;
end

if is8connect
    vE = feval(measure,vertshift,sigma,truncate);
    hE = feval(measure,horshift,sigma,truncate);
    vdE = feval(measure,diagshiftUpLeft,sigma,truncate);
    hdE = feval(measure,diagshiftDownLeft,sigma,truncate);
    energy = sum(vE(:))+sum(hE(:))+sum(vdE(:))+sum(hdE(:));
    emap = vE; emap(1:end+1,1:end-1,2) = hE;
    emap(1:end-1,1:end-1,3) = vdE;emap(1:end-1,1:end-1,4) = hdE;
else
    vE = feval(measure,vertshift,sigma,truncate);
    hE = feval(measure,horshift,sigma,truncate);
    energy = sum(vE(:))+sum(hE(:));
    emap = vE; emap(1:end+1,1:end-1,2) = hE;
end

else
    
img = phaseimg;
phaseimg = mod(phaseimg,1);

% First find the joint probabilities

if is8connect
shifts = zeros([nConnect Ny Nx]);
shifts(2,2:Ny,:) = phaseimg(1:end-1,:)-phaseimg(2:end,:);
shifts(8,:,2:Nx) = phaseimg(:,1:end-1)-phaseimg(:,2:end);
shifts(6,1:Ny-1,:) = phaseimg(2:end,:)-phaseimg(1:end-1,:);
shifts(4,:,1:Nx-1) = phaseimg(:,2:end)-phaseimg(:,1:end-1);
shifts(1,2:Ny,2:Nx) = (phaseimg(1:end-1,1:end-1)-phaseimg(2:end,2:end))*df;
shifts(7,1:Ny-1,2:Nx) = (phaseimg(2:end,1:end-1)-phaseimg(1:end-1,2:end))*df;
shifts(5,1:Ny-1,1:Nx-1) = (phaseimg(2:end,2:end)-phaseimg(1:end-1,1:end-1))*df;
shifts(3,2:Ny,1:Nx-1) = (phaseimg(1:end-1,2:end)-phaseimg(2:end,1:end-1))*df;

shiftmap = zeros([nConnect Ny Nx]);
shiftmap(2,2:Ny,:) = 1;
shiftmap(8,:,2:Nx) = 1;
shiftmap(6,1:Ny-1,:) = 1;
shiftmap(4,:,1:Nx-1) = 1;
shiftmap(1,2:Ny,2:Nx) = 1;
shiftmap(7,1:Ny-1,2:Nx) = 1;
shiftmap(5,1:Ny-1,1:Nx-1) = 1;
shiftmap(3,2:Ny,1:Nx-1) = 1;

else % Use only 4-connect adjacencies 
    
shifts = zeros([nConnect Ny Nx]);
shifts(1,2:Ny,:) = phaseimg(1:end-1,:)-phaseimg(2:end,:);
shifts(4,:,2:Nx) = phaseimg(:,1:end-1)-phaseimg(:,2:end);
shifts(3,1:Ny-1,:) = phaseimg(2:end,:)-phaseimg(1:end-1,:);
shifts(2,:,1:Nx-1) = phaseimg(:,2:end)-phaseimg(:,1:end-1);

shiftmap = zeros([nConnect Ny Nx]);
shiftmap(1,2:Ny,:) = 1;
shiftmap(4,:,2:Nx) = 1;
shiftmap(3,1:Ny-1,:) = 1;
shiftmap(2,:,1:Nx-1) = 1;
end

jointprob = zeros([size(shiftmap) nStates*nStates]);

for p = 0:nStates
for q = 0:nStates
    jointprob(:,:,:,q*nStates+p+1) = feval(measure,(shifts + q - p),sigma,truncate).*shiftmap;
end
end

jointprob = sum(exp(-jointprob),4)*probXs/(sqrt(2*pi)*sigma);

% Now evaluate the probability of the proposed solution

phaseimg = img;
phaseimg = mod(phaseimg,1);

if is8connect
shifts(2,2:Ny,:) = phaseimg(1:end-1,:)-phaseimg(2:end,:);
shifts(8,:,2:Nx) = phaseimg(:,1:end-1)-phaseimg(:,2:end);
shifts(6,1:Ny-1,:) = phaseimg(2:end,:)-phaseimg(1:end-1,:);
shifts(4,:,1:Nx-1) = phaseimg(:,2:end)-phaseimg(:,1:end-1);
shifts(1,2:Ny,2:Nx) = (phaseimg(1:end-1,1:end-1)-phaseimg(2:end,2:end))*df;
shifts(7,1:Ny-1,2:Nx) = (phaseimg(2:end,1:end-1)-phaseimg(1:end-1,2:end))*df;
shifts(5,1:Ny-1,1:Nx-1) = (phaseimg(2:end,2:end)-phaseimg(1:end-1,1:end-1))*df;
shifts(3,2:Ny,1:Nx-1) = (phaseimg(1:end-1,2:end)-phaseimg(2:end,1:end-1))*df;

else % Use only 4-connect adjacencies 
    
shifts(1,2:Ny,:) = phaseimg(1:end-1,:)-phaseimg(2:end,:);
shifts(4,:,2:Nx) = phaseimg(:,1:end-1)-phaseimg(:,2:end);
shifts(3,1:Ny-1,:) = phaseimg(2:end,:)-phaseimg(1:end-1,:);
shifts(2,:,1:Nx-1) = phaseimg(:,2:end)-phaseimg(:,1:end-1);
end

margprob  = exp(-feval(measure,(shifts),sigma,truncate)).*shiftmap*probXs/(sqrt(2*pi)*sigma);

emap = ((margprob + probXd/(nStates+1))./(jointprob+probXd)).*shiftmap;

energy = sum(emap(:))/sum(shiftmap(:));
    
end