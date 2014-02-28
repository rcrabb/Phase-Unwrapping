function [smoothnessterm] = FindSmoothnessTermChoi(phaseimg,maxwrap,phase2dist,nConnect,gamma,tau,beta)
% tau should be in same units as phase2dist
MAX_STATES = 2;
diagFactor = 1;
df = 1/diagFactor;
[Ny, Nx] = size(phaseimg);

% Various arguments

if (nargin < 2)
    K = MAX_STATES;
else
    K = maxwrap; 
end

if (nargin < 3)
    phase2dist = 1;
end

if (nargin < 4)
    nConnect =  4;
end

if (nargin < 5)
    gamma = 1;
end

if (nargin < 6)
    tau = 200;
end

if (nargin < 7)
    beta = 1;
end

% find the penalty

if nConnect == 8
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

Cs = cat(4,shifts-1, shifts, shifts+1);
[~, cstar] = min(abs(Cs),[],4);
[I1 I2 I3]=ndgrid(1:nConnect, 1:Ny,1:Nx);
minshift = Cs(sub2ind(size(Cs),I1,I2,I3,cstar));

penalty = gamma*exp(-beta*(phase2dist*minshift).^2).*double(phase2dist*minshift<tau).*shiftmap;

smoothnessterm = zeros([nConnect Ny Nx K+1 K+1]);

for q = 1:K+1        
    for p = 1:K+1
        %messageIn0(:,:,:,q,p) = feval(measure,(shifts + q - p),sigma,truncate).*shiftmap;
        smoothnessterm(:,:,:,q,p) = double(cstar-2~=q-p).*penalty;
    end
end