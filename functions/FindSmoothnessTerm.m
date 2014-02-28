function [smoothnessterm] = FindSmoothnessTerm(phaseimg,maxwrap,nConnect,measure,sigma,truncate,probXs)

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
    nConnect =  4;
end

if (nargin < 7)
    probXs = 1;
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

smoothnessterm = zeros([nConnect Ny Nx K+1 K+1]);


% 1 2 3
% 8 * 4
% 7 6 5
            
% initialize to uniform distribution

%            min_p ( Pr(p,q) * ?_(s in N(p)\q) m_sp(fs) )
% m_pq(fq) = -----------------------------------------------
%                     ?_fq p(fq)


for q = 1:K+1        
    for p = 1:K+1
        smoothnessterm(:,:,:,q,p) = feval(measure,(shifts + q - p),sigma,truncate).*shiftmap;
        
    end
end

if (probXs ~= 1)
    if (length(probXs) == 1)
        messageIn0 = -log( ... 
            (probXs/(sqrt(2*pi)*sigma)*exp(-messageIn0) + probXd/(K+1) ) ./...
            (probXs/(sqrt(2*pi)*sigma)*repmat(sum(exp(-messageIn0),4),[1 1 1 K+1 1])+probXd)...
            );
    else
        messageIn0 = -log( ...
            (probXs./(sqrt(2*pi)*sigma)*exp(-messageIn0) + probXd./(K+1) ) ./...
            (probXs./(sqrt(2*pi)*sigma)*repmat(sum(exp(-messageIn0),4),[1 1 1 K+1 1])+probXd)...
            );
    end
end