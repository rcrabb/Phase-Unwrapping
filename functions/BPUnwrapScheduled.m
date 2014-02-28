function [WrapState Energy Steps Shifts nRes ] = BPUnwrapScheduled(phaseimg,sigma,nStates,measure,truncate,nConnect)


MAX_ITER = 400;
MAX_STATES = 2;
diagFactor = 1;
df = 1/diagFactor;
endEarly = 0;
[Ny, Nx] = size(phaseimg);
isScheduled = 1;

% Various arguments
if (nargin < 2)
    sigma = (var([vertshift(:)' horshift(:)']))^.5;
    sigmaDiag = sigma;%(var([diagshiftUpLeft(:)' diagshiftDownLeft(:)']))^.5;
else
    sigmaDiag = sigma*diagFactor;
end

if (nargin < 3)
    K = MAX_STATES;
else
    K = nStates; 
end

if (nargin < 4)
    measure = @L2;
end

if (nargin < 5)
    truncate = 0;
end

if (nargin < 6)
    nConnect =  4;
end


if nConnect == 8

scheduleLabels = [{'down'},{'up'},{'downright'},{'upleft'},{'right'},{'left'},{'upright'},{'downleft'}];
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

scheduleLabels = [{'up'},{'right'},{'down'},{'left'}];
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

messageOut = zeros([nConnect Ny Nx K+1]);
messageIn0 = zeros([nConnect Ny Nx K+1 K+1]);
marginalized = zeros([Ny Nx K+1]);
unmarginalized = zeros([Ny Nx K+1]);

% 1 2 3
% 8 * 4
% 7 6 5
            
% initialize to uniform distribution

%            min_p ( Pr(p,q) * ?_(s in N(p)\q) m_sp(fs) )
% m_pq(fq) = -----------------------------------------------
%                     ?_fq p(fq)

for q = 1:K+1        
    for p = 1:K+1
        messageIn0(:,:,:,q,p) = feval(measure,(shifts + q - p),sigma,truncate).*shiftmap;
    end
end

%messageIn0 = min(messageIn0,.75 / (2*sigma));

for iter = 1:MAX_ITER
    
 if (isScheduled && iter > 1)
    for sdx = [3 1 2 4]
        schedule = scheduleLabels{sdx};
    switch(schedule)
        case 'down'
            lineIndex = 1:Ny;
        case 'up'
            lineIndex = Ny:-1:1;
        case 'right'
            lineIndex = 1:Nx;
        case 'left'
            lineIndex = Nx:-1:1;
            
    end
    
    for ldx = lineIndex
        
      switch(schedule)
          case 'down'
              if ldx == Ny
                  mAdd = unmarginalized(ldx,:,:);
              else
                  mAdd = squeeze(unmarginalized(ldx,:,:)) - squeeze(messageOut(1,ldx+1,:,:));
              end
              
              mIn = messageIn0(sdx,ldx,:,:,:) + repmat(reshape(mAdd,[1 1 Nx 1 K+1]),[1 1 1 K+1 1]);
              mIn = min(mIn,[],5);
              mIn = mIn + repmat( log(sum(exp(-mIn),4)), [1 1 1 K+1]);
              
              % Copy messages to outgoing
              if (ldx < Ny)
                  unmarginalized(ldx+1,:,:) = unmarginalized(ldx+1,:,:) - reshape(squeeze(messageOut(sdx,ldx,:,:)) - squeeze(mIn), [1 Nx K+1]);
              end
              messageOut(sdx,ldx,:,:) = mIn;
              
          case 'up'
              if ldx == 1
                  mAdd = unmarginalized(ldx,:,:);
              else
                  mAdd = squeeze(unmarginalized(ldx,:,:)) - squeeze(messageOut(3,ldx-1,:,:));
              end
              
              mIn = messageIn0(sdx,ldx,:,:,:) + repmat(reshape(mAdd,[1 1 Nx 1 K+1]),[1 1 1 K+1 1]);
              mIn = min(mIn,[],5);
              mIn = mIn + repmat( log(sum(exp(-mIn),4)), [1 1 1 K+1]);
              
              % Copy messages to outgoing
              if (ldx >1)
                  unmarginalized(ldx-1,:,:) = unmarginalized(ldx-1,:,:) - reshape(squeeze(messageOut(sdx,ldx,:,:)) - squeeze(mIn), [1 Nx K+1]);
              end
              messageOut(sdx,ldx,:,:) = mIn;
          case 'right'
              if ldx == Nx
                  mAdd = unmarginalized(:,ldx,:);
              else
                  mAdd = squeeze(unmarginalized(:,ldx,:)) - squeeze(messageOut(4,:,ldx+1,:));
              end
              
              mIn = messageIn0(sdx,:,ldx,:,:) + repmat(reshape(mAdd,[1 Ny 1 1 K+1]),[1 1 1 K+1 1]);
              mIn = min(mIn,[],5);
              mIn = mIn + repmat( log(sum(exp(-mIn),4)), [1 1 1 K+1]);
              
              % Copy messages to outgoing
              if (ldx < Nx)
                  unmarginalized(:,ldx+1,:) = unmarginalized(:,ldx+1,:) - reshape(squeeze(messageOut(sdx,:,ldx,:)) - squeeze(mIn),[Ny 1 K+1]);
              end
              messageOut(sdx,:,ldx,:) = mIn;
              
          case 'left'
              if ldx == 1
                  mAdd = unmarginalized(:,ldx,:);
              else
                  mAdd = squeeze(unmarginalized(:,ldx,:)) - squeeze(messageOut(2,:,ldx-1,:));
              end
              
              mIn = messageIn0(sdx,:,ldx,:,:) + repmat(reshape(mAdd,[1 Ny 1 1 K+1]),[1 1 1 K+1 1]);
              mIn = min(mIn,[],5);
              mIn = mIn + repmat( log(sum(exp(-mIn),4)), [1 1 1 K+1]);
              
              % Copy messages to outgoing
              if (ldx > 1)
                  unmarginalized(:,ldx-1,:) = unmarginalized(:,ldx-1,:) - reshape(squeeze(messageOut(sdx,:,ldx,:)) - squeeze(mIn),[Ny 1 K+1]);
              end
              messageOut(sdx,:,ldx,:) = mIn;
      end
    end
    end
    
    marginalized = unmarginalized + repmat( log(sum(exp(-unmarginalized),3)),[1 1 K+1]);   
    
 else
if nConnect == 8
     messageAdd = repmat(reshape(unmarginalized,[1 Ny Nx K+1]),[nConnect 1 1 1]) - [... 
                   reshape(padarray(squeeze(messageOut(5,1:end-1,1:end-1,:)) ,[1 1],0,'pre') ,[1 Ny Nx K+1])  ; ...
                   reshape(padarray(squeeze(messageOut(6,1:end-1,:,:)) ,[1 0],0,'pre') ,[1 Ny Nx K+1])  ; ...
          reshape(padarray(padarray(squeeze(messageOut(7,1:end-1,2:end,:)) ,[0 1],0,'post'), [1 0],'pre'),[1 Ny Nx K+1]) ;...
                   reshape(padarray(squeeze(messageOut(8,:,2:end,:)) ,[0 1],0,'post') ,[1 Ny Nx K+1]) ; ...
                   reshape(padarray(squeeze(messageOut(1,2:end,2:end,:)) ,[1 1],0,'post') ,[1 Ny Nx K+1])  ; ...
                   reshape(padarray(squeeze(messageOut(2,2:end,:,:)) ,[1 0],0,'post') ,[1 Ny Nx K+1])  ; ...
          reshape(padarray(padarray(squeeze(messageOut(3,2:end,1:end-1,:)) ,[0 1],0,'pre'), [1 0],'post') ,[1 Ny Nx K+1])  ; ...
                   reshape(padarray(squeeze(messageOut(4,:,1:end-1,:)) ,[0 1],0,'pre') ,[1 Ny Nx K+1])   ...
          ];
else
    messageAdd = repmat(reshape(unmarginalized,[1 Ny Nx K+1]),[nConnect 1 1 1]) - [... 
                   reshape(padarray(squeeze(messageOut(3,1:end-1,:,:)) ,[1 0],0,'pre') ,[1 Ny Nx K+1])  ; ...
                   reshape(padarray(squeeze(messageOut(4,:,2:end,:)) ,[0 1],0,'post') ,[1 Ny Nx K+1]) ; ...
                   reshape(padarray(squeeze(messageOut(1,2:end,:,:)) ,[1 0],0,'post') ,[1 Ny Nx K+1])  ; ...
                   reshape(padarray(squeeze(messageOut(2,:,1:end-1,:)) ,[0 1],0,'pre') ,[1 Ny Nx K+1])   ...
          ];
end
    
    
    messageIn = messageIn0 + repmat(reshape(messageAdd,[nConnect Ny Nx 1 K+1]),[1 1 1 K+1 1]);    
    messageIn = min(messageIn,[],5);
    
    % Normalize
    messageIn = messageIn + repmat( log(sum(exp(-messageIn),4)), [1 1 1 K+1]);
    messageOut = messageIn;

    % Compute marginal probabilities
if nConnect == 8
    unmarginalized = padarray(squeeze(messageOut(2,2:end,:,:)) ,[1 0],0,'post')   + ...
                   padarray(squeeze(messageOut(8,:,2:end,:)) ,[0 1],0,'post')   + ...
                   padarray(squeeze(messageOut(6,1:end-1,:,:)) ,[1 0],0,'pre')   + ...
                   padarray(squeeze(messageOut(4,:,1:end-1,:)) ,[0 1],0,'pre')   + ...
                   padarray(squeeze(messageOut(1,2:end,2:end,:)) ,[1 1],0,'post')   + ...
          padarray(padarray(squeeze(messageOut(3,2:end,1:end-1,:)) ,[0 1],0,'pre'), [1 0],'post')   + ...
                   padarray(squeeze(messageOut(5,1:end-1,1:end-1,:)) ,[1 1],0,'pre')   + ...
          padarray(padarray(squeeze(messageOut(7,1:end-1,2:end,:)) ,[0 1],0,'post'), [1 0],'pre');
    marginalized = unmarginalized + repmat( log(sum(exp(-unmarginalized),3)),[1 1 K+1]);
else
    unmarginalized = padarray(squeeze(messageOut(1,2:end,:,:)) ,[1 0],0,'post')   + ...
                   padarray(squeeze(messageOut(4,:,2:end,:)) ,[0 1],0,'post')   + ...
                   padarray(squeeze(messageOut(3,1:end-1,:,:)) ,[1 0],0,'pre')   + ...
                   padarray(squeeze(messageOut(2,:,1:end-1,:)) ,[0 1],0,'pre');
    marginalized = unmarginalized + repmat( log(sum(exp(-unmarginalized),3)),[1 1 K+1]);    
end
 end
    
    [ml, map] = min(marginalized, [], 3);
    WrapState(:,:,iter) = map;
    
    vshift(:,:,iter) = map(1:Ny-1,:)-map(2:Ny,:);
    hshift(:,:,iter) = map(:,1:Nx-1)-map(:,2:Nx);
    
    % Compute Energy
    [nrg emap] = BPEnergy(phaseimg+map,sigma,measure,truncate,nConnect==8);
    Energy(iter) = nrg;
    curlsum = vshift(:,1:end-1,iter)-vshift(:,2:end,iter)-hshift(1:end-1,:,iter)+hshift(2:end,:,iter);
    nRes(iter) = sum(abs(curlsum(:))>1e-15);

    % check for convergence
    if iter > 1
        %sum(vec(WrapState(:,:,iter)-WrapState(:,:,iter-1)))
        %sum(abs(marginalized(:)-oldmarg(:)))
        if ( (sum(abs(marginalized(:)-oldmarg(:))) < 1e-10) &&...
             (sum(vec(WrapState(:,:,iter)-WrapState(:,:,iter-1))) == 0) )
            disp( sprintf('Converged at %d iterations',iter));
            endEarly = 1;
            break;
        end
    end
    oldmarg = marginalized;
end

vertshift = phaseimg(2:end,:)-phaseimg(1:end-1,:);
horshift = phaseimg(:,2:end)-phaseimg(:,1:end-1);
if nConnect == 8
    diagshiftUpLeft = phaseimg(2:end,2:end)-phaseimg(1:end-1,1:end-1);
    diagshiftDownLeft = phaseimg(1:end-1,2:end)-phaseimg(2:end,1:end-1);
end
Shifts.v = repmat(vertshift,[1 1 size(vshift,3)])-vshift;
Shifts.h = repmat(horshift,[1 1 size(hshift,3)])-hshift;
Steps.v = vshift;
Steps.h = hshift;
if ~endEarly
    disp( sprintf('No convergence after %d iterations',iter));
end

end