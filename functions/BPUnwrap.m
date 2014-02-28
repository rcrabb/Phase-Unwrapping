function [WrapState Energy Steps Shifts nRes ] = BPUnwrap(phaseimg,sigma,nStates,measure,truncate,nConnect,probXs,dataterm)


MAX_ITER = 1500;
MAX_STATES = 2;
diagFactor = 1;
df = 1/diagFactor;
endEarly = 0;
[Ny, Nx] = size(phaseimg);
vertshift = phaseimg(2:end,:)-phaseimg(1:end-1,:);
horshift = phaseimg(:,2:end)-phaseimg(:,1:end-1);

% Various arguments
if (nargin < 2)
    sigma = (var([vertshift(:)' horshift(:)']))^.5;
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

if (nargin < 7)
    probXs = 1;
    probXd = 0;
else
    if ~(probXs>0 && probXs < 1.0000001)
        return;
    end
    probXd  = 1-probXs;
end

if (nargin < 8)
    dataterm = zeros(nConnect,Ny,Nx,K+1);
else
    dataterm = repmat(reshape(dataterm,[1 Ny Nx K+1]),[nConnect 1 1 1]);
end

% 1 2 3
% 8 * 4
% 7 6 5

%Using proposed method or Choi
if ~isequal(measure,@FindSmoothnessTermChoi)

% initialize to uniform distribution

%            min_p ( Pr(p,q) * ?_(s in N(p)\q) m_sp(fs) )
% m_pq(fq) = -----------------------------------------------
%                     ?_fq p(fq)

    messageIn0 = FindSmoothnessTerm(phaseimg,K,nConnect,measure,sigma,truncate,probXs);

else
    % Use Choi's method
    try %to get parameters
        tau = sigma.tau;
        gamma = sigma.gamma;
        beta = sigma.beta;
        phase2dist = sigma.phase2dist;
    catch

        disp('unable to read parameters for Choi Smoothness function\n');
    end
    if (exist('tau') && exist('gamma') && exist('beta'))
        messageIn0 = FindSmoothnessTermChoi(phaseimg,K,phase2dist,nConnect,gamma,tau,beta);
    else
        disp('Make sure the paramters for tau, gamma and beta are stored in sigma variable\n');
        messageIn0 = FindSmoothnessTermChoi(phaseimg,K,phase2dist,nConnect);
    end
end
   
% NOW INCLUDED IN FindSmoothnessTerm() FUNCTION
% % if we are using formulation from WrapStateProb.pdf
% % -NOTE: phase is scaled to be 1 instead of 2pi, so does not appear in Xd
% if (probXs ~= 1)
%     messageIn0 = -log( ... 
%         (probXs/(sqrt(2*pi)*sigma)*exp(-messageIn0) + probXd/(K+1) ) ./...
%         (probXs/(sqrt(2*pi)*sigma)*repmat(sum(exp(-messageIn0),4),[1 1 1 K+1 1])+probXd)...
%         );
% end
    
messageOut = zeros([nConnect Ny Nx K+1]);
marginalized = zeros([Ny Nx K+1]);
unmarginalized = zeros([Ny Nx K+1]);

for iter = 1:MAX_ITER
    
    
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
    
    messageIn = messageIn0 + repmat(reshape(messageAdd+dataterm,[nConnect Ny Nx 1 K+1]),[1 1 1 K+1 1]) ;
    
    messageIn = min(messageIn,[],5);
    
    % Normalize
    messageIn = messageIn + repmat( log(sum(exp(-messageIn),4)), [1 1 1 K+1]);

    % Copy messages to outgoing
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
    
    [ml, map] = min(marginalized, [], 3);
    WrapState(:,:,iter) = map;
    
    vshift(:,:,iter) = map(1:Ny-1,:)-map(2:Ny,:);
    hshift(:,:,iter) = map(:,1:Nx-1)-map(:,2:Nx);
    
    % Compute Energy
    if ~isequal(measure,@FindSmoothnessTermChoi)
        [nrg emap] = BPEnergy(phaseimg+map,sigma,measure,truncate,nConnect==8,probXs,K);
    else
        nrg = 0; emap = 0;
    end
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
    if iter > 3
          if ( (sum(vec(WrapState(:,:,iter)~=WrapState(:,:,iter-3))) == 0) ...
            && (sum(vec(WrapState(:,:,iter)~=WrapState(:,:,iter-2))) == 0) ...
            && (sum(vec(WrapState(:,:,iter)~=WrapState(:,:,iter-1))) == 0) )
            disp( sprintf('Stable at %d iterations',iter));
            endEarly = 1;
            break;
          end
    end
    oldmarg = marginalized;
end

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
