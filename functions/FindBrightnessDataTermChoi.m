function [dataterm] = FindBrightnessDataTermChoi(phaseimg,ab,badpixelthresh,maxwrap,phase2dist,tol,maxiter,useGrabCut)
% Data term computation by Choi et al's method
% Phase should be normalized to range from 0 to 1 instead of -pi to pi or 0 to 2pi

if nargin < 6
    tol = 1e-7;
end
if nargin < 7
    maxiter = 500;
end
if nargin < 8
    useGrabCut = 0;
end

numStates = maxwrap+1;
[M N] = size(phaseimg);

r = phaseimg;
Rmod = phase2dist;


% Whole image 
 % in case phase hasn't been wrapped and modded
 if (max(r(:)) > 1)
     RWrap = mod(r/Rmod,1);
 else
     RWrap = r;
 end
abc = ab.*RWrap.^2;


% Masked pixels
if length(badpixelthresh) == 1
    abmask = ab>badpixelthresh;
elseif length(badpixelthresh) == 2
    abmask = ab>badpixelthresh(1) & ab<badpixelthresh(2);
else
    abmask = ones(size(r));
end

 % in case phase hasn't been modded and wrapped
X = r(abmask);
if (max(r(:) > 100))
    wD = mod(X/Rmod,1);
else
    wD = mod(X,1);
end
    
x = ab(abmask).*wD.^2;
maxab = max(ab(:));
p = 2;

%[p,m,sigma,pkn,niter]=em(x,p,m,sigma,tol,maxiter)
[p,m,sigma,pkn,niter]=em(x',p,[],[],tol,maxiter);

if isnan(p)
    dataterm = zeros([size(ab) maxwrap+1]);
    return
end

if m(1) > m(2)
    h = 1;
    l = 2;
else
    h = 2;
    l = 1;
end

H = x(pkn(h,:)>=pkn(l,:));
L = x(pkn(h,:)<pkn(l,:));
Hthresh = (min(H)+max(L))/2;
% For the record, this is how probabilities are determined:
% p1 = p(1)*exp(-(x(1)-m(1))^2/(2*sigma(1)^2)) / (sigma(1)*sqrt(2*pi));
% p2 = p(2)*exp(-(x(1)-m(2))^2/(2*sigma(2)^2)) / (sigma(2)*sqrt(2*pi));

% GrabCut to refine threshold results
if (useGrabCut)
    diffThreshold = 1e-3;
    Beta = .1;
    maxIterations = 10;
    G = 50;
    Kclusters = 5;

    finalLabel = GCAlgo( abc, abc>Hthresh,  Kclusters, G, maxIterations, Beta, diffThreshold, [] );


    % Build dataterm
    % label the image U, H, L
    Umask = ~abmask;  % to be expanded by adding borders and other regions
    Lmask = ~finalLabel & ~Umask;
    Hmask = ~(Umask | Lmask);
 else

    % Build dataterm
    % label the image U, H, L
    Umask = ~abmask;  % to be expanded by adding borders and other regions
    Lmask = abc < Hthresh & ~Umask;
    Hmask = ~(Umask | Lmask);
end

% Find P(H|I_i) and P(L|I_i) for all pixels
dataterm = zeros([size(abc) maxwrap+1]);
pH = p(h)*exp(-(abc-m(h)).^2/(2*sigma(h)^2)) / (sigma(h)*sqrt(2*pi));
pL = p(l)*exp(-(abc-m(l)).^2/(2*sigma(l)^2)) / (sigma(l)*sqrt(2*pi));
pI = pH + pL;
pH = pH./pI; 
pL = pL./pI;

% Create the data terms based on H, L and U, then combine them using masks
% For H, dt = 1-pH if n=0, 1 otherwise
dtH = 1 - pH;
dtH(:,:,2:maxwrap+1) = 1;
% For L, dt = 1 if n=0, 1-pL/maxwrap otherwise
dtL = ones(size(abc));
dtL(:,:,2:maxwrap+1) = repmat(1-pL/maxwrap,[1 1 maxwrap]);
% For U, dt = 1
dtU = ones([size(abc) maxwrap+1]);

dataterm = repmat(Umask,[1 1 maxwrap+1]).*dtU + ...
           repmat(Hmask,[1 1 maxwrap+1]).*dtH + ...
           repmat(Lmask,[1 1 maxwrap+1]).*dtL;
       
%dataterm = dataterm./repmat(sum(dataterm,3),[1 1 maxwrap+1]);