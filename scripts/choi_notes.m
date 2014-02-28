%
% Load data
dataname = 'visionlab01';
[R AB] = pgm2mat(['C:\data\CVPRdata\' dataname]);
r = median(R,3);
ab = median(AB,3);
clear R AB

%% parameter settings
phase2dist = 2500;
abthresh = 50;
maxwrap = ceil(max(r(:))/phase2dist);%3;
tol = .0001;
maxiter = 500;


%% Choi notes (now contained in FindBrightnessDataTermChoi()
% Whole image 
phaseimg = mod(r/phase2dist,1);
WrapState = floor(r/phase2dist);
abc = ab.*phaseimg.^2;

% Masked pixels
abmask = ab>abthresh;
X = r(abmask);
wD = mod(X/phase2dist,1);
xwrapstate = floor(X/phase2dist);
x = ab(abmask).*wD.^2;
maxab = max(ab(:));
p = 2;
m = double(sum(xwrapstate>0))/length(x);
m(2) = 1-m(1);

%[p,m,sigma,pkn,niter]=em(x,p,m,sigma,tol,maxiter)
[p,m,sigma,pkn,niter]=em(x',p,[],[],tol,maxiter);
% or
%[label, model, llh] = emgm(x',2);
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

if (mean(L) > mean(H))
    disp('uhoh - mean Low is greater than High');
end

%% Use GrabCut to further segment
diffThreshold = 1e-3;
Beta = .1;
maxIterations = 10;
G = 50;
Kclusters = 5;

finalLabel = GCAlgo( abc, abc>Hthresh,  Kclusters, G, maxIterations, Beta, diffThreshold, [] );
%% Build dataterm
% label the image U, H, L
Umask = ~abmask;  % to be expanded by adding borders and other regions
Lmask = ~finalLabel & ~Umask;%abc < Hthresh & ~Umask;
Hmask = ~(Umask | Lmask);

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
       
dtnormalized = dataterm./repmat(sum(dataterm,3),[1 1 maxwrap+1]);

%% Display the H, L labeling
binsize = 20;
hbins = binsize/2:binsize:roundn(max(H)+binsize,2);
figure;hist(H,hbins);
 hold on;
 hist(L,hbins); 
hobj = findobj(gca,'Type','patch');
%display(h);
 
set(hobj(1),'FaceColor','b','EdgeColor','k');
set(hobj(2),'FaceColor','g','EdgeColor','k');
title('Histogram of corrected amplitude values');

% create color image showing low and high
hsvimg = ones([size(ab) 3]);
hsvimg(:,:,3) = ((ab+maxab/4).*abmask + ab.*~abmask)/(maxab*1.25);
hsvimg(:,:,2) = .75;
hsvimg(:,:,1) = 1*double(Hmask) + .3*double(Lmask);
labeledab = hsv2rgb(hsvimg);
figure;imshow(labeledab);
title({'High amplitude: red,  Low amplitude: green','Overlayed on intensity image'})
figure;seq(mod(r/phase2dist,1));
title('Wrapped phase image (input)');
figure;seq(ab);
title('Original intensity image (input)');
figure;seq(abc);
title('Corrected intensity image (input)');


hsvwrap(:,:,3) = mod(r/phase2dist,1)*3/4+.25;
hsvwrap(:,:,2) = .5;
hsvwrap(:,:,1) = 1*double(Hmask) + .3*double(Lmask);
labeledwrap = hsv2rgb(hsvwrap);
figure;imshow(labeledwrap);
title({'High amplitude: red,  Low amplitude: green','Overlayed on wrapped phase image'})

figure;seq(WrapState);
title('Ground truth number of wraps');
colorbar('YTick',0:max(WrapState(:)));
correctpixels = double((Lmask == WrapState>0).*abmask)*2+double(~abmask);
figure;seq(correctpixels);
title('Correct pixels: white, Incorrect: black, Excluded: gray');

figure;zseq(phaseimg);
title('Wrapped Phase Image');

%% Compare results of data term between mine and Choi's method

% decrease image size for faster testing
rfactor = .5;

load('BrightnessConst20130222.mat');

% Load data
[R AB] = pgm2mat(['C:\data\ChoiTest\' dataname]);
r = median(R,3);
ab = median(AB,3);

clear R AB

%%
% parameter settings
phase2dist = 2500;
abthresh = 50;
maxwrap = ceil(max(r(:))/phase2dist);%3;
tol = .0001;
maxiter = 500;
nConnect = 8;
measure = @L2;
sigma = 1;
truncate = .5;
gamma = 1;
tau = 2000;
beta = 1;
eta = 1e-6;
dtweight = 1;

phaseimg = mod(r/phase2dist,1);
WrapStateGT = floor(r/phase2dist);

dataterm = FindBrightnessDataTerm(phaseimg,ab,BrightnessConst,maxwrap,phase2dist);
dataterm = -log((dataterm+eta)./repmat(sum(dataterm+eta,3),[1 1 maxwrap+1]));
datatermChoi = FindBrightnessDataTermChoi(phaseimg,ab,abthresh,maxwrap,phase2dist,tol,maxiter);

% resize inputs
dataterm=imresize(dataterm,rfactor,'nearest');
datatermChoi=imresize(datatermChoi,rfactor,'nearest');
phaseimg =  imresize(phaseimg,rfactor,'nearest');
WrapStateGT = imresize(WrapStateGT,rfactor,'nearest');

%% Test the data and smoothness terms alone
[~, dtsol] = min(dataterm,[],3);
dtsol = dtsol-1;
pctcorrect = sum(vec(dtsol==WrapStateGT))/length(WrapStateGT(:));

[~, dtsolchoi] = min(datatermChoi,[],3);
dtsolchoi = dtsolchoi-1;
pctcorrectChoi = sum(vec(dtsolchoi==WrapStateGT))/length(WrapStateGT(:));

smoothnessterm = FindSmoothnessTerm(phaseimg,maxwrap,nConnect,measure,sigma,truncate);

%smoothnesstermChoi = FindSmoothnessTermChoi(phaseimg,maxwrap,phase2dist,nConnect,gamma,tau,beta);
%smoothnesstermChoi = FindSmoothnessTermChoi(phaseimg,maxwrap,phase2dist,nConnect);
smoothnesstermChoi = FindSmoothnessTermChoi(phaseimg,maxwrap,phase2dist,nConnect,ChoiParam.gamma,ChoiParam.tau,ChoiParam.beta);


%BPUnwrap(phaseimg,sigma,nStates,measure,truncate,nConnect,probXs,dataterm)
[WrapState Energy Steps Shifts nRes ] = BPUnwrap(phaseimg,ChoiParams,maxwrap,@FindSmoothnessTermChoi,truncate,nConnect,1,dtweight*datatermChoi);
pctcorrectChoi = sum(vec(WrapState(:,:,end)==WrapStateGT))/length(WrapStateGT(:));


%% parameter loop for Choi
dtterms = 1:2;
dtweights = 10.^(-2:2);
taus = [.001 1 5 15 50 500];%10.^(2:2:8);
gammas = 10.^(-4:2:2);
betas = 10.^(-3:2);
% %change this if you want to start over%
if 1
    iter = 0;
    clear choiendstate ChoiParams
else
    iter = length(ChoiParams);
end

inititer = iter;
totaltests = length(taus)*length(gammas)*length(betas)*length(dtweights)*length(dtterms);
h = waitbar(0,'Running loopy belief propagation for Choi');

for dtterm = dtterms
    if dtterm == 1
        datatermvar = dataterm;
    elseif dtterm == 2
        datatermvar = datatermChoi;
    end
for dtweight = dtweights
for tau = taus
for gamma =  gammas
for beta = betas
    iter = iter+1;
    ChoiParams(iter).tau = tau;
    ChoiParams(iter).gamma = gamma;
    ChoiParams(iter).beta = beta;
    ChoiParams(iter).phase2dist = phase2dist;
    ChoiParams(iter).maxwrap = maxwrap;
    ChoiParams(iter).dtterm = dtterm;
    command = '[WrapState,~,~,~,~] = BPUnwrap(phaseimg,ChoiParams(iter),maxwrap,@FindSmoothnessTermChoi,truncate,nConnect,1,dtweight*datatermvar);';
    T = evalc(command);
    choiendstate(:,:,iter) = WrapState(:,:,end);
    ChoiParams(iter).iterations = size(WrapState,3);
    ChoiParams(iter).States = uint8(WrapState);
    ChoiParams(iter).pctcorrect = sum(vec(WrapState(:,:,end)==WrapStateGT+1))/length(WrapStateGT(:));
    waitbar((iter-inititer)/totaltests,h);
end
end
save(['workspace_ChoiResults_' dataname num2str(phase2dist) '_' date '_mid.mat'],'ChoiParams','choiendstate','phaseimg','WrapStateGT','datatermChoi');
end
end
end
clear inititer totaltests T command
save(['workspace_ChoiResults_' dataname num2str(phase2dist) '_' date '.mat'],'ChoiParams','choiendstate','phaseimg','WrapStateGT','datatermChoi');
close(h);

%%
% parameter loop for Ours

sigmas = [ .1  1 ];
dtweights = [ .0001 .1 1 ];
measures = [{@L1} ];
m = 1;
trunc = [.02 .1 .5];%.5;%[ 0 .01 1 ];h

totaltests = length(measures)*length(trunc)*length(sigmas)*length(dtweights)/2;
h = waitbar(0,'Running loopy belief propagation');
% %change this if you want to start over%
if 0
    iter = 0;
    clear endstates s
else
    iter = length(s);
end

for dtweight = dtweights
for m = 1:1%length(measures)
for truncate = trunc
for sigma = sigmas
    iter = iter+1;   
    
    command = '[WrapState,~,~,~,~] = BPUnwrap(phaseimg,sigma,maxwrap,measures{m},truncate,nConnect,1,dataterm*dtweight)';
    T = evalc(command);
    
    s(iter).sigma = sigma;
    s(iter).maxwrap = maxwrap;
    s(iter).States = uint8(WrapState);
    s(iter).dtflag = dtweight>0;
    s(iter).m = measures{m};
    s(iter).truncate = truncate;
    s(iter).dtweight = dtweight;  
    endstates(:,:,iter) =  uint8(WrapState(:,:,end));
    s(iter).pctcorrect = sum(vec(WrapState(:,:,end)==WrapStateGT+1))/length(WrapStateGT(:));
    
    waitbar(iter/totaltests,h);
end
end
save(['workspace_ourresults_mid_' dataname num2str(phase2dist) date '.mat'],'s','endstates','phaseimg','WrapStateGT');

end
end
close(h);
save(['workspace_ourresults_' dataname num2str(phase2dist) date '.mat'],'s','endstates','phaseimg','WrapStateGT');

%% take a look at the highest performing results
iter = 0;
pctcorrectChoi = [ChoiParams(:).pctcorrect];
for i = 1:length(ChoiParams)
    if (  pctcorrectChoi(i) > .16)
        iter = iter+1;
        tmp.iterations = i;
        goodresults(iter) = ChoiParams(i);
    end
end

%% plot percentage correct vs beta
figure;plot([ChoiParams(:).pctcorrect], [ChoiParams(:).beta],'r.');
figure;plot([ChoiParams(:).pctcorrect], [ChoiParams(:).gamma],'r.');
figure;plot([ChoiParams(:).pctcorrect], [ChoiParams(:).tau],'r.');

%% take a look at the highest performing results
iter = 0;
for i = 1:length(s)
    if s(i).pctcorrect > .98
        iter = iter + 1;
        sgood(iter) = s(i);
        tind(iter) = i;
    end
end

%% Retest best cases
iter = 0;
for i = 1:length(sgood)
    iter = iter+1;
    sigma = sgood(i).sigma;
    dtweight = sgood(i).dtweight;
    meas = sgood(i).m;
    truncate = sgood(i).truncate;
    
    command = '[WrapState,~,~,~,~] = BPUnwrap(phaseimg,sigma,maxwrap,meas,truncate,nConnect,1,dataterm*dtweight)';
    T = evalc(command);
    s(iter) = sgood(i);
    s(iter).maxwrap = maxwrap;
    s(iter).States = uint8(WrapState);
    s(iter).dtflag = dtweight>0;
    endstates(:,:,iter) =  uint8(WrapState(:,:,end));
    s(iter).pctcorrect = sum(vec(WrapState(:,:,end)==WrapStateGT+1))/length(WrapStateGT(:));
end


%% Test many frequencies

% Load data
load('BrightnessConst20130222.mat');
[R AB] = pgm2mat('C:\data\ChoiTest\Dojo_debug1');
r = median(R,3);
ab = median(AB,3);

% parameter settings
rfactor = .25;
abthresh = 50;
tol = .0001;
nConnect = 8;
eta = 1e-6;

modRange = [100 500 1000 1500 2500];


for phase2dist = modRange

maxwrap = ceil(max(r(:))/phase2dist);

phaseimg = mod(r/phase2dist,1);
WrapStateGT = floor(r/phase2dist);

dataterm = FindBrightnessDataTerm(phaseimg,ab,BrightnessConst,maxwrap,phase2dist);
dataterm = -log((dataterm+eta)./repmat(sum(dataterm+eta,3),[1 1 maxwrap+1]));
datatermChoi = FindBrightnessDataTermChoi(phaseimg,ab,abthresh,maxwrap,phase2dist,tol,maxiter);

% resize inputs
dataterm=imresize(dataterm,rfactor,'nearest');
datatermChoi=imresize(datatermChoi,rfactor,'nearest');
phaseimg =  imresize(phaseimg,rfactor,'nearest');
WrapStateGT = imresize(WrapStateGT,rfactor,'nearest');

    % parameter loop for Choi
    dtweights = 1;
    taus = 10.^(0:2:6);
    gammas = 10.^(-4:2:0);
    betas = 10.^(-6:-4);
    % %change this if you want to start over%
    if 1
        iter = 0;
        clear choiendstate ChoiParams
    else
        iter = length(ChoiParams);
    end

    inititer = iter;
    totaltests = length(taus)*length(gammas)*length(betas)*length(dtweights);
    h = waitbar(0,'Running loopy belief propagation for Choi');

    for dtweight = dtweights
    for tau = taus
    for gamma =  gammas
    for beta = betas
        iter = iter+1;
        ChoiParams(iter).tau = tau;
        ChoiParams(iter).gamma = gamma;
        ChoiParams(iter).beta = beta;
        ChoiParams(iter).phase2dist = phase2dist;
        ChoiParams(iter).maxwrap = maxwrap;
        command = '[WrapState,~,~,~,~] = BPUnwrap(phaseimg,ChoiParams(iter),maxwrap,@FindSmoothnessTermChoi,truncate,nConnect,1,dtweight*datatermChoi);';
        T = evalc(command);
        choiendstate(:,:,iter) = WrapState(:,:,end);
        ChoiParams(iter).iterations = size(WrapState,3);
        ChoiParams(iter).States = uint8(WrapState);
        ChoiParams(iter).pctcorrect = sum(vec(WrapState(:,:,end)==WrapStateGT+1))/length(WrapStateGT(:));
        ChoiParams(iter).phaseimg = phaseimg;
        ChoiParams(iter).GT = WrapStateGT;
        waitbar((iter-inititer)/totaltests,h);
    end
    end
    save(['workspace_ChoiResults_' dataname num2str(phase2dist) '_' date '_mid.mat'],'ChoiParams','choiendstate','phaseimg','WrapStateGT');
    end
    end
    clear inititer totaltests T command
    save(['workspace_ChoiResults_' dataname num2str(phase2dist) '_' date '.mat'],'ChoiParams','choiendstate','phaseimg','WrapStateGT','datatermChoi');
    close(h);

    % parameter loop for Ours

    sigmas = [.001 .01 .1  1  10 100];
    dtweights = [0 .0001 .1 1 ];
    measures = [{@L1} ];
    m = 1;
    trunc = .5;%[ 0 .01 1 ];h

    totaltests = length(trunc)*length(sigmas)*length(dtweights);
    h = waitbar(0,'Running loopy belief propagation');
    % %change this if you want to start over%
    if 1
        iter = 0;
        clear endstates s
    else
        iter = length(s);
    end

    for dtweight = dtweights
    for m = 1:1%length(measures)
    for truncate = trunc
    for sigma = sigmas
        iter = iter+1;   

        command = '[WrapState,~,~,~,~] = BPUnwrap(phaseimg,sigma,maxwrap,measures{m},truncate,nConnect,1,dataterm*dtweight)';
        T = evalc(command);

        s(iter).sigma = sigma;
        s(iter).maxwrap = maxwrap;
        s(iter).States = uint8(WrapState);
        s(iter).dtflag = dtweight>0;
        s(iter).m = measures{m};
        s(iter).truncate = truncate;
        s(iter).dtweight = dtweight;
        s(iter).phaseimg = phaseimg;
        endstates(:,:,iter) =  uint8(WrapState(:,:,end));
        s(iter).pctcorrect = sum(vec(WrapState(:,:,end)==WrapStateGT+1))/length(WrapStateGT(:));

        waitbar(iter/totaltests,h);
    end
    end
    end
    end
    close(h);
    save(['workspace_ourresults_' dataname num2str(phase2dist) date '.mat'],'s','endstates','phaseimg','WrapStateGT','dataterm');

    
end