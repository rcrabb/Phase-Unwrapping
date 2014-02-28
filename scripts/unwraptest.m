% Generate test depth images
%% Test image with brightness term
tic
% mm per 1 phase
phase2dist = 2;
dist2phase = 1/phase2dist;
sigmas = [ .1 .5 1   5 100];
measures = [{@L1} {@L2}];
trunc = [ 0 .01 .5 1 ];
nConnect = 8;
probXs = 1;

nStates = 5;

A = gethlut();
[M N ~] = size(A);
% BrightnessConst = load(fullfile('C:\data\BrightnessCalSVC','collected_stats.mat'),'Clse');
% BrightnessConst = BrightnessConst.Clse;
load('BrightnessConst20130222.mat');

% testdir = 'C:\data\BrightnessPriorData\home2';
% Brightness = loadpgm(fullfile(testdir,'ConfidenceImage.pgm'));
% Z = loadpgm(fullfile(testdir,'ZImg.pgm'));
% AB = Brightness(:,:,1);
% Distance = Z(:,:,1)./A(:,:,3);
% phaseimg = Distance*dist2phase;

% Create DUCK test image
 % SEE CODE AT AROUND LINE 470
 load('ducksynth.mat');
% 
% ds =  load('ducksynth.mat');
% sf = 10;
% ds.depth = ds.phaseimg(1:sf:end,1:sf:end)*phase2dist;
% ds.XYZ = repmat(ds.depth,[1 1 3]).*A(1:sf:end,1:sf:end,:);
%   figure;mesh(ds.XYZ(:,:,1),ds.XYZ(:,:,2),10-ds.XYZ(:,:,3));

[dataterm] = FindBrightnessDataTerm(phaseimgwrap,AB,BrightnessConst,nStates,phase2dist);

% % %
% figure;
% for i = 1:6
%     subplot(6,2,i*2-1)
%     imshow(finstateAB(:,:,i));
%     set(gcf, 'Position', get(0,'Screensize'));
%    title(['Probability Wrap State is ' num2str(i-1)]);
% end

eta = .001;
dataterm = -log((dataterm+eta)./repmat(sum(dataterm+eta,3),[1 1 nStates+1]));
% dataterm = zeros(size(dataterm));
% to save time let's reduce the size
piw = imresize(phaseimgwrap,.5,'nearest');dt=imresize(dataterm,.5,'nearest');
pi =  imresize(phaseimg,.5,'nearest');gt = pi-piw;
[tmp dtgt] = min(dt,[],3);
tic

for i = 1:length(sigmas)
   
    [WrapState{i} Energy Steps Shifts nRes ] = BPUnwrap(piw,sigmas(i),nStates,measures{2},trunc(1),nConnect,probXs,zeros(size(dt)));

fnm = ['ws' num2str(i) '.png'];
ws = WrapState{i};
sum(vec(ws(:,:,end)==gt))/length(vec(ws(:,:,end)))
%imwrite(uint8(ws(:,:,end)*20),fnm);
Finstate(:,:,i) = ws(:,:,end);

    [WrapStateAB{i} Energy Steps Shifts nRes ] = BPUnwrap(piw,sigmas(i),nStates,measures{2},trunc(1),nConnect,probXs,dt);
toc
fnm = ['wsAB' num2str(i) '.png'];
ws = WrapStateAB{i};
sum(vec(ws(:,:,end)==gt))/length(vec(ws(:,:,end)))
% imwrite(uint8(ws(:,:,end)*20),fnm);
FinstateAB(:,:,i) = ws(:,:,end);


end
save workspace20130226


figure;
for i = 1:length(WrapState)
        subplot(length(WrapStateAB),2,i*2-1)
    imshow(finstateAB(:,:,i)/7);
    set(gcf, 'Position', get(0,'Screensize'));
   title(['Sigma: ' num2str(sigmas(i)) ' ; with convergence after' num2str(size(WrapStateAB{i},3)) ' frames']);

    subplot(length(WrapState),2,i*2)
    imshow(finstate(:,:,i)/7);
    set(gcf, 'Position', get(0,'Screensize'));
   title(['Sigma: ' num2str(sigmas(i)) ' ; with convergence after' num2str(size(WrapState{i},3)) ' frames']);
end

%% Load series of depth images
dimgs = dir('C:\Users\rcrabb\Documents\MATLAB\depthimages\depth*');
for i = 1:length(dimgs)
    img = imread(['C:\Users\rcrabb\Documents\MATLAB\depthimages\' dimgs(i).name]);
    if (size(img,3) == 3)
        img = double(rgb2gray(img));
    else
        img = double(img);
    end
    img = img./max(img(:));
    dimg(i) = {img};
end
save('C:\Users\rcrabb\Documents\MATLAB\depthimages\dimgs.mat','dimg');
%% Create data site
load('C:\Users\rcrabb\Documents\MATLAB\depthimages\dimgs.mat');
MaxPixelSize = 100;
phaserange = [ 1.5 3];
for i = 1:length(dimg)
    img = dimg{i};
    [m n] = size(img);
    factor = ceil(max(m/100,n/100));
    pi.scale = rand()*(phaserange(2)-phaserange(1))+phaserange(1);
    pi.gt = pi.scale*img(1:factor:m,1:factor:n);
    pi.phaseimg = mod(pi.gt,1);
    pi.gtstate = pi.gt-pi.phaseimg;
    %[Shifts Energy nRes Steps WrapState] = BPStateUnwrapQMinLog8connectFast(pi.phaseimg);
    %pi.States = WrapState;
    %pi.energy = BPEnergy(WrapState(:,:,end)+pi.phaseimg);
    %result = pi.gt;
    %[m n ] = size(pi.gt);
    %result(:,n+1:2*n) = pi.gtstate;
    %result(m+1:2*m,1:n) = pi.phaseimg + WrapState(:,:,end);
    %result(m+1:2*m,n+1:2*n) = WrapState(:,:,end);
    %pi.result = result;
    deptharray(i) = {pi};
end
clear MaxPixelSize phaserange i img factor pi m n
%%
for l = 1:length(d)
    pi = d{l};
    clear results
    for i = 1:length(pi.s)
        %eRatio(i,l) = pi.s(i).energy/pi.s(i).gtEnergy;
        result = pi.gt;
        [m n ] = size(pi.gt);
        result(:,n+1:2*n) = pi.gtstate;
        ws = double( pi.s(i).States(:,:,end));
        result(m+1:2*m,1:n) = pi.phaseimg + ws;
        result(m+1:2*m,n+1:2*n) = ws;
        results(:,:,i) = result;
    end
    figure;seq(results);
    pause;
end

%% Test a number of settings for multiple test cases
clear s
L = length(deptharray);
N = length(nState);
%sigmas = [.001 .01 .05 .1 .4 1 3];
sigmas = [.001 .01 .1 1 ];
S = length(sigmas);
measures = [{@L1} {@L2}];
M = length(measures);
trunc = [ 0 .01 .5 1 ];
T = length(trunc);
nState = [2 3 ];
nIterMax = N*L*S;
h = waitbar(0,'Running loopy belief propagation');
for l = 1:L
    iter = 0;
    clear s endstates results
    pi = deptharray{l};
    phaseimg = pi.phaseimg;
for m = 1:length(measures)
for truncate = trunc
for nStates = nState
for i = 1:length(sigmas)
    iter = iter+1;    
    [Shifts Energy nRes Steps States] = BPUnwrap(phaseimg,sigmas(i),nStates,measures{m},truncate);
    
    bpe = BPEnergy(pi.gt,sigmas(i),measures{m},truncate);

    s(iter).sigma = sigmas(i);
    s(iter).nStates = nStates;
%     s(iter).Shifts.h = Shifts.h;
%     s(iter).Shifts.v = Shifts.v;
    s(iter).Energy = Energy;
    s(iter).energy = BPEnergy(phaseimg+States(:,:,end),1,1);
    s(iter).gtEnergy = bpe;
    s(iter).nRes = nRes;
%     s(i).Steps = Steps;
    s(iter).States = uint8(States);
    %figure;seq(s(iter).States(:,:,end));
    
    % Show results image
    result = pi.gt;
    [m n ] = size(pi.gt);
    result(:,n+1:2*n) = pi.gtstate;
    ws = double( States(:,:,end));
    result(m+1:2*m,1:n) = pi.phaseimg + ws;
    result(m+1:2*m,n+1:2*n) = ws;
    results(:,:,iter) = result;
    pi.results = results;
        
    endstates(:,:,iter) = s(iter).States(:,:,end);
    pi.endstates = endstates;
    pi.s = s;
    deptharray(l) = {pi};
    waitbar(((l-1)*N*S+iter)/nIterMax,h);
    %title(sprintf('sigma: %f  - nStaes: %d',sigmas(i),nStates));
    clear Shifts Energy nRes Steps States
end
end 
end
end
end
close(h);
%save multitest.mat deptharray
% save in parts
beep;
%%
for j = 1:length(deptharray)
    d= deptharray{j};
    savename = ['depthimages\testresults' num2str(j,'%02d') '.mat'];
    save(savename,'d');
end
clear deptharray

%%
for j = 1:12
    savename = ['depthimages\testresults' num2str(j,'%02d') '.mat'];
    load(savename);
    for i = 1:length(d.s)
        eRatio(i,j) = d.s(i).energy/d.s(i).gtEnergy;
    end
    results{j} = d.result;
end
figure;plot(eRatio);
%% Test a number of settings 
clear s params prcCorrect sol
%sigmas = [.3 1 ];
sigmas = [.001 .01 .05 .1 .3 .7 1 2 3];
pXs = [.5 .9 .99 .9999 1];
%sigmas = [.15 .25 .35 .45 .55 .65 .75 .9];
%sigmas = [[.25 .35 .45 .55 .65 .75 .9],1+ [.25 .35 .45 .55 .65 .75 .9]];
nState = [1 2 ];
iter = 0;
for nStates = nState
for i = 1:length(sigmas)
for probXs = pXs    
    iter = iter+1;
    [WrapState Energy Steps Shifts nRes ] = BPUnwrap(phaseimg,sigmas(i),nStates,@L2,0,8,probXs);

    sol = WrapState(:,:,end);
    prcCorrect(iter) = max(max(max(sum(vec(sol==sol0)),sum(vec(sol==sol1))),...
        max(sum(vec(sol==sol3)),sum(vec(sol==sol2)))),...
        max(sum(vec(sol==sol4)),sum(vec(sol==sol5)))) / length(sol(:));
    params(iter,:) = [sigmas(i) nStates probXs prcCorrect(iter) Energy(end)];
%     s(iter).Shifts.h = Shifts.h;
%     s(iter).Shifts.v = Shifts.v;
%     s(iter).Energy = Energy;
%     s(iter).nRes = nRes;
%     %s(i).Steps = Steps;
%     s(iter).States = uint8(States);
%     clear Shifts Energy nRes Steps States
%     %figure;seq(s(iter).States(:,:,end));
%     endstates(:,:,iter) = s(iter).States(:,:,end);
%     %title(sprintf('sigma: %f  - nStaes: %d',sigmas(i),nStates));
end
end 
end
%save con8.mat s sigmas nState
beep;
%%
prms = params;
pc= prcCorrect;
prms(:,4) = pc'
prms(find(prcCorrect==max(prcCorrect)),:)
%%
for i = 1:iter 
    %unwrap = PU2D(double(s(i).Shifts.v(:,:,end)),double(s(i).Shifts.h(:,:,end)));
    %figure;seq(unwrap);
    figure;seq(s(i).States(:,:,end));
    title(sprintf('sigma: %f  - nStaes: %d',sigmas(mod(i-1,length(sigmas))+1),ceil((i-.1)/length(sigmas))));
    disp(sprintf('sigma: %f  - nStaes: %d',sigmas(mod(i-1,length(sigmas))+1),ceil((i-.1)/length(sigmas))));
end
beep;
%figure;plot(eRatio);
%%
for i = 1:length(deptharray)
    pi = deptharray{i};
    figure;seq(pi.result)
end

%% Gaussian bump
[X,Y] = meshgrid(-2:.2:2, -2:.2:2);
Z = X .* exp(-X.^2 - Y.^2);
Z = 3*(Z+.5);
surf(X,Y,Z);
gbumpsGT = Z;
gbumpsshiftGT = Phase2PhaseShifts(gbumpsGT);
gbumpswrap = mod(gbumpsGT,1);
tmp = Phase2PhaseShifts(gbumpswrap);
gbumpssol.h = gbumpsshiftGT.h - tmp.h;
gbumpssol.v = gbumpsshiftGT.v - tmp.v;
bumpsol = gbumpsGT-gbumpswrap;
%figure;seq(gbumpswrap);
figure;surf(X,Y,gbumpswrap);
[shiftgbumps gbE gbRes stepsgbumps] = loopyunwrapping(gbumpswrap);
[msShifts msEnergy msnRes msSteps Mc2s Ms2c] = BPUnwrap(gbumpswrap);
% simga .05-.3 seem to work
[usShifts usEnergy usnRes usSteps uStates] = BPStateUnwrap(gbumpswrap);
% simga .05-10 seem to work
[ppShifts ppEnergy ppnRes ppSteps ppStates] = BPStateUnwrapPMin(gbumpswrap, .1);
% simga .05-.3 seem to work
[pqShifts pqEnergy pqnRes pqSteps pqStates] = BPStateUnwrapQMin(gbumpswrap,.05);
[pqShifts pqEnergy pqnRes pqSteps pqStates] = BPStateUnwrapQMinLog8connect(gbumpswrap);
% simga .05-.3 seem to work
[pqsShifts pqsEnergy pqsnRes pqsSteps pqsStates] = BPStateUnwrapQSum(gbumpswrap,.05);

if (gbRes(end) == 0)
    gbumpunwrap = PU2D(shiftgbumps.v(:,:,end), shiftgbumps.h(:,:,end));
    figure;surf(X,Y,gbumpunwrap);
end

if (msnRes(end) == 0)
    gbumpunwrap = PU2D(msShifts.v(:,:,end), msShifts.h(:,:,end));
    figure;surf(X,Y,gbumpunwrap);
end


if (msnRes(end) == 0)
    munwrap = PU2D(msShifts.v(:,:,end), msShifts.h(:,:,end));
    figure;seq(munwrap);
end

if (usnRes(end) == 0)
    uunwrap = PU2D(usShifts.v(:,:,end), usShifts.h(:,:,end));
    figure;seq(uunwrap);
end

if (ppnRes(end) == 0)
    ppunwrap = PU2D(ppShifts.v(:,:,end), ppShifts.h(:,:,end));
    figure;seq(ppunwrap);
end

if (pqnRes(end) == 0)
    pqunwrap = PU2D(pqShifts.v(:,:,end), pqShifts.h(:,:,end));
    figure;seq(pqunwrap);
end

if (pqsnRes(end) == 0)
    pqunwrap = PU2D(pqsShifts.v(:,:,end), pqsShifts.h(:,:,end));
    figure;seq(pqunwrap);
end

%save states in a single images
bumpstates = uStates(:,:,1:35);
bumpstates(:,22:42,1:21) = ppStates;
bumpstates(:,22:42,22:35) = repmat(ppStates(:,:,end),[1 1 14]);
bumpstates(22:42,1:21,1:35) = pqsStates(:,:,1:35);
bumpstates(22:42,22:42,1:22) = pqStates;
bumpstates(22:42,22:42,23:35) = repmat(pqStates(:,:,end),[1 1 13]);

save bumpoutput.mat usShifts usEnergy usnRes usSteps uStates ppShifts ppEnergy ppnRes ppSteps ppStates pqShifts pqEnergy pqnRes pqSteps pqStates pqsShifts pqsEnergy pqsnRes pqsSteps pqsStates bumpstates bumpsol
clear tmp X Y Z

%% Simple Slope
slopeGT = repmat([0:.2:1.8]',[1 10]);
slopeGT = (slopeGT'+slopeGT)/2;
% %reverse!
% slopeGT = max(slopeGT(:))-slopeGT;
sshiftGT = Phase2PhaseShifts(slopeGT);
swrap = mod(slopeGT,1);
phaseimg = swrap;
tmp = Phase2PhaseShifts(swrap);
ssol.h = -(sshiftGT.h - tmp.h);
ssol.v = -(sshiftGT.v - tmp.v);
figure;seq(slopeGT);
[sShifts sEnergy snRes sSteps] = loopyunwrapping(phaseimg);
[msShifts msEnergy msnRes msSteps Mc2s Ms2c] = BPUnwrap(phaseimg);
[usShifts usEnergy usnRes usSteps uStates] = BPStateUnwrap(phaseimg);
[ppShifts ppEnergy ppnRes ppSteps ppStates] = BPStateUnwrapPMin(phaseimg);
[pqShifts pqEnergy pqnRes pqSteps pqStates] = BPStateUnwrapQMin(phaseimg);
[lpqShifts lpqEnergy lpqnRes lpqSteps lpqStates] = BPStateUnwrapQMinLog(phaseimg);
[lpqShifts lpqEnergy lpqnRes lpqSteps lpqStates] = BPStateUnwrapQMinLog8connect(phaseimg);
[lpqShifts lpqEnergy lpqnRes lpqSteps lpqStates] = BPStateUnwrapQMinLog8connectFast(phaseimg);
[pqsShifts pqsEnergy pqsnRes pqsSteps pqsStates] = BPStateUnwrapQsum(phaseimg);

if (msnRes(end) == 0)
    munwrap = PU2D(msShifts.v(:,:,end), msShifts.h(:,:,end));
    figure;seq(munwrap);
end

if (usnRes(end) == 0)
    uunwrap = PU2D(usShifts.v(:,:,end), usShifts.h(:,:,end));
    figure;seq(uunwrap);
end

if (ppnRes(end) == 0)
    ppunwrap = PU2D(ppShifts.v(:,:,end), ppShifts.h(:,:,end));
    figure;seq(ppunwrap);
end

if (pqnRes(end) == 0)
    pqunwrap = PU2D(pqShifts.v(:,:,end), pqShifts.h(:,:,end));
    figure;seq(pqunwrap);
end

if (pqsnRes(end) == 0)
    pqunwrap = PU2D(pqsShifts.v(:,:,end), pqsShifts.h(:,:,end));
    figure;seq(pqunwrap);
end
%% DUCK SET UP
duck = double(imread('duck.png'))/255;
minduck = min(duck(duck>0));
duck(duck>0) = duck(duck>0)-minduck+.03;
duck = duck*3;
% subsample
scale = 5;
duckGT = duck(1:scale:end,1:scale:end);
duckGT2 = duckGT;
duckGT1 = duckGT;
duckGT1(duckGT==0) = 1;
duckGT2(duckGT==0) = 2;
dshiftGT = Phase2PhaseShifts(duckGT);
dwrap = mod(duckGT,1);

phaseimg = dwrap;
sol0 = duckGT-dwrap;
ducky = duckGT;
ducky(duckGT==0)=1;
sol1 = ducky-dwrap;
ducky(duckGT==0)=2;
sol2 = ducky-dwrap;
sol3 = sol2+1;
sol4 = sol0+1;
sol5=sol1+1;


tmp = Phase2PhaseShifts(dwrap);
dsol.h = dshiftGT.h - tmp.h;
dsol.v = dshiftGT.v - tmp.v;
% [Y X] = size(dwrap);
% figure;surf(1:X,1:Y,duckGT);
figure;seq(duckGT)
figure;seq(dwrap);
phaseimg = dwrap;

%% DUCK RUN
[sShifts sEnergy snRes sSteps] = loopyunwrapping(phaseimg);
[msShifts msEnergy msnRes msSteps Mc2s Ms2c] = BPUnwrap(phaseimg);
[usShifts usEnergy usnRes usSteps uStates] = BPStateUnwrap(phaseimg);
[ppShifts ppEnergy ppnRes ppSteps ppStates] = BPStateUnwrapPMin(phaseimg);
[pqShifts pqEnergy pqnRes pqSteps pqStates] = BPStateUnwrapQMin(phaseimg);
[pqShifts pqEnergy pqnRes pqSteps pqlStates] = BPStateUnwrapQMinLog(phaseimg);
[pqShifts pqEnergy pqnRes pqSteps pql8States] = BPStateUnwrapQMinLog8connect(phaseimg);
[pqsShifts pqsEnergy pqsnRes pqsSteps pqsStates] = BPStateUnwrapQsum(phaseimg);

if (msnRes(end) == 0)
    munwrap = PU2D(msShifts.v(:,:,end), msShifts.h(:,:,end));
    figure;seq(munwrap);
end

if (usnRes(end) == 0)
    uunwrap = PU2D(usShifts.v(:,:,end), usShifts.h(:,:,end));
    figure;seq(uunwrap);
end

if (ppnRes(end) == 0)
    ppunwrap = PU2D(ppShifts.v(:,:,end), ppShifts.h(:,:,end));
    figure;seq(ppunwrap);
end

if (pqnRes(end) == 0)
    pqunwrap = PU2D(pqShifts.v(:,:,end), pqShifts.h(:,:,end));
    figure;seq(pqunwrap);
end

if (pqsnRes(end) == 0)
    pqunwrap = PU2D(pqsShifts.v(:,:,end), pqsShifts.h(:,:,end));
    figure;seq(pqunwrap);
end

%%
% test 8connect vs 8connectFast
tic;
[pqShifts pqEnergy pqnRes pqSteps pql8States] = BPStateUnwrapQMinLog8connect(phaseimg,.2,2);
toc;
tic;
[pqfShifts pqfEnergy pqfnRes pqfSteps pql8fStates] = BPStateUnwrapQMinLog8connectFast(phaseimg,.2,2);
toc;
%% Runway image
rw = imread('DepthMap.jpg');
rw = rgb2gray(double(rw)/255);
rw = rw*3;

rwGT = rw(1:5:end,1:5:end);
delphiGT = Phase2PhaseShifts(phiGT);
rwwrap = mod(duckGT,1);
tmp = Phase2PhaseShifts(dwrap);
dsol.h = dshiftGT.h - tmp.h;
dsol.v = dshiftGT.v - tmp.v;
% [Y X] = size(dwrap);
% figure;surf(1:X,1:Y,duckGT);
figure;seq(duckGT)
figure;seq(dwrap);
[dShifts dEnergy dnRes dSteps Mc2s Ms2c] = BPUnwrap(dwrap);
if (dnRes(end) == 0)
    dunwrap = PU2D(dShifts.v(:,:,end), dShifts.h(:,:,end));
    figure;seq(dunwrap);
end

%% Debug shortcuts

% final energy
vE = (dshiftGT.v).^2./(2*sigma);
hE = (dshiftGT.h).^2./(2*sigma);
EnergyGT = sum(vE(:))+sum(hE(:));

% energy of proposed solution
vE = (dShifts.v(:,:,end)).^2./(2*sigma);
hE = (dShifts.h(:,:,end)).^2./(2*sigma);
EnergyGT = sum(vE(:))+sum(hE(:));

% initial energy
vE = (vertshift).^2./(2*sigma);
hE = (horshift).^2./(2*sigma);
EnergyInit = sum(vE(:))+sum(hE(:));

%% Fill in BrightnessConst and fix up Duck test

BrightnessConst = load(fullfile('C:\data\BrightnessCalSVC','collected_stats.mat'),'Clse');
BrightnessConst = BrightnessConst.Clse;
% Fill in BrightnessConst
bc = BrightnessConst;
bc(:,1)=bc(:,2);bc(:,319)=bc(:,318);bc(:,320)=bc(:,318);
bc(2,:)=bc(3,:);bc(1,:)=bc(2,:);bc(199,:)=bc(198,:);bc(200,:)=bc(198,:);
newbc = zeros(size(bc));
for col = 1:N
    for row = ceil(M/2):-1:1
        if (bc(row,col)~=0)
            bval = bc(row,col);
        else
            newbc(row,col) = bval/2+newbc(row,col);
        end
    end
    for row = ceil(M/2):M
        if (bc(row,col)~=0)
            bval = bc(row,col);
        else
            newbc(row,col) = bval/2+newbc(row,col);
        end
    end
end

for row = 1:M
    for col = 165:-1:1
        if (bc(row,col)~=0)
            bval = bc(row,col);
        else
            newbc(row,col) = bval/2+newbc(row,col);
        end
    end
    for col = 165:N
        if (bc(row,col)~=0)
            bval = bc(row,col);
        else
            newbc(row,col) = bval/2+newbc(row,col);
        end
    end
end
bc  =bc+newbc;
BrightnessConst = bc;
clear col row bval bc newbc


if ~exist('ducksynth.mat');
bg = 5; fg = 4; stretch = 3.3;
duck = double(imread('C:\Users\rcrabb\Documents\MATLAB\depthimages\depthDuck.png'));
duckresize = imresize(duck,.64);
duckresize = duckresize(1:200,:)-100;
duckresize(duckresize<=74) = 0;
layers(:,:,1) = duckresize >74;
layers(:,:,2) = duckresize <= 74;
duckresize = fg-duckresize.*stretch/ max(duckresize(layers(:,:,1)));
duckresize(layers(:,:,2)) = bg;
phaseimg = duckresize./A(:,:,3);
nStates = ceil(max(phaseimg(:)));
albedo = layers(:,:,1)*.9;
albedo(layers(:,:,2)) = .6;
[AB] = synthAB(phaseimg*phase2dist,albedo,BrightnessConst,A,'hack');
ABtmp = AB;
duckback = bg*ones(size(AB))./A(:,:,3);
[ABback] = synthAB(duckback*phase2dist,.6*ones(size(AB)),BrightnessConst,A,'hack');
 se = [ 0 0 0 1 0 0 0; 0 1 1 1 1 1 0; 0 1 1 1 1 1 0; 1 1 1 1 1 1 1; 0 1 1 1 1 1 0; 0 1 1 1 1 1 0; 0 0 0 1 0 0 0];
%se = strel('disk',2);%[ 0 1 1 0; 1 1 1 1; 1 1 1 1; 0 1 1 0];%[0 1 0; 1 1 1; 0 1 0];
layers(:,:,1) = imerode(layers(:,:,1),se);
layers(:,:,2) = ~layers(:,:,1);
phaseimg(layers(:,:,2)) = duckback(layers(:,:,2));
AB(layers(:,:,2)) = ABback(layers(:,:,2));
phaseimgwrap = mod(phaseimg,1);
save('ducksynth.mat','AB','phaseimg','phaseimgwrap');
else
    load('ducksynth.mat');
end