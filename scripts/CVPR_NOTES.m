%% CVPR Experiment %% 
cvprdata = 'C:\data\cvprdata\';
dirlist = getdir(cvprdata);


% decrease image size for faster testing
rfactor = .5;

load('BrightnessConst20130222.mat');

%%
% Load data
bigs = [];
for testnum = 1:length(dirlist)
dataname = dirlist{testnum};
[R AB] = pgm2mat([cvprdata dataname]);
r = median(R,3);
ab = median(AB,3);

clear R AB

%
% implementation parameter settings
for phase2dist = [ 3000]
for dtnorm = [ 1]
abthresh = 50;
maxwrap = ceil(max(r(:))/phase2dist);%
tol = .0001;
maxiter = 500;

pXss = [ 1];
sigmas = [ 1e-3 1e-2 1 ];
dtweights = [0  1e-4 1];
trunc = [ .5];%.5;%[ 0 .01 1 ];h


settings.phase2dist = phase2dist;
settings.nConnects = 8;
settings.measures = {@L1};

settings.dtnorm = dtnorm;
settings.trunc = trunc;
settings.sigmas = sigmas;
settings.dtweights = dtweights;
settings.dataname = dataname;
settings.maxwrap = maxwrap;
settings.pXss = pXss;

phaseimg = mod(r/phase2dist,1);
WrapStateGT = floor(r/phase2dist);
dataterm = FindBrightnessDataTerm(phaseimg,ab,BrightnessConst,maxwrap,phase2dist);
eta = min(dataterm(dataterm~=0))/100;
if dtnorm
    dataterm = -log((dataterm+eta)./repmat(sum(dataterm+eta,3),[1 1 maxwrap+1]));
else
    dataterm = -log((dataterm+eta));%./repmat(sum(dataterm+eta,3),[1 1 maxwrap+1]));
end
% resize inputs
if rfactor < 1
    dataterm=imresize(dataterm,rfactor,'nearest');
    phaseimg =  imresize(phaseimg,rfactor,'nearest');
    WrapStateGT = imresize(WrapStateGT,rfactor,'nearest');    
end

settings.rfactor = rfactor;
settings.WrapStateGT = WrapStateGT;

% Run the loop to test the best combination of parameters
[s bestsettings] = findbestparameters(settings, phaseimg, dataterm);
bigs = [bigs s];
end
end
end


%% REPEAT FOR CHOI
% Load data
choidata = 'C:\data\ChoiTest\';
dirlist = getdir(choidata);

bigs = [];
for testnum = 1:length(dirlist)
dataname = dirlist{testnum};
[R AB] = pgm2mat([choidata dataname]); 
r = median(R,3);
ab = median(AB,3);

clear R AB

%
% implementation parameter settings
for phase2dist = [ 3000]
abthresh = 50;
maxwrap = ceil(max(r(:))/phase2dist);%
tol = .0001;
maxiter = 500;

dtweights = [1.0000e-003 1];
taus =  [0.0100 1 100 10000];%10.^(2:2:8);
gammas = [1.0000e-006 1.0000e-004 0.0100 1 100];%10.^(-4:2:2);
betas =  [1.0000e-003 0.0100 0.1000 1 10 100 1000];


settings.phase2dist = phase2dist;
settings.nConnects = 8;
settings.measures = {@FindSmoothnessTermChoi};

settings.taus = taus;
settings.gammas = gammas;
settings.betas = betas;
settings.dtnorm = 1;
settings.trunc = [];
settings.sigmas = [];
settings.dtweights = dtweights;
settings.dataname = dataname;
settings.maxwrap = maxwrap;
settings.pXss = [];

phaseimg = mod(r/phase2dist,1);
WrapStateGT = floor(r/phase2dist);
dataterm = FindBrightnessDataTermChoi(phaseimg,ab,abthresh,maxwrap,phase2dist,tol,maxiter);

% resize inputs
if rfactor < 1
    dataterm=imresize(dataterm,rfactor,'nearest');
    phaseimg =  imresize(phaseimg,rfactor,'nearest');
    WrapStateGT = imresize(WrapStateGT,rfactor,'nearest');    
end

settings.rfactor = rfactor;
settings.WrapStateGT = WrapStateGT;

% Run the loop to test the best combination of parameters
[s bestsettings] = findbestparametersChoi(settings, phaseimg, dataterm);
bigs = [bigs s];
end
end


%% Test output
r2 = [];
testdir = dir('C:\data\tmpchoi\w*');
dirlist = {testdir(:).name};
for testnum = 1:length(dirlist)
    dataname = ['C:\data\tmpchoi\' dirlist{testnum}];
    %load(['workspace_ourresults_' dataname '_' num2str(phase2dist) '_' num2str(rfactor) '_' date '.mat'],'s');
    load(dataname,'s');
    r2 = [r2 s];
end


%% Find best data terms
cvprdata = 'C:\data\CVPRData\';
dirlist = getdir(cvprdata);
clear pcts
phases = [ 500:500:3500];
for testnum = 1:length(dirlist)
        dataname = dirlist{testnum};
        [R AB] = pgm2mat([cvprdata dataname]);
        r = median(R,3);
        ab = median(AB,3);
        clear R AB

    for j = 1:length(phases)
        phase2dist = phases(j);
        maxwrap = ceil(max(r(:))/phase2dist);%
%         phaseimg = mod(r/phase2dist,1);
%         gt = floor(r/phase2dist);
%         dt = FindBrightnessDataTerm(phaseimg,ab,BrightnessConst,maxwrap,phase2dist);
%         [Conf gtWrapStateDT] = max(dt,[],3);
%         gtcorrect = sum(vec(gtWrapStateDT-1==gt))/length(vec(gt));
%         pcts(testnum,j) = gtcorrect;
         maxw(testnum,j) = maxwrap;
    end
    disp('Hey!')
end





















%% find best choi result
r2 = [];
testdir = dir('C:\data\tmpallchoi\w*');
dirlist = {testdir(:).name};
for testnum = 1:length(dirlist)
    dataname = ['C:\data\tmpallchoi\' dirlist{testnum}];
    %load(['workspace_ourresults_' dataname '_' num2str(phase2dist) '_' num2str(rfactor) '_' date '.mat'],'s');
    load(dataname,'ChoiParams');
    if max([ChoiParams.pctcorrect]) > .9
        disp([num2str(max([ChoiParams.pctcorrect])) ' - ' dirlist{testnum}]);
    end
end






























%% RESULTS ANALASYS %%
%% Test output
r1 = [];
r2 = [];
dirdt = 'C:\data\10_29_visionlab_dt\';
dirst = 'C:\data\10_29_visionlab_st\';
filenames = dir([dirst 'w*']);
dirlist = {filenames(:).name};
solst = [];
solboth = [];
soldt = [];
for testnum = 1:length(dirlist)
    dataname = [dirst dirlist{testnum}];
    load(dataname,'s');
    r1 = [r1 s];
    [bestst(testnum), sti(testnum)] = max([s(:).pctcorrect]);
    solst(:,:,testnum) = s(sti(testnum)).endstate(:,:,end);
    dataname = [dirdt dirlist{testnum}];
    load(dataname,'s','WrapStateGT','dataterm');
    r2 = [r2 s];
    [bestboth(testnum), stb(testnum)] = max([s(:).pctcorrect]);
    solboth(:,:,testnum) = s(stb(testnum)).endstate(:,:,end);
    
    [Conf gtWrapStateDT] = min(dataterm,[],3);
    soldt(:,:,testnum) = gtWrapStateDT;
    bestdt(testnum) = sum(vec(gtWrapStateDT-1==WrapStateGT))/length(vec(WrapStateGT));
    baseline(testnum) = sum(WrapStateGT(:)==0);
end
 
%%  Compare data terms alone Ours to Choi
cvprdata = 'C:\data\CVPRData\';
p2drange = [1000:250:3500];
[pctsOurs pctsChoi] = comparebrightnessterms(cvprdata, p2drange);

xax = 299792458./(2*p2drange*1e3);
xax = xax(end:-1:1);
xla = round(xax([1 3 5 6 7 8 9 10 11]));

pOurs = mean(pctsOurs);
pOurs = pOurs(end:-1:1);
pChoi = mean(pctsChoi);
pChoi = pChoi(end:-1:1);

figure('Color',[1.0 1.0 1.0]);
plot(xax,pOurs,':',xax,pChoi,'--');
xlabel('Modulation Frequency (Mhz)');
ylabel('Percent Correct Pixel Labels');
set(gca,'XTick',xla)
for i = 1:length(xla)
    xlab{i} = num2str(xla(i));
end
set(gca,'XTickLabel',xlab)
legend('Proposed method','Choi et al.');


%% Compare full terms Ours vs Choi
dirchoi = 'C:\data\10_29_visionlab_choi\';
dirours = 'C:\data\10_29_visionlab_dt\';
% dirchoi = 'C:\data\10_31_dojo_choi\';
% dirours = 'C:\data\10_31_dojo_ours\';
[pctsOurs pctsChoi solboth solchoi GTs] = comparefullterms(dirours,dirchoi);

%% Compare full terms Choi vs Choi
dirchoi = 'C:\data\10_31_dojo_choi\';
dirours = 'C:\data\10_31_dojo_choi2\';
% dirchoi = 'C:\data\10_31_dojo_choi\';
% dirours = 'C:\data\10_31_dojo_ours\';
[pctsChoi0 pctsChoi1 solchoi0 solchoi1 GTs] = comparefulltermschoi(dirours,dirchoi);