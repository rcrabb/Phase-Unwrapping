
x = 100;
y = [ 50 100 160 220 270 ];
y = [ 50 75 100 130 160 190 220 270 ];

y = [  100 160 220  ];
figure;plot(squeeze(D(x,y,:))',squeeze(AB(x,y,:).^-.5)','.');
%figure;plot(squeeze(AB(x,y,:).^-1),'.r');
%figure;plot(squeeze(D(x,y,:).^2),'.r');

%%
d1 = [];
b1 = [];
[C masks Coss norms D R AB] = BrightnessExpMoving(A, 'aset01');
d1 = [d1 squeeze(D(100,:,:))];
b1 = [b1 squeeze(AB(100,:,:))];
[C masks Coss norms D R AB] = BrightnessExpMoving(A, 'aset02');
d1 = [d1 squeeze(D(100,:,:))];
b1 = [b1 squeeze(AB(100,:,:))];
[C masks Coss norms D R AB] = BrightnessExpMoving(A, 'aset03');
d1 = [d1 squeeze(D(100,:,:))];
b1 = [b1 squeeze(AB(100,:,:))];
[C masks Coss norms D R AB] = BrightnessExpMoving(A, 'aset04');
d1 = [d1 squeeze(D(100,:,:))];
b1 = [b1 squeeze(AB(100,:,:))];
[C masks Coss norms D R AB] = BrightnessExpMoving(A, 'aset05');
d1 = [d1 squeeze(D(100,:,:))];
b1 = [b1 squeeze(AB(100,:,:))];
[C masks Coss norms D R AB] = BrightnessExpMoving(A, 'aset06');
d1 = [d1 squeeze(D(100,:,:))];
b1 = [b1 squeeze(AB(100,:,:))];
%%

[Clse Dmse corrs stdD] = BrightnessExpMoving2(basedir,goodframes);

%% November 19, 2012
basedir = 'C:\data\BrightnessProb\dojo2\moving';
minBrightness = 10;

r100 = load('C:\data\BrightnessProb\dojo2\moving\row100.mat');
goodframes = squeeze(r100.d(160,:)>1820);

[Clse Dmse corrs stdD] = LoadBrightnessData(basedir,goodframes);
[M N] = size(Clse);
% find prior distributions of distance, brightness, and calibrated brightness
dists = [];
brights = [];
cbrights = [];
for row = 1:M
    rname = strcat('row',num2str(row,'%03d'),'.mat');
    load(fullfile(basedir,rname),'d','b');
    tmask = repmat(goodframes,[320 1]) & (b>=minBrightness) & repmat(Clse(row,:)'>0,[1 4300]);
    if (sum(tmask(:))==0)
        continue
    end
    dists = [dists; d(tmask)];
    brights = [brights; b(tmask)];
    R = b./repmat(squeeze(Clse(row,:)'),[1 size(d,2)]);
    cbrights = [cbrights; R(tmask)];
end
[N F] = size(d);

load(fullfile(basedir,'goodpixels.mat'),'dists','brights','cbrights');
calbrightvals = [.1 .25 .5 .75 1 1.5]*1e-7;
distrib = BuildDistanceDistribution(basedir,Clse,calbrightvals,goodframes);


%% nov 26

% alternate to using BuildDistanceDistribution() if
% 'dists','brights','cbrights' are available

%load(fullfile(basedir,'goodpixels.mat'),'dists','brights','cbrights');
calbrightvals = [.1 .25 .5 .75 1 1.5]*1e-7;
delta0 = 7.6415e-010;
maxDistance = 0;
minhistsize = 1.2e6;
histsizes = zeros(size(calbrightvals));
deltas = histsizes;

figure('Color',[1.0 1.0 1.0]);

for cbv = 1:length(calbrightvals)
    % min brightness and good frames have been verified above
    delta = delta0;
    while (histsizes(cbv)< minhistsize)
        histmask =   cbrights < calbrightvals(cbv) + delta & ...
                     cbrights > calbrightvals(cbv) - delta;
        histsizes(cbv) = sum(histmask);
        distrib{cbv}.d = dists(histmask)/1000;
        maxDistance = max(maxDistance,max(distrib{cbv}.d));
        delta = delta*1.1;
    end
    deltas(cbv) = delta;
end

maxDistance =  ceil(maxDistance);
histstep = .25;
histbins = 0:histstep:maxDistance;
for cbv = 1:length(calbrightvals)
    % plot
    subplot(2,3,cbv);hist(distrib{cbv}.d,histbins);
    xlim([0 maxDistance]);
    ylim([0 max(histsizes)]);
    %ylim([0 length(distrib{cbv}.d)/(2*length(histbins))]);
    xlabel('D (m)');
    %ylabel('P(D|B) = P(B|D) / \int_D P(B|D)');

    sTitle = ['Histogram of Distances for Calibrated Brightness val ' num2str(calbrightvals(cbv))];
    title(sTitle);
end

%% Dec 6
%% Dec 9

% using a subsample of the total data
%datafile = 'workspace12062012.mat';

% Or try using all data
datafile = 'workspace11202012.mat';
%datafile = 'C:\data\BrightnessPriorData\goodpixels_sample.mat';

%load(datafile,'cbrights','dists');

% histogram  the data
nbbin = 120;
ndbin = 150;
numbins = [ndbin nbbin];
tic
histmat = bighist(datafile,numbins);
toc
%save('bighistogram.mat','histmat');

% sample value of calibration constant I
I = 14.7113;

% create a 2D Histogram
maxD = 12;
maxCB = 2.5;
dbins = linspace(0,maxD,numbins(1));
bbins = linspace(0,maxCB,numbins(2));

%histmat = hist2(dists, cbrights, dbins, hbins);

bhist = sum(histmat,2)/sum(histmat(:));
dhist = sum(histmat,1)'/sum(histmat(:));
ymax = ceil(max([bhist; dhist])*20)/20;

figure('Color',[1.0 1.0 1.0]);
subplot(1,2,1);
bar(bbins',bhist,'histc');
xlabel('Calibrated Brightness');
title('Histogram of Calibrated Brightness Values');
xlim([0 maxCB]);
ylim([0 ymax]);

subplot(1,2,2);
bar(dbins',dhist,'histc');
title('Histogram of Distances');
xlabel('Distance (m)');
xlim([0 maxD]);
ylim([0 ymax]);


figure('Color',[1.0 1.0 1.0]);
seq(flipud(log(histmat'+1)));
%plot(dists,cbrights,'r.')
title('2D Log Histogram of Distance vs Brightness');
ylabel('Distance (m)')
xlabel('Calibrated Brightness');
set(gca,'XTick',[nbbin/5:nbbin/5:nbbin]);
set(gca,'XTickLabel',{'0.5','1','1.5','2','2.5'})
set(gca,'YTick',[0:ndbin/6:ndbin*5/6]);
set(gca,'YTickLabel',{'12','10','8','6','4','2'})


calbrightvals = [.05 .1 .25 .75 1 1.4 1.75 2];
for i = 1:length(calbrightvals)
    cidx(i) = find(bbins>=calbrightvals(i),1,'first');
end
% % 
% % % first make an uncorrected histogram
% % figure('Color',[1.0 1.0 1.0]);
% % for cbv = 1:length(calbrightvals)
% %     % plot
% %     subplot(2,4,cbv);
% %     B = calbrightvals(cbv);
% %     condhist = squeeze(histmat(cidx(cbv),:));
% %     bar(dbins',condhist/sum(condhist),'histc');
% %     xlim([0 maxD]);
% %     ylim([0 .4]);
% %     xlabel('D (m)');
% % 
% %     sTitle = {'Histogram of Distances'; ['for Calibrated Brightness val ' num2str(B)]};
% %     title(sTitle);
% %     
% %     % ADD PLOT OF HYPOTHESIZED DISTRIBUTION
% %     % distribution based on uniform prior for albedo and
% %     % cos(beta)sin(beta) for surface normal Beta
% %     H = @(D)( 2*D.^2./I .* (1 - B.*D.^2./I) );
% %     dplot = dbins([0 dbins(1:end-1)]<=(I/B)^.5);
% %     hold;plot(dplot,H(dplot)/sum(H(dplot)),'r');
% % 
% % end
% % 
% % % now make a corrected histogram
% % figure('Color',[1.0 1.0 1.0]);
% % for cbv = 1:length(calbrightvals)
% %     % plot
% %     subplot(2,4,cbv);
% %     B = calbrightvals(cbv);
% %     chist = histmat(cidx(cbv),:)'./dhist;
% %     chist(isnan(chist)|isinf(chist)) = 0;
% %     bar(dbins',chist/sum(chist),'histc');
% %     xlim([0 maxD]);
% %     ylim([0 .4]);
% %     xlabel('D (m)');
% % 
% %     sTitle = {'IP Weighted Histogram of Distances'; ['for Calibrated Brightness val ' num2str(B)]};
% %     title(sTitle);
% %     
% %     
% %     % ADD PLOT OF HYPOTHESIZED DISTRIBUTION
% %     % distribution based on uniform prior for albedo and
% %     % cos(beta)sin(beta) for surface normal Beta
% %     H = @(D)( 2*D.^2./I .* (1 - B.*D.^2./I) );
% %     dplot = dbins([0 dbins(1:end-1)]<=(I/B)^.5);
% %     hold;plot(dplot,H(dplot)/sum(H(dplot)),'r');
% % 
% % end


% Plot comparison of histogram, IPW histogram and predicted PDF

% now make a corrected histogram
figure('Color',[1.0 1.0 1.0]);
for cbv = 1:length(calbrightvals)
    % plot
    subplot(2,4,cbv);
    B = calbrightvals(cbv);
    H = @(D)( 2*D.^2./I .* (1 - B.*D.^2./I) );
    pidx = [0 dbins(1:end-1)]<=(I/B)^.5;
    dplot = dbins(pidx);
    condhist = squeeze(histmat(cidx(cbv),pidx));
    chist = histmat(cidx(cbv),pidx)'./dhist(pidx);
    chist(isnan(chist)|isinf(chist)) = 0;
    
    plot(...
        dplot,condhist/sum(condhist),'k',...
        dplot,H(dplot)/sum(H(dplot)),'r');
%        dplot,chist/sum(chist),'g',...
    xlim([0 maxD]);
    ylim([0 .4]);
    xlabel('D (m)');
    legend('Histogram of Data','p(\beta) \propto sin(\beta)cos(\beta)');
    %legend('Histogram of Data','IP Weighted Histogram','p(\beta) \propto sin(\beta)cos(\beta)');
    sTitle = {'Comparison of distributions'; ['for Calibrated Brightness val ' num2str(B)]};
    title(sTitle);
    
    
    % ADD PLOT OF HYPOTHESIZED DISTRIBUTION
    % distribution based on uniform prior for albedo and
    % cos(beta)sin(beta) for surface normal Beta
    H = @(D)( 2*D.^2./I .* (1 - B.*D.^2./I) );
    dplot = dbins([0 dbins(1:end-1)]<=(I/B)^.5);
    hold;plot(dplot,H(dplot)/sum(H(dplot)),'r');

end



%% Jan 18, 2013

% Create vector matrix A using script in BrightnessProbNotes
A = gethlut();
[M N ~] = size(A);

% Create values for brightness calibration coefficient using calibration
% data
basedir = 'C:\data\BrightnessCalSVC';
calsubdirs = dir(fullfile(basedir,'caldim*'));

% select foreground
for i = 1:length(calsubdirs)
    BrightnessExpMoving(A, calsubdirs(i).name);
end

% Examine the calibration images. look for AGC effects.
datafiles = dir(fullfile(basedir,'data_*'));
for i = 1:length(datafiles)
    tmp = load(fullfile(basedir,datafiles(i).name),'AB');
    frameMaxAB(:,i) = max(reshape(tmp.AB,[M*N size(tmp.AB,3)]))';
end

for i = 1:length(datafiles)
    tmp = load(fullfile(basedir,datafiles(i).name),'C','AB');
    medC = median(tmp.C,3);
    gpmask = repmat(medC-5e8,[1 1 size(tmp.C,3)]) < tmp.C & ...
             repmat(medC+5e8,[1 1 size(tmp.C,3)]) > tmp.C;
    MC(:,:,i) = medC;
    ngoodpix(:,i) = sum(reshape(gpmask,[M*N size(gpmask,3)]),1)';
end

goodframes = ngoodpix(:)>50000;
% 
% [Clse Dmse corrs stdD] = BrightnessExpMoving2(basedir,goodframes);
% 
% % script to recompute with distance in m instead of mm
% [Clse Dmse corrs stdD] = BrightnessExpMoving2(basedir,goodframes);

% looks like there's an alternative computation
[Clse Dmse corrs stdD] = LoadBrightnessData(basedir,goodframes');

% Let's check a few pixels to plot 
smppts = [100 100; 70 140; 140 100];
albedo = .15;
for nn = 1:size(smppts,1)
    n = smppts(nn,2);
load(fullfile(basedir,strcat('row',num2str(smppts(nn,1),'%03d'),'.mat')),...
    'd','b','cb','mask');
    idxs = logical(squeeze(mask(n,:)) & goodframes');
    adotn = cb(n,idxs);    bright = b(n,idxs);    dist = d(n,idxs);
    X = albedo*adotn./dist.^2;
    figure;plot(bright,X,'g.');
    idxs = ~logical(squeeze(mask(n,:)) & goodframes');
    adotn = cb(n,idxs);    bright = b(n,idxs);    dist = d(n,idxs);
    X = albedo*adotn./dist.^2;
    hold;plot(bright,X,'ro');
%     idxs = ~logical(squeeze(mask(n,:)));
%     adotn = cb(n,idxs);    bright = b(n,idxs);    dist = d(n,idxs);
%     X = albedo*adotn./dist.^2;
%     hold;plot(bright,X,'ro');
%     idxs = ~logical(goodframes');
%     adotn = cb(n,idxs);    bright = b(n,idxs);    dist = d(n,idxs);
%     X = albedo*adotn./dist.^2;
%     plot(bright,X,'kx');
end

load(fullfile(basedir,'collected_stats.mat'));

basedir = 'C:\data\BrightnessPriorData';
minAB = 100;
maxAB = 2700;
badframepct = .05;

% Filter out bad frames and pixels using brightness thresholds
% Compute the calibrated brightness of remaining pixels and save colletion
% of data
BrightnessExp4(basedir,Clse,maxAB,minAB,badframepct);

subsample = .001;
BrightnessExp5(basedir,subsample);

% View which pixels are included in histograms
tmp = load(fullfile(basedir,subdirs(3).name,'Brightness.mat'));
tmp.R = load(fullfile(basedir,subdirs(3).name,'Distance.mat'),'R');
tmp.R=tmp.R.R;
for f = 1:size(tmp.R,3)
    fmask(:,:,f) = ~imdilate(~(tmp.AB(:,:,f)<=maxAB),ones(5)) & tmp.AB(:,:,f) > minAB;
end
figure;zseq(tmp.R.*fmask*700);

calbrightvals = [.1 .25 .5 .75 1 1.5]*1e-7;
distrib = BuildDistanceDistribution(basedir,(Clse/1e6),calbrightvals,goodframes);

%% Feb 18 2013

% check values for brightness const
ABgt = loadpgm('C:\data\BrightnessCalSVC\caldim1\ConfidenceImage.pgm');
% just using a subset of good pixels, up to the 200th column
ABgt = ABgt(:,1:200,1);
A = gethlut();

BrightnessConst = load(fullfile('C:\data\BrightnessCalSVC','collected_stats.mat'),'Clse');
BrightnessConst = BrightnessConst.Clse;

Z =  loadpgm('C:\data\BrightnessCalSVC\caldim1\ZImg.pgm');
depth = Z(:,:,1)./A(:,:,3);
depth = depth(:,1:200);
ABtest = synthAB(depth(:,1:200)/1000,.15*ones(200,200),BrightnessConst(:,1:200),A(:,1:200,:));