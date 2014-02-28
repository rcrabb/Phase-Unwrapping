function [hists, badhists snaps disthists] = BrightnessExp3(basedir)

usedistance = 1;
savework = 1;
showsnaps = 0;
%basedir = 'C:\data\BrightnessPriorData';
subdirs = dir(basedir);
subdirs = subdirs(4:end);

nsubplotcol = ceil(length(subdirs)/3);

badframepct = .05;
maxAB = 2730;
hists = zeros(maxAB+1,length(subdirs));
badhists = zeros(maxAB+1,length(subdirs));
maxD = 10;
distbins = 201;
Dedges = 1/(distbins-1):maxD/(distbins-1):maxD;

disthists = zeros(200,length(subdirs));
confthresh = 20;
figure;

for i = 1:length(subdirs);
    
    fname = fullfile(basedir,subdirs(i).name,'ConfidenceImage.pgm');
    abmatname = fullfile(basedir,subdirs(i).name,'Brightness.mat');
    distancename= fullfile(basedir,subdirs(i).name,'Distance.mat');
    
    if usedistance
    if ~exist(distancename)
        X = loadpgm(fullfile(basedir,subdirs(i).name,'XImg.pgm'));
        Y = loadpgm(fullfile(basedir,subdirs(i).name,'YImg.pgm'));
        Z = loadpgm(fullfile(basedir,subdirs(i).name,'ZImg.pgm'));
        R = sqrt(X.^2 + Y.^2 + Z.^2)/1000;
        clear X Y Z
        save(distancename,'R');
    else
        load(distancename);
    end
    end
    
    
    if ~exist(abmatname)
        AB = loadpgm(fname);
        save(abmatname,'AB');
    else
        load(abmatname);
    end
    [M N F] = size(AB);
    snaps(:,:,i) = AB(:,:,ceil(F/2));
    % check for frames with too many saturated pixels
    for f = 1:F
        curhist = hist(vec(AB(:,:,f)),maxAB+1)';
        if (curhist(maxAB+1) > badframepct*N*M)
            badhists(:,i) = badhists(:,i) + curhist;
        else
            hists(:,i) = hists(:,i) + curhist;
        end
    end
    
    if (usedistance)
        disthists(:,i) = histc(R(AB>confthresh & AB < maxAB-confthresh),Dedges);
    end
    
    subplot(3,nsubplotcol,i);
    bar(hists(:,i));
end

if savework
if usedistance
    save(fullfile(basedir,'a_histograms.mat'),'hists','badhists','snaps','Dedges','disthists');
else
    save(fullfile(basedir,'a_histograms.mat'),'hists','badhists','snaps');
end
end

if showsnaps
figure;
for i = 1:size(snaps,3)
    subplot(3,nsubplotcol,i);
    imshow(snaps(:,:,i)/255);
end
end

if usedistance
    figure;
    for i = 1:size(disthists,2)
    subplot(3,nsubplotcol,i);
    bar(Dedges,disthists(:,i))
    end
end
return
%%
 distsums = sum(disthists,1);
 disthist = zeros(200,1);
for i = 1:size(disthists,2)
disthist = disthist + disthists(:,i)/distsums(i);
end
figure;bar(Dedges,disthist)

absums = sum(hists,1);
 abhist = zeros(maxAB+1,1);
for i = 1:size(disthists,2)
abhist = abhist + hists(:,i)/absums(i);
end
figure;bar(abhist)


