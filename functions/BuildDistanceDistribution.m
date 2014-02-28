function [distrib] = BuildDistanceDistribution(basedir,Clse,calbrightvals,goodframes)

% Load in saved data
%[M N] = size(A);
M = 200; N = 320;
%basedir = 'C:\data\BrightnessProb\dojo2\moving';

maxDistance = 12;
minhistsize = 500;
histstep = .25;
histbins = 0:histstep:maxDistance;

albedo = .5;
delta = 7.6415e-010;
minBrightness = 10;

if exist(fullfile(basedir,'distrib.mat'))
    load(fullfile(basedir,'distrib.mat'));
else
    
for cbv = 1:length(calbrightvals)
    distrib{cbv}.d= [];
end
histsizes = zeros(6,1);

for row = 1:M

I = Clse(row,:);
if (sum(I>0)==0)
    continue;
end
rname = strcat('row',num2str(row,'%03d'),'.mat');
load(fullfile(basedir,rname));
[N F] = size(d);

for cbv = 1:length(calbrightvals)
    histmask =  (repmat(goodframes,[N 1]) & ...
                 b > minBrightness) & ...
                 b./repmat(I',[1 F]) < calbrightvals(cbv) + delta & ...
                 b./repmat(I',[1 F]) > calbrightvals(cbv) - delta;
    histsizes(cbv) = histsizes(cbv) + sum(histmask);
    distrib{cbv}.d = [distrib{cbv}.d; d(histmask)/1000];
end

end


if (1)
    save(fullfile(basedir,'distrib.mat'),'distrib','calbrightvals','delta','histsizes');
end

end

figure('Color',[1.0 1.0 1.0]);

for cbv = 1:length(calbrightvals)
    

subplot(2,3,cbv);hist(distrib{cbv}.d,histbins);
xlim([0 maxDistance]);
ylim([0 max(histsizes)]);
%ylim([0 length(distrib{cbv}.d)/(2*length(histbins))]);
xlabel('D (m)');
%ylabel('P(D|B) = P(B|D) / \int_D P(B|D)');

sTitle = ['Histogram of Distances for Calibrated Brightness val ' num2str(calbrightvals(cbv))];
title(sTitle);
end