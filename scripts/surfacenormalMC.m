% Monte Carlo simulation for surface normals
%% Generate points uniformly randomly on a hemisphere
ntrials = 50000;
normals = randn(ntrials,3);
normals(:,3) = abs(normals(:,3));
normals = normals./repmat((sum(normals.^2,2).^.5),[1 3]);
pov = [0;0;1];

% TH = pi*(rand(1,ntrials)-.5);
% PHI = asin(2*rand(1,ntrials)-1);
% [normals(:,1),normals(:,2),normals(:,3)] = sph2cart(TH,PHI,1);
% pov = [1;0;0];
% 
% % are these points evenly distributed on the hemisphere?
% figure;plot3(normals(:,1),normals(:,2),normals(:,3),'.');
% axis equal; title('Uniformly distributed points on a hemisphere')

% find the angle between surface normals and illumination source
costheta = normals*pov;
theta = acos(costheta);

% take a look at the distribution of angles
figure;hist(theta,100); 
title('Histogram of surface normal angle');
xlabel('Normal angle (rad)');

% place into bins with weights proportional to projected size
nbins = 100;
bins = zeros(nbins,1);
binsnoweight = bins;
% populate bins
for t = 1:ntrials
    binnum = ceil(theta(t)/(.5*pi) * nbins);
    bins(binnum) = bins(binnum) + costheta(t);    
    binsnoweight(binnum) = binsnoweight(binnum) + 1;
end

% take a look at the histogram
binsize = pi/(2*nbins);
angles = binsize:binsize:pi/2;
figure;plot(angles,bins);
% compare it with cosine(theta) * sine(theta)
hold;
plot(angles,.5*pi*sin(angles).*cos(angles)*ntrials/nbins,'r');
title('Histogram of surface normal angle weighted by projected size');
xlabel('Normal angle (rad)');
legend('histogram bins','graph of sin(\beta)cos(\beta)');
%clear ntrials TH PHI pov costheta theta nbins bins binsize binnum normals


figure;plot(angles,binsnoweight);
hold;
plot(angles,.5*pi*sin(angles)*ntrials/nbins,'r');
title('Histogram of surface normal angle from uniformly populated hemisphere');
xlabel('Normal angle (rad)');
legend('histogram bins','graph of sin(\beta)');

% Monte Carlo simulation for surface normals
%% Generate points uniformly randomly in angle, normal, distance
ntrials = 50000;
maxDistance = 12;
maxBeta = pi/2;
maxAlbedo = 1;
minhistsize = 500;
histstep = .25;
histbins = 0:histstep:maxDistance;

sample = rand(ntrials,3).*repmat([maxDistance maxBeta maxAlbedo], [ntrials 1]);
I = 14.7;

% compute brightness
sample(:,4) = I*sample(:,2).*sample(:,3)./(sample(:,1).^2);

% [s i] = sort(sample(:,4));
% sample = sample(i,:);

% Range of values for B is [0 2.5]
Bs = [ .1 .2 .5 1 1.5 2.5 ];
figure('Color',[1.0 1.0 1.0]);
for iB = 1:length(Bs)
B = Bs(iB);

% increase delta until we have enough samples
delta = .0001;
histidx = sample(:,4)>B-delta & sample(:,4)<B+delta;
histsize = sum(histidx);
while (histsize<minhistsize)
    delta = delta*10;
    histidx = sample(:,4)>B-delta & sample(:,4)<B+delta;
    histsize = sum(histidx);
    
end

subplot(2,3,iB);hist(sample(histidx,1),histbins);xlim([0 maxDistance]);ylim([0 ntrials/100]);
xlabel('D (m)');
%ylabel('P(D|B) = P(B|D) / \int_D P(B|D)');

sTitle = ['Histogram of Distances for Brightness val ' num2str(B)];
title(sTitle);
end

%%
