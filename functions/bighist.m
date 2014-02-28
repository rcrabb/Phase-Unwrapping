function bhist = bighist(datafile,numbins)
if length(numbins) == 1
    numbins(2) = numbins(1);
end

% Or try using all data

%load('workspace11202012.mat','cbrights','dists');

load(datafile,'cbrights','dists');

bhist = zeros(numbins)';

% create a 2D Histogram
maxD = 12;
maxCB = 2.5;
dbins = linspace(0,maxD,numbins(1));
bbins = linspace(0,maxCB,numbins(2));

datamax = length(dists);
chunksize = 1e6;
nsplit = ceil(datamax/chunksize);
for i = 0:nsplit-1
bhist = bhist + hist2(...
    dists(i*chunksize+1:min((i+1)*chunksize,datamax))/1000,...
    cbrights(i*chunksize+1:min((i+1)*chunksize,datamax))*1e7,...
    dbins,bbins);
end