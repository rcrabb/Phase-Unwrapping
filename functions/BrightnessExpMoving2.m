function [Clse Dmse corrs stdD] = BrightnessExpMoving2(basedir,goodframes)

% Load in saved data
%[M N] = size(A);
M = 200; N = 320;
%basedir = 'C:\data\BrightnessProb\dojo2\moving';
Clse = zeros(M,N);
Dmse = Clse;
stdD = Clse;
corrs = Clse;

% Check if row files have already been created, and if so skip to loading
if ~exist(fullfile(basedir,'row001.mat'))

for row = 1:M

rname = strcat('row',num2str(row,'%03d'),'.mat');
d = [];
b = [];
cb = [];
mask = [];

 % possible saved data: 'C','masks','Coss','norms','R','AB','subdir','basedir'
dirlist = dir(fullfile(basedir,'data_*'));
for i = 1:length(dirlist)
    load(fullfile(basedir,dirlist(i).name),'R');
    d = [d squeeze(R(row,:,:))];
    clear R;
    load(fullfile(basedir,dirlist(i).name),'AB');
    b = [b squeeze(AB(row,:,:))];
    clear AB;    
    load(fullfile(basedir,dirlist(i).name),'Coss');
    cb = [cb squeeze(Coss(row,:,:))];
    clear Coss;    
    load(fullfile(basedir,dirlist(i).name),'masks');
    mask = [mask squeeze(masks(row,:,:))];
    clear masks;
    
end
[N F] = size(d);

for n = 1:N 
    if (sum(mask(n,:))>300)
    idxs = logical(squeeze(mask(n,:)) & goodframes);
    ns = sum(idxs);
    adotn = cb(n,idxs);
    bright = b(n,idxs);
    dist = d(n,idxs);
    X = sqrt(adotn./bright);
    Clse(row,n) = sum(X.*dist)/sum(X.*X);
    % Compute mean square error
    Dmse(row,n) = mean((dist-Clse(row,n)*X).^2);
    % Compute standard deviation
    stdD(row,n) = sqrt(Dmse(row,n)*ns/(ns-2));
    % Find per frame error
  %  Derr(row,n,idxs) = dist-Clse(row,n)*X;
    correlation = corrcoef(dist,X);
    corrs(row,n) = correlation(2,1);
    else
    Clse(row,n) = 0;
    Dmse(row,n) = 0;
    corrs(row,n) = 0;
    end
end
save(fullfile(basedir,rname),'d','b','cb','mask');

end

% row files have been created, just load and recompute
% FOLLOW UP TO:
% if ~exist(fullfile(basedir,'row001.mat'))
else
    
for row = 1:M

 % possible saved data: 'C','masks','Coss','norms','R','AB','subdir','basedir'
rname = strcat('row',num2str(row,'%03d'),'.mat');
load(fullfile(basedir,rname),'d','b','cb','mask');

% check that distance is in meters
if nanmean(d) > 500
    d = d/1000;
    save(fullfile(basedir,rname),'d','b','cb','mask');
end
if size(goodframes,2) == 1
    goodframes = goodframes';
end

[N F] = size(d);

for n = 1:N 
    if (sum(mask(n,:))>300)
    idxs = logical(squeeze(mask(n,:) & goodframes));
    ns = sum(idxs);
    adotn = cb(n,idxs);
    bright = b(n,idxs);
    dist = d(n,idxs);
    X = sqrt(adotn./bright);
    Clse(row,n) = sum(X.*dist)/sum(X.*X);
    % Compute mean square error
    Dmse(row,n) = mean((dist-Clse(row,n)*X).^2);
    % Compute standard deviation
    stdD(row,n) = sqrt(Dmse(row,n)*ns/(ns-2));
    % Find per frame error
  %  Derr(row,n,idxs) = dist-Clse(row,n)*X;
    correlation = corrcoef(dist,X);
    corrs(row,n) = correlation(2,1);
    else
    Clse(row,n) = 0;
    Dmse(row,n) = 0;
    corrs(row,n) = 0;
    end
end

end
    
end
    
save(fullfile(basedir,'collected_stats.mat'),'corrs','Clse','stdD','Dmse');
