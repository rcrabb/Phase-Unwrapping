function [dists brights cbrights] = BrightnessExp5(basedir,subsample)
% When all of the good pixels have been selected, combine the data from all datasets
% Optionally, combine only a subsample
subdirs = getdir(basedir);
dists = [];
brights = [];
cbrights = [];

if (subsample == 1) || isempty(subsample)
    % When all of the good pixels have been selected, combine the data from all datasets
% Optionally, combine only a subsample
h = waitbar(0,'Loading FULL brightness  data.');
for i = 1:length(subdirs)
    tmp = load(fullfile(basedir,subdirs(i).name,'goodpixels.mat'),'dists','brights','cbrights');
    dists = [dists; tmp.dists(:)];
    brights = [brights; tmp.brights(:)];
    cbrights = [cbrights; tmp.cbrights(:)];
    waitbar(i/length(subdirs));
end
close(h);
save(fullfile(basedir,'goodpixels.mat'),'dists','brights','cbrights');

else % use subsample
    
h = waitbar(0,'Loading  brightness  data.');
for i = 1:length(subdirs)
    tmp = load(fullfile(basedir,subdirs(i).name,'goodpixels.mat'),'dists','brights','cbrights');
    sampmask = randperm(length(tmp.dists));
    sampmask = sampmask(1:floor(subsample*length(sampmask)));
    dists = [dists; tmp.dists(sampmask)];
    brights = [brights; tmp.brights(sampmask)];
    cbrights = [cbrights; tmp.cbrights(sampmask)];
    waitbar(i/length(subdirs));
end
close(h);
save(fullfile(basedir,'goodpixels_sample.mat'),'dists','brights','cbrights');


end