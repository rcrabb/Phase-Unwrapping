function [pctsOurs pctsChoi] = comparebrightnessterms(cvprdata, p2drange, rfactor)

% Compare data terms alone Ours to Choi
dirlist = getdir(cvprdata);

if nargin < 3
    rfactor = 1;
end

load('BrightnessConst20130222.mat');

%%
% Load data
pctsOurs = [];
pctsChoi = [];
for testnum = 1:length(dirlist)
    dataname = dirlist{testnum};
    [R AB] = pgm2mat([cvprdata dataname]);
    r = median(R,3);
    ab = median(AB,3);

    clear R AB

    %
    % implementation parameter settings
    for p = 1:length(p2drange)
        phase2dist = p2drange(p);
        abthresh = 50;
        maxwrap = ceil(max(r(:))/phase2dist);%
        tol = .0001;
        maxiter = 1500;

        phaseimg = mod(r/phase2dist,1);
        WrapStateGT = floor(r/phase2dist);
        dataterm = FindBrightnessDataTerm(phaseimg,ab,BrightnessConst,maxwrap,phase2dist);
        datatermChoi = FindBrightnessDataTermChoi(phaseimg,ab,abthresh,maxwrap,phase2dist,tol,maxiter);

        % resize inputs
        if rfactor < 1
            dataterm=imresize(dataterm,rfactor,'nearest');
            datatermChoi=imresize(dataterm,rfactor,'nearest');
            WrapStateGT = imresize(WrapStateGT,rfactor,'nearest');    
        end

        [Conf gtWrapStateDT] = max(dataterm,[],3);
        pctsOurs(testnum,p) = sum(vec(gtWrapStateDT-1==WrapStateGT))/length(vec(WrapStateGT));
        
        [Conf gtWrapStateDTChoi] = min(datatermChoi,[],3);
        pctsChoi(testnum,p) = sum(vec(gtWrapStateDTChoi-1==WrapStateGT))/length(vec(WrapStateGT));

    end
    testnum/length(dirlist)
end