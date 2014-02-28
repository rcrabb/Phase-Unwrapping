% ducktest

%% Create Duck Test Image
load('BrightnessConst20130222.mat');
A = gethlut();
% mm per 1 phase
phase2dist = 2;
dist2phase = 1/phase2dist;
sigmas = [.01 .1 1 10 ];
dtweights = [0 .0001 .001 .01 .1 1 ];
measures = [{@L1} {@L2}];
m = 1;
trunc = .5;%[ 0 .01 1 ];
truncate = trunc(1);
nConnect = 8;
probXs = 1;
afg = .855; abg = .56;
bg = 3.5; fg = 3; stretch = 1.5;
duck = double(imread('C:\Users\rcrabb\Documents\MATLAB\depthimages\depthDuck.png'));
duckresize = imresize(duck,.64,'nearest');
duckresize = duckresize(1:200,:)-100;
duckresize(duckresize<=0) = 0;
%load('depthDuck.albedo.mat');
layers(:,:,1) = duckresize >0;
layers(:,:,2) = duckresize <= 0;
duckresize = fg-duckresize.*stretch/ max(duckresize(layers(:,:,1)));
duckresize(layers(:,:,2)) = bg;
phaseimg = duckresize./A(:,:,3);
nStates = ceil(max(phaseimg(:)));
albedo = layers(:,:,1)*afg;
albedo(layers(:,:,2)) = abg;
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

%clear se ABtmp layers duckback ABback


%%

[dataterm] = FindBrightnessDataTerm(phaseimgwrap,AB,BrightnessConst,nStates,phase2dist);

eta = .001;
dataterm = -log((dataterm+eta)./repmat(sum(dataterm+eta,3),[1 1 nStates+1]));
% dataterm = zeros(size(dataterm));
% to save time let's reduce the size
rfactor = .5;
piw = imresize(phaseimgwrap,rfactor,'nearest');dt=imresize(dataterm,rfactor,'nearest');
pi =  imresize(phaseimg,rfactor,'nearest');gt = pi-piw;
[tmp dtgt] = min(dt,[],3);

totaltests = length(measures)*length(trunc)*length(sigmas)*length(dtweights)/2;
h = waitbar(0,'Running loopy belief propagation');
% %change this if you want to start over%
if 1
    iter = 0;
    clear endstates s
else
    iter = length(s);
end
    tic
for dtweight = dtweights
for m = 1:1%length(measures)
for truncate = trunc
for sigma = sigmas
    iter = iter+1;   
   
    [WrapState Energy Steps Shifts nRes ] = BPUnwrap(piw,sigma,nStates,measures{m},truncate,nConnect,probXs,dt*dtweight);

    s(iter).sigma = sigma;
    s(iter).nStates = nStates;
    s(iter).Energy = Energy;
 
    s(iter).energy = BPEnergy(piw+WrapState(:,:,end),sigma,measures{m},truncate,0,probXs,nStates);
   % s(iter).gtEnergy = bpe;
    s(iter).nRes = nRes;
%    s(iter).States = uint8(WrapState);
    s(iter).dtflag = dtweight>0;
    s(iter).m = measures{m};
    s(iter).trunc = truncate;
    s(iter).dtweight = dtweight;
    
    
   endstates(:,:,iter) =  uint8(WrapState(:,:,end));
    s(iter).energyAB = BPEnergyAB(piw+WrapState(:,:,end),sigma,measures{m},truncate,0,probXs,nStates,dt*dtweight);    
    ws = WrapState;
    s(iter).pctcorrect = sum(vec(ws(:,:,end)==gt+1))/length(vec(ws(:,:,end)));
    
    waitbar(iter/totaltests,h);
    toc
    if mod(iter,4)==0
save workspace20130425_looping
    end
end
end
end
end
close(h);
save workspace20130321
%% Find cases where AB prior has improved results
bii = 0;
%ws = load workspace20130506;
for iter = 1:length(ws.s)
    % look for cases without prior term
    if (ws.s(iter).dtweight == 0)
        si = find( (logical([ws.s(:).trunc] == ws.s(iter).trunc) & ...
                    logical([ws.s(:).sigma] == ws.s(iter).sigma) )');
        bi = 0;
        for i = 1:length(si)
            if isequal(ws.s(iter).m,ws.s(si(i)).m)
                improv = ws.s(si(i)).percentcorrect - ws.s(iter).percentcorrect;
                if improv > 0
                    bi = improv;
                    bii = si(i);
                end
            end
        end
        if bi > 0
            disp(['Improvement of ' num2str(bi) ', from ' num2str(ws.s(iter).percentcorrect) ' to ' num2str(ws.s(bii).percentcorrect) ...
                ', with sigma: ' num2str(ws.s(bii).sigma) ', dtweight: ' num2str(ws.s(bii).dtweight) ', trunc: ' num2str(ws.s(bii).trunc) ', measure: ' func2str(ws.s(bii).m)]);
        end
    end
end
%% Find cases where AB prior has improved results, from decent results
bii = 0;
%ws = load workspace20130506;
dtweights = [.001 .01 .03 .05 .1 1];
%iter = length(s);
for i = 1:length(s)
    % look for cases without prior term
    if (s(i).pctcorrect > .3 && ~s(i).dtflag)
        bi = 0;
        bii = 0;
        for dtweight = dtweights
            iter = iter+1; 
%         % First generate new tests based on decent control settings
%             [WrapState Energy Steps Shifts nRes ] = BPUnwrap(piw,s(i).sigma,nStates,s(i).m,s(i).trunc,nConnect,probXs,dt*dtweight);
% 
%             s(iter).sigma = s(i).sigma;
%             s(iter).nStates = nStates;
%             s(iter).Energy = Energy;
% 
%             s(iter).energy = BPEnergy(piw+WrapState(:,:,end),s(i).sigma,s(i).m,s(i).trunc,0,probXs,nStates);
%             s(iter).nRes = nRes;
%             %s(iter).States = uint8(WrapState);
%             s(iter).dtflag = dtweight>0;
%             s(iter).m = s(i).m;
%             s(iter).trunc = s(i).trunc;
%             s(iter).dtweight = dtweight;
% 
% 
%             endstates(:,:,iter) = WrapState(:,:,end);
%             s(iter).energyAB = BPEnergyAB(piw+double(WrapState(:,:,end)),s(i).sigma,s(i).m,s(i).trunc,0,probXs,nStates,dt*dtweight);    
%             es = WrapState(:,:,end);
%             s(iter).pctcorrect = sum(vec(es==gt+1))/length(vec(es));

%             % next let's find the Biggest Improvement
%             improv = s(iter).pctcorrect - s(i).pctcorrect;
%             if improv > 0
%                 bi = improv;
%                 bii = iter;
%             end
            % Find all cases better than data term alone
            if s(iter).pctcorrect - gtdt > .01 && s(iter).sigma == 1
                bi = s(iter).pctcorrect - gtdt; bii = iter;
                disp(['Improvement of ' num2str(bi) ', from ' num2str(s(i).pctcorrect) ' to ' num2str(s(bii).pctcorrect) ...
                ', with sigma: ' num2str(s(bii).sigma) ', dtweight: ' num2str(s(bii).dtweight) ', trunc: ' num2str(s(bii).trunc) ', measure: ' func2str(s(bii).m)]);
            end
            
        end
%         if bi > 0
%             disp(['Improvement of ' num2str(bi) ', from ' num2str(s(i).pctcorrect) ' to ' num2str(s(bii).pctcorrect) ...
%                 ', with sigma: ' num2str(s(bii).sigma) ', dtweight: ' num2str(s(bii).dtweight) ', trunc: ' num2str(s(bii).trunc) ', measure: ' func2str(s(bii).m)]);
%         end
    end
end
%% Find cases where smoothness+data beats data alone
bii = 0;
[Conf gtWrapStateDT] = min(dt,[],3);
gtcorrect = sum(vec(gtWrapStateDT-1==gt))/length(vec(gt));
visited = zeros([length(ws.s) 1]);
vsig = zeros(size(sigmas));
vtrunc = zeros(size(trunc));
vdt = zeros(size(un));
%ws = load workspace20130506;
for iter = 1:length(ws.s)
    
    improv = ws.s(iter).percentcorrect - gtcorrect;
    if improv > 0
        vsig = vsig + double(sigmas == ws.s(iter).sigma);
        vtrunc = vtrunc + double(trunc == ws.s(iter).trunc);
        vdt = vdt + double(un == ws.s(iter).dtweight);
    end
%     
%     % look for cases with the prior term
%     if (ws.s(iter).dtweight > 0 && ~visited(iter))
%         si = find( (logical([ws.s(:).trunc] == ws.s(iter).trunc) & ...
%                     logical([ws.s(:).sigma] == ws.s(iter).sigma) )');
%         bi = ws.s(iter).percentcorrect - gtcorrect;
%         bii = iter;
%         for i = 1:length(si)
%             if isequal(ws.s(iter).m,ws.s(si(i)).m)
%                 visited(si(i)) = 1;
%                 improv = ws.s(si(i)).percentcorrect - gtcorrect;
%                 if improv > bi
%                     bi = improv;
%                     bii = si(i);
%                 end
%             end
%         end
%         if bi > 0
%             disp(['Improvement of ' num2str(bi) ...
%                 ', with sigma: ' num2str(ws.s(bii).sigma) ', dtweight: ' num2str(ws.s(bii).dtweight) ', trunc: ' num2str(ws.s(bii).trunc) ', measure: ' func2str(ws.s(bii).m)]);
%         end
%     end

    % Look at cases where energy has gone down with 
    %if ws.s(iter).energy > ws.s(iter).energyAB
end
%%

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

%%
 i = 1;
[tmpnrg tmpmap] = BPEnergy(piw+double(endstates(:,:,i)),s(i).sigma,s(i).m,s(i).trunc,0,probXs,nStates);
[tmpnrgAB tmpmapAB] = BPEnergyAB(piw+double(endstates(:,:,i)),s(i).sigma,s(i).m,s(i).trunc,0,probXs,nStates,dt);


%% Double check that pctcorrect is comparing the right values
for iter = 1:length(s) 
    s(iter).percentcorrect = sum(vec(endstates(:,:,iter)==gt+1))/length(vec(endstates(:,:,iter)));
end

%% Try varying the albedos
load('BrightnessConst20130222.mat');
A = gethlut();
% mm per 1 phase
phase2dist = 2;
dist2phase = 1/phase2dist;
sigma = 1;%[1 5 100];
dtweights = [.005 .03 1];% 0];%;[.001 .01 .03 .05 .1 1];
measures = {@L1};%[{@L1} {@L2}];
m = 1;
truncate = .5;%[ 0 .01 1 ];
nConnect = 8;
probXs = 1;
bg = 3.5; fg = 3; stretch = 1.5;
afg = .855; abg = .56;
resizefactor = .64;
duck = double(imread('C:\Users\rcrabb\Documents\MATLAB\depthimages\depthDuck.png'));
duckresize = imresize(duck,resizefactor);
duckresize = duckresize(1:200,:)-100;
duckresize(duckresize<=0) = 0;
load('depthDuck.albedo.mat');
layersresize = imresize(ducklayers,resizefactor,'nearest');
layersresize = layersresize(1:200,:);
layers = layersresize >0;
layers(:,:,2) = layersresize ==0;
duckresize = fg-duckresize.*stretch/ max(duckresize(layers(:,:,1)));
duckresize(layers(:,:,2)) = bg;
phaseimg = duckresize./A(:,:,3);
nStates = ceil(max(phaseimg(:)));
duckback = bg*ones(size(A(:,:,3)))./A(:,:,3);
 se = [ 0 0 0 1 0 0 0; 0 1 1 1 1 1 0; 0 1 1 1 1 1 0; 1 1 1 1 1 1 1; 0 1 1 1 1 1 0; 0 1 1 1 1 1 0; 0 0 0 1 0 0 0];
%se = strel('disk',2);%[ 0 1 1 0; 1 1 1 1; 1 1 1 1; 0 1 1 0];%[0 1 0; 1 1 1; 0 1 0];
layers(:,:,1) = imerode(layers(:,:,1),se);
layers(:,:,2) = ~layers(:,:,1);
phaseimg(layers(:,:,2)) = duckback(layers(:,:,2));
phaseimgwrap = mod(phaseimg,1);
nlayers = max(layersresize(:));

%%
% %change this if you want to start over%
if 0
    iter = 0;
    clear endstates s albedo
end

tic
totaltests = 1;
h = waitbar(0,'Running loopy belief propagation');
for testnum = 1:totaltests
% 
% albedovalues = rand([1 nlayers+1]);
albedo = layers(:,:,2)*albedovalues(1);
for i = 1:nlayers
    albedo(layersresize==i) = albedovalues(i+1);
end

[AB] = synthAB(phaseimg*phase2dist,albedo,BrightnessConst,A,'hack');
[ABback] = synthAB(duckback*phase2dist,albedovalues(1)*ones(size(AB)),BrightnessConst,A,'hack');
AB(layers(:,:,2)) = ABback(layers(:,:,2));

[dataterm] = FindBrightnessDataTerm(phaseimgwrap,AB,BrightnessConst,nStates,phase2dist);

eta = .001;
dataterm = -log((dataterm+eta)./repmat(sum(dataterm+eta,3),[1 1 nStates+1]));
% dataterm = zeros(size(dataterm));
% to save time let's reduce the size
rfactor = .5;
piw = imresize(phaseimgwrap,rfactor,'nearest');dt=imresize(dataterm,rfactor,'nearest');
pi =  imresize(phaseimg,rfactor,'nearest');gt = pi-piw;
[tmp dtgt] = min(dt,[],3);

for dtweight = dtweights
    
    iter = iter+1;   
   
    [WrapState Energy Steps Shifts nRes ] = BPUnwrap(piw,sigma,nStates,measures{m},truncate,nConnect,probXs,dt*dtweight);

    s(iter).sigma = sigma;
    s(iter).nStates = nStates;
    s(iter).Energy = Energy;
 
    s(iter).energy = BPEnergy(piw+WrapState(:,:,end),sigma,measures{m},truncate,0,probXs,nStates);
   % s(iter).gtEnergy = bpe;
    s(iter).nRes = nRes;
    s(iter).States = uint8(WrapState);
    s(iter).dtflag = dtweight>0;
    s(iter).m = measures{m};
    s(iter).trunc = truncate;
    s(iter).dtweight = dtweight;
    s(iter).albvals = albedovalues;
    
    
   endstates(:,:,iter) =  WrapState(:,:,end);
    s(iter).energyAB = BPEnergyAB(piw+WrapState(:,:,end),sigma,measures{m},truncate,0,probXs,nStates,dt*dtweight);    
    ws = WrapState;
    s(iter).pctcorrect = sum(vec(ws(:,:,end)==gt+1))/length(vec(ws(:,:,end)));
    
    waitbar(testnum/totaltests/3,h);
    toc
end
if mod(iter,9)==0
 %   save workspace20130418running.mat s
end
end
close(h);
%save workspace20130418
%% evaluate some of the new albedo combos
% histograms of percent correct for each dataterm weight
for dtw = dtweights
figure;hist([s([s(:).dtweight]==dtw).pctcorrect],25);
title(['DataTerm weight is: ' num2str(dtw)]);
end
% How many times did each dtweight do a decent job (70%+)
for dtw = dtweights
sum([s([s(:).dtweight]==dtw).pctcorrect]>.7)
end
% eh, useless to try to look at a correlation for dtweight and pctcorrect
figure;plot([s(:).dtweight],[s(:).pctcorrect],'r.');

% now let's try to see which albedos get it wrong
lr = imresize(layersresize,.5,'nearest');
wrong = [];
total = [];
for j = 1:length(s)
    albedovalues = s(j).albvals;
    albedo = (lr==0).*albedovalues(1);
    for i = 1:nlayers
        albedo(lr==i) = albedovalues(i+1);
    end
    wrong = [wrong; albedo(gt+1==endstates(:,:,j))];
    total = [total; albedo(:)];
end
% ooops wrong should be right
 figure;hist(wrong,100)
figure;hist(total,100)

%%
load('BrightnessConst20130222.mat');
A = gethlut();
% mm per 1 phase
phase2dist = 2;
dist2phase = 1/phase2dist;
bg = 3.5; fg = 3; stretch = 1.5;
duck = double(imread('C:\Users\rcrabb\Documents\MATLAB\depthimages\depthDuck.png'));
duckresize = imresize(duck,.64);
duckresize = duckresize(1:200,:)-100;
duckresize(duckresize<=74) = 0;
layers(:,:,1) = duckresize >74;
layers(:,:,2) = duckresize <= 74;
 se = [ 0 0 0 1 0 0 0; 0 1 1 1 1 1 0; 0 1 1 1 1 1 0; 1 1 1 1 1 1 1; 0 1 1 1 1 1 0; 0 1 1 1 1 1 0; 0 0 0 1 0 0 0];
%se = strel('disk',2);%[ 0 1 1 0; 1 1 1 1; 1 1 1 1; 0 1 1 0];%[0 1 0; 1 1 1; 0 1 0];
layers(:,:,3) = imerode(layers(:,:,1),se);
layers(:,:,4) = ~layers(:,:,3);
duckresize = fg-duckresize.*stretch/ max(duckresize(layers(:,:,1)));
duckresize(layers(:,:,2)) = bg;
phaseimg = duckresize./A(:,:,3);
nStates = ceil(max(phaseimg(:)));
itera = 0;
afgn = .855:.005:.865;
abgn = .56:.0025:.65;
clear alb pctalb albws albab
for afgi = 1:length(afgn)
    afg = afgn(afgi);
  for abgi = 1:length(abgn)
    abg = abgn(abgi);
    
    itera = itera +1;
    albedo = layers(:,:,1)*afg;
    albedo(layers(:,:,2)) = abg;
    [AB] = synthAB(phaseimg*phase2dist,albedo,BrightnessConst,A,'hack');
    ABtmp = AB;
    duckback = bg*ones(size(AB))./A(:,:,3);
    [ABback] = synthAB(duckback*phase2dist,abg*ones(size(AB)),BrightnessConst,A,'hack');

    phaseimg(layers(:,:,4)) = duckback(layers(:,:,4));
    AB(layers(:,:,4)) = ABback(layers(:,:,4));
    phaseimgwrap = mod(phaseimg,1);

    [dataterm] = FindBrightnessDataTerm(phaseimgwrap,AB,BrightnessConst,nStates,phase2dist);
    dt=imresize(dataterm,rfactor,'nearest');
    [Conf gtWrapStateDT] = max(dt,[],3);

    alb(itera).afg = afg;
    alb(itera).abg = abg;
    alb(itera).WrapState = gtWrapStateDT;
    albws(:,:,itera) = alb(itera).WrapState;
    albab(:,:,itera) = AB;
    alb(itera).percentcorrect = sum(vec(gtWrapStateDT==gt+1))/length(vec(gt));
    pctalb(afgi,abgi) = alb(itera).percentcorrect;
  end
end
%clear se ABtmp layers duckback ABback AB afgn abgn afgi abgi abg afg phaseimgwrap

%% Use a random white-noise type albedo for foreground.
 % not starting from scratch, so use one of the above sections to set it up, then change the
 % albedo
 % %change this if you want to start over%
if 1
    iter = 0;
    clear endstates s albedo
end

tic
 
randdots = imresize( rand(size(duckresize)/4),4,'nearest'); 
albedo = .9*layers(:,:,2) + randdots.*layers(:,:,1);
[AB] = synthAB(phaseimg*phase2dist,albedo,BrightnessConst,A,'hack');
[ABback] = synthAB(duckback*phase2dist,.9*ones(size(AB)),BrightnessConst,A,'hack');
AB(layers(:,:,2)) = ABback(layers(:,:,2));

[dataterm] = FindBrightnessDataTerm(phaseimgwrap,AB,BrightnessConst,nStates,phase2dist);

eta = .001;
dataterm = -log((dataterm+eta)./repmat(sum(dataterm+eta,3),[1 1 nStates+1]));
% dataterm = zeros(size(dataterm));
% to save time let's reduce the size
rfactor = .5;
piw = imresize(phaseimgwrap,rfactor,'nearest');dt=imresize(dataterm,rfactor,'nearest');
pi =  imresize(phaseimg,rfactor,'nearest');gt = pi-piw;
[tmp dtgt] = min(dt,[],3);

% h = waitbar(0,'Running loopy belief propagation');
for dtweight = dtweights
for sigma = sigmas    
    iter = iter+1;   
   
    [WrapState Energy Steps Shifts nRes ] = BPUnwrap(piw,sigma,nStates,measures{m},truncate,nConnect,probXs,dt*dtweight);

    s(iter).sigma = sigma;
    s(iter).nStates = nStates;
    s(iter).Energy = Energy;
 
    s(iter).energy = BPEnergy(piw+WrapState(:,:,end),sigma,measures{m},truncate,0,probXs,nStates);
   % s(iter).gtEnergy = bpe;
    s(iter).nRes = nRes;
    s(iter).States = uint8(WrapState);
    s(iter).dtflag = dtweight>0;
    s(iter).m = measures{m};
    s(iter).trunc = truncate;
    s(iter).dtweight = dtweight;
%    s(iter).albvals = albedovalues;
    
    
   endstates(:,:,iter) =  WrapState(:,:,end);
    s(iter).energyAB = BPEnergyAB(piw+WrapState(:,:,end),sigma,measures{m},truncate,0,probXs,nStates,dt*dtweight);    
    ws = WrapState;
    s(iter).pctcorrect = sum(vec(ws(:,:,end)==gt+1))/length(vec(ws(:,:,end)));
    
    save workspace20130501_looping.mat s 
%     waitbar(iter/length(dtweights),h);
    toc
end
end
% close(h);
 
save workspace20130501.mat s AB dataterm piw pi albedo endstates