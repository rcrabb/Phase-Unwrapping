function [s bestsettings] =  findbestparameters(settings, phaseimg, dataterm, s)
    
% Save every X seconds
savetime = 1600;
bestscore = 0;

phase2dist = settings.phase2dist;
nConnects = settings.nConnects;
measures = settings.measures;
trunc = settings.trunc;
sigmas = settings.sigmas;
dtweights = settings.dtweights;
dataname = settings.dataname;
maxwrap = settings.maxwrap;
WrapStateGT = settings.WrapStateGT;
pXss = settings.pXss;
rfactor = settings.rfactor;

filename = [dataname '_' num2str(phase2dist) '_' num2str(rfactor) '_dtnorm' num2str(settings.dtnorm) '_' date '.mat'];

totaltests = length(nConnects)*length(measures)*length(trunc)*length(sigmas)*length(dtweights)*length(pXss);
h = waitbar(0,'Running loopy belief propagation');
% %change this if you want to start over%
if nargin < 4
    iter = 0;
    clear endstates s
else
    iter = length(s);
end

tic
for pXs = pXss
for nConnect = nConnects
for dtweight = dtweights
for m = 1:length(measures)
for truncate = trunc
for sigma = sigmas
    iter = iter+1;   
    
    command = '[WrapState,~,~,~,~] = BPUnwrap(phaseimg,sigma,maxwrap,measures{m},truncate,nConnect,pXs,dataterm*dtweight)';
    T = evalc(command);
    
    s(iter).sigma = sigma;
    s(iter).dtnorm = settings.dtnorm;
    s(iter).pXs = pXs;
    s(iter).maxwrap = maxwrap;
    s(iter).States = uint8(WrapState);
    s(iter).dtflag = dtweight>0;
    s(iter).m = measures{m};
    s(iter).truncate = truncate;
    s(iter).dtweight = dtweight;  
    endstates(:,:,iter) =  uint8(WrapState(:,:,end));
    s(iter).endstate = uint8(WrapState(:,:,end));
    s(iter).pctcorrect = sum(vec(WrapState(:,:,end)==WrapStateGT+1))/length(WrapStateGT(:));
    
    % update timekeeping
    waitbar(iter/totaltests,h);
    telapsed = toc
    if telapsed > savetime
        save(['workspace_ourresults_mid_' filename],'s','endstates','phaseimg','WrapStateGT','dataterm','settings','bestscore','bestsettings');
        tic;
    end
    
    % update best score
    if s(iter).pctcorrect > bestscore
        bestscore = s(iter).pctcorrect;
        bestsettings = s(iter);
    end
end
end
end
end
end
end
close(h);
save(['workspace_ourresults_' filename],'s','endstates','phaseimg','WrapStateGT','dataterm','settings','bestscore','bestsettings');
s = rmfield(s,'States');