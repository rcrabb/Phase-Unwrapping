function [s bestsettings] =  findbestparametersChoi(settings, phaseimg, dataterm, s)
    
% Save every X seconds
savetime = 900;
bestscore = 0;

taus = settings.taus;
gammas = settings.gammas;
betas = settings.betas;
phase2dist = settings.phase2dist;
nConnect = settings.nConnects;
measures = settings.measures;
truncate = 0;
dtweights = settings.dtweights;
dataname = settings.dataname;
maxwrap = settings.maxwrap;
WrapStateGT = settings.WrapStateGT;
rfactor = settings.rfactor;

filename = [dataname '_' num2str(phase2dist) '_' num2str(rfactor) '_' date '.mat'];

%totaltests = length(nConnects)*length(measures)*length(trunc)*length(sigmas)*length(dtweights)*length(pXss);
%h = waitbar(0,'Running loopy belief propagation');
% %change this if you want to start over%
if nargin < 4
    iter = 0;
    clear endstates s
else
    iter = length(s);
end

tic
for dtweight = dtweights
for tau = taus
for gamma =  gammas
for beta = betas
    iter = iter+1;   
    s(iter).tau = tau;
    s(iter).gamma = gamma;
    s(iter).beta = beta;
    s(iter).phase2dist = phase2dist;
    s(iter).maxwrap = maxwrap;
    s(iter).dtterm = [];
    s(iter).sigma = [];
    s(iter).pXs = [];
    s(iter).dtnorm = [];
    command = '[WrapState,~,~,~,~] = BPUnwrap(phaseimg,s(iter),maxwrap,@FindSmoothnessTermChoi,truncate,nConnect,1,dtweight*dataterm);';
    T = evalc(command);
    
    s(iter).maxwrap = maxwrap;
    s(iter).States = uint8(WrapState);
    s(iter).dtflag = dtweight>0;
    s(iter).m = @FindSmoothnessTermChoi;
    s(iter).truncate = truncate;
    s(iter).dtweight = dtweight;  
    endstates(:,:,iter) =  uint8(WrapState(:,:,end));
    s(iter).endstate = uint8(WrapState(:,:,end));
    s(iter).pctcorrect = sum(vec(WrapState(:,:,end)==WrapStateGT+1))/length(WrapStateGT(:));
    
    % update timekeeping
 %   waitbar(iter/totaltests,h);
    telapsed = toc
    if telapsed > savetime
        save(['workspace_choi_mid_' filename],'s','endstates','phaseimg','WrapStateGT','dataterm','settings','bestscore','bestsettings');
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
%close(h);
save(['workspace_choi_' filename],'s','endstates','phaseimg','WrapStateGT','dataterm','settings','bestscore','bestsettings');
s = rmfield(s,'States');