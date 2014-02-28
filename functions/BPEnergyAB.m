function [energy emap] = BPEnergyAB(phaseimg,sigma,measure,truncate,is8connect,probXs,nStates,dataterm,dtweight)

[M N K] = size(dataterm);

if nargin<9
    dtweight = 1;
end

wrapstate = vec(ceil(phaseimg));
dtv = reshape(dataterm, [M*N K]);
dmap = reshape(dtv(wrapstate),[M N])*dtweight;

[energy emap] = BPEnergy(phaseimg,sigma,measure,truncate,is8connect,probXs,nStates);

emap(5+4*is8connect,:,:) = dmap;

energy = energy + sum(dmap(:))/(M*N);