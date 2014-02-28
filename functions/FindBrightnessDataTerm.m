function [dataterm] = FindBrightnessDataTerm(phaseimg,ab,brightnessConst,maxWraps,phase2dist)
% Phase should be normalized to range from 0 to 1 instead of -pi to pi or 0 to 2pi
numStates = maxWraps+1;
[M N] = size(phaseimg);
eta = 0;
% check phase2dist for units (m vs mm)
if phase2dist > 100
    phase2dist = phase2dist / 1000;
end

% D.^2./I .* (1 - B.*D.^2./I)

% P(WrapState|Phase,Brightness) = (Distance)^2 * (1 - Brightness * Distance^2 / BrightnessConst)
%                          -----------------------------------------------------------------------
%                          \Sum_Distances (Distance)^2 * (1 - Brightness * Distance / BrightnessConst)

Dist = (repmat(phaseimg,[1 1 numStates]) + repmat(reshape(0:numStates-1,[1 1 numStates]),[M N 1]))*phase2dist;
% Recall only valid distances have 1 >= Brightness * Distance^2 / BrightnessConst
mask = Dist.^2 .* repmat(ab ./ brightnessConst, [1 1 numStates]);
dataterm = Dist.^2 .* (1-mask) ./ repmat(brightnessConst, [1 1 numStates]);
mask = mask <= 1;
%dataterm = dataterm.*mask ./repmat(sum(dataterm.*mask,3),[1 1 numStates]);
dataterm = dataterm.*mask + ~mask.*eta;

% single line version
%dataterm = (Dist^2 .* (1-(Dist.^2 .* ab ./ brightnessConst))).*((Dist.^2 .* ab ./ brightnessConst) <= 1) ./repmat(sum((Dist^2 .* (1-(Dist.^2 .* ab ./ brightnessConst))).*((Dist.^2 .* ab ./ brightnessConst) <= 1),3),[1 1 numStates]);
