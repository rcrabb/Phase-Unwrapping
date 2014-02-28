function [Shifts] = Phase2PhaseShifts(phaseimg)
% [Shifts] = Phase2PhaseShifts(phaseimg)
% 
Shifts.v = phaseimg(2:end,:)-phaseimg(1:end-1,:);
Shifts.h = phaseimg(:,2:end)-phaseimg(:,1:end-1);