function [bprob] = BrightnessProb(PhaseImg, Brightness, nStates, C, Phase2Dist)

[M N] = size(PhaseImg);
bprob = zeros(M,N,nStates);
int_precision = 10000;
if (nargin < 5)
    Phase2Dist = 1;
end

% create or load lookup-table 
if (exist('brightness_integral_table.mat'))
    load('brightness_integral_table.mat');
else
    iter = 0;
    for c = 1/(int_precision*2):1/(int_precision):1
        iter = iter+1;
        F = @(x)(1+(sec(x).*tan(x).*c).^2).^.5;
        Q(iter) = quad(F,0,acos(c));
    end
    save('brightness_integral_table.mat','int_precision','Q');
end

bprob = ((repmat(PhaseImg,[1 1 nStates]) + ...
          repmat(reshape(0:nStates-1,[1 1 nStates]),[M N 1] ) ).* ...
                                                   Phase2Dist).^2 ...
        *repmat(Brightness./C,[1 1 nStates]);
bprob = Q(ceil(bprob*int_precision+.000001)).*(bprob<=1);
bprob = bprob./repmat(sum(bprob,3),[1 1 nStates]);
