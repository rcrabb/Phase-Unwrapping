load('workspace20130322','BrightnessConst');
A = gethlut();

rhead = .3;
centerhead = [0 -.8 3];

centerbody = [0 0 3];
sidelengthbody = [.8 1 .5];
centerarm = [0 -.25 3];
sidelengtharm = [2 .15 .3];
centerleg1 = [.25 .9 3];
sidelengthleg1 = [.2 .8 .3];
centerleg2 = [-.25 .9 3];
sidelengthleg2 = [.2 .8 .3];


centerroom = [0 0 0];
sidelengthroom = [10 2.6 10];

bgdist = 5;
albedofg = .8;
albedobg = .6;

[AB1 R1 fgmask1] = synthABsphere(BrightnessConst,rhead,centerhead,A,albedofg, bgdist, albedobg);

[AB2 R2 fgmask2] = synthABbox(BrightnessConst,sidelengthbody,centerbody,A,albedofg, bgdist, albedobg);
[AB R fgmask] = combinesynth(AB1,AB2,R1,R2,fgmask1,fgmask2);

[AB2 R2 fgmask2] = synthABbox(BrightnessConst,sidelengtharm,centerarm,A,albedofg, bgdist, albedobg);
[AB R fgmask] = combinesynth(AB,AB2,R,R2,fgmask,fgmask2);

[AB2 R2 fgmask2] = synthABbox(BrightnessConst,sidelengthleg1,centerleg1,A,albedofg, bgdist, albedobg);
[AB R fgmask] = combinesynth(AB,AB2,R,R2,fgmask,fgmask2);

[AB2 R2 fgmask2] = synthABbox(BrightnessConst,sidelengthleg2,centerleg2,A,albedofg, bgdist, albedobg);
[AB R fgmask] = combinesynth(AB,AB2,R,R2,fgmask,fgmask2);
% 
% [AB2 R2 fgmask2] = synthABbox(BrightnessConst,sidelengthroom,centerroom,A,albedofg, bgdist, albedobg);
% [AB R fgmask] = combinesynth(AB,AB2,R,R2,fgmask,fgmask2);