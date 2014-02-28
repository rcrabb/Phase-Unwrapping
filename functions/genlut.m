%% get pixel angle lookup table 
function [pixth, M] = genlut(imgdim)
%pixth=camera_pixel_angles(Nr,Nc,fc,cc,kc,alpha_c,dist_flag)

Nr = imgdim(1);
Nc = imgdim(2);


% RGB Intrinsic Parameters
fx_rgb = 5.1885790117450188e+02;
fy_rgb = 5.1946961112127485e+02;
cx_rgb = 3.2558244941119034e+02;
cy_rgb = 2.5373616633400465e+02;

% RGB Distortion Parameters
k1_rgb =  2.0796615318809061e-01;
k2_rgb = -5.8613825163911781e-01;
p1_rgb = 7.2231363135888329e-04;
p2_rgb = 1.0479627195765181e-03;
k3_rgb = 4.9856986684705107e-01;

fc = [fx_rgb, fy_rgb];
cc = [cx_rgb, cy_rgb];
kc = [k1_rgb, k2_rgb, p1_rgb, p2_rgb, k3_rgb];
alpha_c = 0;
dist_flag = 1;
[pixth, M]=camera_pixel_angles(Nr,Nc,fc,cc,kc,alpha_c,dist_flag);