[Bhist Rhist] = surface_discontiniuity();
ab = abs(Bhist);
r = abs(Rhist);
ablog = log(ab+1);
rlog = log(r+1);
abedges = logspace(0,log10(max(ab+1)),64)-1;
redges = logspace(0,log10(max(r+1)),64)-1;
histmat = hist2(r, ab, redges, abedges);
figure; pcolor(redges,abedges,histmat'); colorbar ; axis square tight ;

abedges = linspace(0,8,32);
redges = linspace(0,9,32);
histmat = hist2(rlog, ablog, redges, abedges);
figure('Color',[1.0 1.0 1.0]); pcolor(redges,abedges,histmat'); colorbar ; axis square tight ;



%%
load('C:\data\nyu_depth_v2_labeled.mat');
clear rawDepths

rf = 4;
maxFrames = 200;
fname = ['C:\data\nyu_mini_' num2str(rf) '_' num2str(maxFrames) '.mat'];
depths = depths(1:rf:end,1:rf:end,1:maxFrames);
instances = instances(1:rf:end,1:rf:end,1:maxFrames);
labels = labels(1:rf:end,1:rf:end,1:maxFrames);
images = images(1:rf:end,1:rf:end,:,1:maxFrames);
fname = ['C:\data\nyu_mini_' num2str(rf) '_' num2str(maxFrames) '.mat'];
save(fname,'depths','images','labels','instances');

%% try detecting edges


edges = single(zeros(size(images,1),size(images,2),size(images,4)));
 for i = 1:size(images,4)
     edges(:,:,i) = edgesDetect(images(:,:,:,i),model);
 end
 save(fname,'depths','images','labels','instances','edges');
 
 %% split up the NYU Dataset
 dfields = [{'depths'},{'images'},{'labels'},{'instances'}]; 
 fpath = 'C:\data\nyu_depth_v2_labeled.mat';
 splitdataset(fpath,dfields);
 
 %% Synthesize Active Brightness image for NYU dataset
  fittingmethod = 'hack';% 'triangles';
  load('BrightnessConst20130222.mat');
  BrightnessConst = imresize(BrightnessConst,size(depths),'lanczos3');
  %  A = gethlut();
  [pixth, A] = genlut([480 640]);
  A = shiftdim(reshape(A,[3 480 640]),1);
  A(:,:,1) = fliplr(A(:,:,1));
  datapath = 'C:\data\nyu_depth_v2_labeled';
  filenames = dir([datapath '\n*']);
 for p = 1:length(filenames)
     fpath = fullfile(datapath,filenames(p).name);
     load(fpath,'depths','images','labels');
     bw = rgb2gray(images);
     albedo = zeros(size(labels));
     for l = 1:max(labels(:))
         albedo(labels==l) = mean(bw(labels==l));
     end
     synthab = synthAB(depths,albedo,BrightnessConst,A,fittingmethod);
     
 end
 
%% Synthesize Active Brightness image for NYU dataset from Scene-SIRFS
nyudir = 'C:\Users\Ryan\Documents\MATLAB\Scene-SIRFS_release1.1\output\NYU';
for x = 1:3
    name = sprintf('scene_%02d.mat', x);
    load(fullfile(nyudir,name));
    [M N] = size(state.height);
    fittingmethod = 'hack';% 'triangles';
    load('BrightnessConst20130222.mat');
    BrightnessConst = imresize(BrightnessConst,[M N],'lanczos3');
    depth = -5*state.height/1000;
    depth(isnan(depth)) = inf;
    albedo = nanmin(state.albedo(:,:,1),[],1);
    [pixth, A] = genlut([M N]);
    A = shiftdim(reshape(A,[3 M N]),1);
    A(:,:,1) = fliplr(A(:,:,1));
    A = A./repmat( sum(A.^2,3),[1 1 3]);
    adotn = abs(dot(A,state.normal,3));
    AB = BrightnessConst .* albedo .* adotn ./ depth.^2;
     synthab(:,:,x) = AB;
  %  synthab(:,:,x) = synthAB(depth,albedo,BrightnessConst,A,fittingmethod);
end
clear depth A pixth albedo fittingmethod name x BrightnessConst AB
%% just testing the different bilateral filter effects
b_uselabels = 0;
w = 4;
GAMMA_C = .001;
GAMMA_S = 4;
[H, W] = size(depths);
enhancedDepth = zeros(size(depths));

% prepare albedo
albedo = zeros(size(labels));

bw = rgb2gray(images);
if ~b_uselabels
    ROI = ones(size(depths));
    enhancedDepth  = blf(depths, w, GAMMA_C, GAMMA_S, [1 1 H W], ROI);
else
    for l = 1:max(labels(:))
        ROI = mean(bw(labels==l));
        albedo(labels==l) = mean(bw(labels==l));
   %     enhancedDepth_layer  = blf(depths, w, GAMMA_C, GAMMA_S, [1 1 H W], ROI);
    %    enhancedDepth(bw(labels==l)) = enhancedDepth_layer(bw(labels==l));
    end
end
albedomedfilt = medfilt2(albedo,[4 4]);albedo(albedo==0) = albedomedfilt(albedo==0);
clear albedomedfilt

 synthimg = synthAB(depths,albedo,BrightnessConst,A,fittingmethod);
     
