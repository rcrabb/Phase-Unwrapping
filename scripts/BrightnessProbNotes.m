% Brightness Probability experiment

%% Figure out pixel vectors from calibration tables
hlut = load('C:\Users\rcrabb\Documents\MATLAB\CaptureScripts\hlut_temp.mat');
N = hlut.cal.header.ncols;
M = hlut.cal.header.nrows;
Xtable = hlut.cal.table.X_Table;
Xtable = reshape(Xtable,[N M])';
Ytable = hlut.cal.table.Y_Table;
Ytable = reshape(Ytable,[N M])';

A = double(Xtable);
A(:,:,2) = Ytable;
A(:,:,3) = 1024;
A = A./repmat((sum(A.^2,3)).^.5,[1 1 3]);
%A = reshape(A,[N*M 3]);
clear Xtable Ytable

%%
addpath('C:\Users\rcrabb\Documents\MATLAB\Robust');
basedir = 'C:\data\BrightnessProb\vlab';
subdir = 'all';
dirlist = dir(fullfile(basedir,strcat(subdir,'*')));
foregroundfile = fullfile(basedir,strcat('mask_',subdir,'.mat'));
workspacefile = fullfile(basedir,strcat('ws_',subdir,'.mat'));
datafile = fullfile(basedir,strcat('data_',subdir,'.mat'));

%%

clear D Bright C masks cs errs Coss Bstd norms
if exist(foregroundfile)
    load(foregroundfile);
    savefg = 0;
else
    fg = [];
    savefg = 1;
end
sample0 = [];

inlierMM = 50;
minGoodFrames = 3;
%figure; hold;

for d = 1:length(dirlist)
    
    % Get median value of frames at each pose
    curdir = dirlist(d).name;
    X = loadpgm(fullfile(basedir,curdir,'XImg.pgm')); x = nanmedian(X,3);
    Y = loadpgm(fullfile(basedir,curdir,'YImg.pgm')); y = nanmedian(Y,3);
    Z = loadpgm(fullfile(basedir,curdir,'ZImg.pgm')); z = nanmedian(Z,3);
	[M N F] = size(X);
    R = (X.^2+Y.^2+Z.^2).^.5;
    r2 = (x.^2+y.^2+z.^2);

%     if (isempty(A))
%         XYZ = reshape(X,[M*N F 1]);
%         XYZ(:,:,2) = reshape(Y,[M*N F 1]);
%         XYZ(:,:,3) = reshape(Z,[M*N F 1]);
%         A = XYZ./repmat((sum(XYZ.^2,3).^.5),[1 1 3]);
%         A = reshape(mean(A,2),[M N 3]);
%     end

    % Keep track of poses
    D(:,:,d) = r2;
    Bs = double(loadpgm(fullfile(basedir,curdir,'ConfidenceImage.pgm')));
    B = nanmedian(Bs,3);
    Bstd(:,:,d) = nanstd(Bs,3);

    % Find the surface norm for this frame
    %   Skip if this has already been done
    if (savefg)
        h = figure;seq(B);
        sr = max(1,round(getrect(h)));
        fg(:,d) = sr;
        close(h);
    else
        sr = fg(:,d);
    end
    sample0 = zeros(size(B));
    for f  =1:floor(size(fg,1)/4)
        sample1 = zeros(size(B));
        sample0rect = sr(f*4-3:f*4);
        sample1(sample0rect(2):sample0rect(2)+sample0rect(4)-1,sample0rect(1):sample0rect(1)+sample0rect(3)-1) = 1;
        sample0 = sample0 | sample1;
    end
        
    nanmap = ~(~sample0|isnan(x)|isnan(y)|isnan(z)|isnan(A(:,:,1)));
    
    XYZ = [x(nanmap)' ; y(nanmap)' ; z(nanmap)' ];
    sample = sort(randsample(length(XYZ),max(round(length(XYZ)/10),min(2200,length(XYZ)))));
    [coefs, P, inliers] = ransacfitplane(XYZ(:,sample), 30, 1);
    norms(:,d) = coefs/norm(coefs(1:3))*sign(coefs(3));
    CosB = abs(sum(A.*repmat(reshape(norms(1:3,d),[1 1 3]),[M N 1]),3));
    Coss(:,:,d) = CosB;
    Rcalc = -norms(4,d)./CosB;
    dR = Rcalc-(r2.^.5);
    
    % Go back and find outlying single frames
    inliermask = abs(R-repmat(Rcalc,[1 1 size(R,3)]))<inlierMM+10 & repmat(sample0,[1 1 size(R,3)]);
    x = sum(X.*inliermask,3)./max(sum(inliermask,3),1);
    y = sum(Y.*inliermask,3)./max(sum(inliermask,3),1);
    z = sum(Z.*inliermask,3)./max(sum(inliermask,3),1);
    mask = sum(inliermask,3)>=minGoodFrames & sample0;
    masks(:,:,d) = mask;
    XYZ = [x(mask)' ; y(mask)' ; z(mask)' ];
  % refit with the better points
    sample = sort(randsample(length(XYZ),round(length(XYZ)/10)));
    [coefs, P, inliers] = ransacfitplane(XYZ(:,sample), 20, 1);
    coefs = coefs/norm(coefs(1:3))*sign(coefs(3));    
    norms(:,d) = coefs/norm(coefs(1:3));
  %  arrow([0 0 0],norms(:,d));
    CosB = abs(sum(A.*repmat(reshape(norms(1:3,d),[1 1 3]),[M N 1]),3));
    Coss(:,:,d) = CosB;
    Rcalc = -norms(4,d)./CosB;
    
    % Usings assumption that B = I*alpha*cos(beta)/D^2
    %Want to see that each pixel has a contant value of I*alpha over frames
    inliermask = abs(R-repmat(Rcalc,[1 1 size(R,3)]))<inlierMM;
    Bmean = sum(Bs.*inliermask,3)./max(sum(inliermask,3),1);
    Bright(:,:,d) = Bmean;
    Dmean = sum(R.*inliermask,3)./max(sum(inliermask,3),1);
    D(:,:,d) = Dmean;
    %C(:,:,d) = (B.*Rcalc.^2)./CosB;
    C(:,:,d) = (Bmean.*Dmean.^2)./CosB;
    mask = sum(inliermask,3)>40&sample0;
    masks(:,:,d) = mask;
    sample0 = mask;
    
    
end
if (savefg)
  %  save(foregroundfile,'fg');
end
c = C;
c(~masks) = NaN;
stds = nanstd(c,3);
Cave = sum(C.*masks,3)./max(sum(masks,3),1);
Cmean = nanmean(c,3);
Cdev = abs(C-repmat(Cave,[1 1 size(C,3)]))./repmat(stds,[1 1 size(C,3)]);
topcdev = sort(Cdev(masks));
Cdev = min(1,Cdev/topcdev(floor(length(topcdev)*.86)));
Cdev = min(1,Cdev/2);

clear B P inliers sample nanmap curdir d basedir subdir dirlist Bs X Y Z x y z sample0 sample1 topcdev 

% save(workspacefile);
%%
% Test the C values as described in document
minSamples = 4;
trainsize = 2/3;
ntrials = 10;
[M N F] = size(masks);
 
% resize variables
% D, Bright, C, masks
GoodPxlMask = reshape(masks, [M*N F]);
GoodPxlIdx = sum(GoodPxlMask,2)>=minSamples;
GoodPxl = logical(GoodPxlMask(GoodPxlIdx,:));

dists = reshape(D, [size(D,1)*size(D,2) size(D,3)]);
dists = dists(GoodPxlIdx,:);

constants = reshape(C, [size(C,1)*size(C,2) size(C,3)]);
constants = constants(GoodPxlIdx,:);

AB = reshape(Bright, [size(Bright,1)*size(Bright,2) size(Bright,3)]);
AB = AB(GoodPxlIdx,:);

an = reshape(Coss, [size(Coss,1)*size(Coss,2) size(Coss,3)]);
an = an(GoodPxlIdx,:);

errs = zeros(size(AB,1),ntrials);
cs = errs;
idx = squeeze(1:size(AB,2));
for i = 1:size(AB,1)
    samples = idx(squeeze(GoodPxl(i,:)));
    nsamples = length(samples);
    ntrain = floor(nsamples*trainsize);
    for t = 1:ntrials
        samples = samples(randperm(length(samples)));
        train = samples(1:ntrain);
        test = samples(ntrain+1:end);
        c = mean(constants(i,train));
        b = AB(i,test);
        d = dists(i,test);
        cosb = an(i,test);
        errs(i,t) = mean( abs((cosb.*ones(size(b))*c./b).^.5-d)./d);
        cs(i,t) = c;
    end
end

% save(workspacefile);

% save(datafile,'D','Bright','C','masks','cs','errs','Coss','Bstd');

%%
% Display results
disp('Distance prediction');
merrs = mean(errs,2);
mean(D(logical(masks)))
h=figure('Color',[1.0 1.0 1.0]);hist(merrs,100)
xlabel('Average Percent Error per Pixel')
ylabel('Pixels')
% figure('Color',[1.0 1.0 1.0]);hist(D(logical(masks)),200);
% xlabel('Distance')
% ylabel('Pixels');
length(merrs)
mean(merrs)
sum(merrs<.1)/length(merrs)
sum(merrs<.05)/length(merrs)
% 
% %show used pixels
% gpmask = reshape(GoodPxlIdx, [M N]);
% dispframe = Bright.*masks.*repmat(gpmask,[1 1 F]);
% dispframe = Bright;
% figure;seq(dispframe);
% 
% %color code error
% Cmed = C;
% Cmed(~logical(masks)) = NaN;
% Cmed = nanmedian(Cmed,3);
%  % put the error into hue values 0-.667
% pcterr = min(1,abs(sqrt(Coss.*repmat(Cmed,[1 1 F])./Bright)-D)./D)*-.667+.667;
% 
%  % see how off the estimate of C is
% cnan = C; cnan(~masks) = NaN;
% stds = nanstd(cnan,3);
% Cave = sum(C.*masks,3)./max(sum(masks,3),1);
% Cdev = abs(C-repmat(Cave,[1 1 size(C,3)]));%./repmat(stds,[1 1 size(C,3)]);
% topcdev = sort(Cdev(masks));
% Cdev = min(1,Cdev/topcdev(floor(length(topcdev)*.99)));
% %Cdev = min(1,Cdev/2);
% pcterr = Cdev*-.667+.667;
%  % color unused frames magenta
% pcterr(~(masks.*repmat(gpmask,[1 1 F]))) = .848;
% %pcterr = -1*pcterr+1;
% dispframecolor = reshape(pcterr,[M N 1 F]);
% dispframecolor(:,:,2,:) = 1;
% dispframecolor(:,:,3,:) = reshape(dispframe/max(dispframe(:)),[M N 1 F])*.75+.25;
% for f = 1:F
%     dfc(:,:,:,f) = uint8(hsv2rgb(dispframecolor(:,:,:,f))*255);
% end
% figure;seq(dfc)


% Find correlation between Distance and sqrt(cosB/Brightness)
[M N F] = size(masks);
minSamples = min(ceil(F*2/3),10);

% resize variables
% D, Bright, C, masks
GoodPxlMask = reshape(masks, [M*N F]);
GoodPxlIdx = sum(GoodPxlMask,2)>=minSamples;
GoodPxl = logical(GoodPxlMask(GoodPxlIdx,:));

dists = reshape(D, [size(D,1)*size(D,2) size(D,3)]);
dists = dists(GoodPxlIdx,:);

constants = reshape(C, [size(C,1)*size(C,2) size(C,3)]);
constants = constants(GoodPxlIdx,:);

AB = reshape(Bright, [size(Bright,1)*size(Bright,2) size(Bright,3)]);
AB = AB(GoodPxlIdx,:);

an = reshape(Coss, [size(Coss,1)*size(Coss,2) size(Coss,3)]);
an = an(GoodPxlIdx,:);

idx = squeeze(1:size(AB,2));
corrs = zeros(size(AB,1),1);
for i = 1:size(AB,1)
    goodCs = constants(i,GoodPxl(i,:));
    adotn = an(i,GoodPxl(i,:));
    b = AB(i,GoodPxl(i,:));
    d = dists(i,GoodPxl(i,:));
    correlation = corrcoef(d,sqrt(adotn./b));
    corrs(i) = correlation(2,1);
end

disp('number of pixels: ');
length(corrs)
mean(corrs)
sum(corrs>.9)/length(corrs)

clear goodCs adotn b d correlation idx an AB constants i GoodPxl GoodPxlIdx GoodPxlMask dists minSamples

%% Perform evaluation without reshaping

% Find least squares estimate of constant
[M N F] = size(masks);
minSamples = min(ceil(F*2/3),10);

% D, Bright, C, masks
GoodPxlMask = masks;%reshape(masks, [M*N F]);
GoodPxlIdx = sum(GoodPxlMask,3)>=minSamples;
%GoodPxl = logical(GoodPxlMask(GoodPxlIdx,:));

idx = squeeze(1:size(masks,3));
corrs = zeros(M,N);
Clse = corrs;
Dmse = corrs;
stdD = corrs;
Derr = zeros(M,N,F);

for m = 1:M
for n = 1:N 
    if (GoodPxlIdx(m,n))
    idxs = idx(squeeze(masks(m,n,:)));
    ns = length(idxs);
    adotn = Coss(m,n,idxs);
    b = Bright(m,n,idxs);
    d = D(m,n,idxs);
    X = sqrt(adotn./b);
    Clse(m,n) = sum(X.*d)/sum(X.*X);
    % Compute mean square error
    Dmse(m,n) = mean((d-Clse(m,n)*X).^2);
    % Compute standard deviation
    stdD(m,n) = sqrt(Dmse(m,n)*ns/(ns-2));
    % Find per frame error
    Derr(m,n,idxs) = d-Clse(m,n)*X;
    correlation = corrcoef(d,X);
    corrs(m,n) = correlation(2,1);
    else
    Clse(m,n) = 0;
    Dmse(m,n) = 0;
    corrs(m,n) = 0;
    end
end
end

clear ns goodCs adotn b d correlation idx an AB constants i GoodPxl GoodPxlIdx GoodPxlMask dists minSamples

%% Display errors in color
errlim = 1000;
cmap = CanestaColormap('StdColorMap')/255;
cmap(2:end-1,:) = cmap(end-1:-1:2,:);

% just so we don't throw away the original
derr = Derr;
zeromap = (derr==0);
% 
% % Separate positive and negative errors
% derr(abs(derr)>errlim) = errlim*sign(derr(abs(derr)>errlim));
% derr(zeromap) = -errlim-10;
% figure;seq('Colormap',cmap,'CLim',[-errlim-1 errlim], derr);
% set(gcf,'Color',[1 1 1]);

% or use absolute
derr = abs(derr);
 maxerr = max(derr(:));
expon = floor(log(maxerr)/log(10));
errlim = ceil((maxerr/10^expon))*10^expon;
figure;seq('Colormap',cmap,'CLim',[0 errlim], derr);
set(gcf,'Color',[1 1 1]);

% maxerr = max(Dmse(:).^.5);
% expon = floor(log(maxerr)/log(10));
% errlim = ceil((maxerr/10^expon))*10^expon;
% figure;seq('Colormap',cmap,'CLim',[0 errlim], Dmse.^.5);
% set(gcf,'Color',[1 1 1]);

figure;hist(corrs(corrs>0).^2,1000)
set(gcf,'Color',[1 1 1]);
xlabel('coefficient of determination')
ylabel('number of pixels')
title('Histogram of Coefficient of Determination');

figure;hist(stdD(stdD>0).^2,200)
set(gcf,'Color',[1 1 1]);
xlabel('Estimate of Variance')
ylabel('number of pixels')
title('Histogram of Estimate of Variance of Error');

clear cmap derr errlim zeromap
%% Plot fitted line and points for some pixel

% D, Bright, C, masks
[M N F] = size(masks);
minSamples = min(ceil(F*1/4));
GoodPxlMask = masks;%reshape(masks, [M*N F]);
GoodPxlIdx = sum(GoodPxlMask,3)>=minSamples;
rep = 1;
while(rep)
    m = ceil(rand*M);
    n = ceil(rand*N);
    rep = ~(GoodPxlIdx(m,n));
end

idx = squeeze(1:size(masks,3));
idxs = idx(squeeze(masks(m,n,:)));
ns = length(idxs);
adotn = squeeze(Coss(m,n,idxs));
b = squeeze(Bright(m,n,idxs));
d = squeeze(D(m,n,idxs));
X = sqrt(adotn./b);
c1 = squeeze(Clse(m,n));
% Y = d;
% figure;plot(X, Y, 'r. ',X,c1*X, 'k-');
% xlabel('\sqrt(cos \beta / B)')
% ylabel('D')
% Check out squared version
Y = d.^2;
X = X.^2;
c2 = sum(X.*Y)/sum(X.*X);
figure;plot(X, Y, 'r. ',X,c2*X, 'k-');
set(gcf,'Color',[1 1 1]);
xlabel('cos(\beta) / B')
ylabel('D^2 (mm^2)')

% Check out unsquared version after finding squared constant
Y = d;
X = sqrt(adotn./b);
figure;plot(X, Y, 'r. ',X,sqrt(c2)*X, 'k-',X,c1*X, 'b-');
set(gcf,'Color',[1 1 1]);
xlabel('$$\sqrt{cos(\beta) / B}$$','Interpreter','latex','FontSize',13)
ylabel('D (mm)')
legend('Measured data','Uses constant from D^2','Uses constant from D','Location','NorthWest')

% Look at errors vs different measurements
err = d - c1*sqrt(adotn./b);
figure;plot(d,err,'r.')
set(gcf,'Color',[1 1 1]);
xlabel('D (mm)')
ylabel('error (mm)')
title('Error vs Distance');

figure;plot(b,err,'r.')
set(gcf,'Color',[1 1 1]);
xlabel('Brightness')
ylabel('error (mm)')
title('Error vs Brightness');

figure;plot(adotn,err,'r.')
set(gcf,'Color',[1 1 1]);
xlabel('cos \beta')
ylabel('error (mm)')
title('Error vs cos \beta');



%% compile data sets

basedir = 'C:\data\BrightnessProb\vlab';
dirlist = dir(fullfile(basedir,strcat('data_','*.mat')));

loadvars = {'D','Bright','C','masks','Coss'};

for v = 1:length(loadvars)
    eval([loadvars{v}, '=[];']);
end
for d = 1:length(dirlist)
    datafile = fullfile(basedir,dirlist(d).name);
for v = 1:length(loadvars)
    lvar = loadvars{v};
    df = load(datafile);
    eval([lvar '= shiftdim([shiftdim(',...
          lvar, ',2);shiftdim(df.',...
          lvar, ',2)],1);']);
    
end
end


% save(fullfile(basedir,strcat('data_all.mat')),'D','Bright','C','masks','Coss');

%% combine all raw data and save as matlab file
basedir = 'C:\data\BrightnessProb\dojo2\static';
subdir = 'all';
dirlist = dir(fullfile(basedir,'*set*'));
foregroundfile = fullfile(basedir,strcat('mask_',subdir,'.mat'));
workspacefile = fullfile(basedir,strcat('ws_',subdir,'.mat'));
datafile = fullfile(basedir,strcat('data_',subdir,'.mat'));

%% Load in all the datas!

clear D Bright C masks cs errs Coss Bstd norms
if exist(foregroundfile)
    load(foregroundfile);
    savefg = 0;
else
    fg = [];
    savefg = 1;
end
sample0 = [];

if (exist(datafile))
    load(datafile)
    loaddata = 1;
else
    loaddata = 0;
end

inlierMM = 50;
minGoodFrames = 30;
%figure; hold;
%Rs = zeros([M N F length(dirlist)]);
xs = zeros([M N length(dirlist)]);
ys = xs;
zs = xs;

for d = 1:length(dirlist)
    
    % Get median value of frames at each pose
    curdir = dirlist(d).name;
  if(~loaddata)
    X = loadpgm(fullfile(basedir,curdir,'XImg.pgm')); x = nanmedian(X,3);
    Y = loadpgm(fullfile(basedir,curdir,'YImg.pgm')); y = nanmedian(Y,3);
    Z = loadpgm(fullfile(basedir,curdir,'ZImg.pgm')); z = nanmedian(Z,3);
    AB = double(loadpgm(fullfile(basedir,curdir,'ConfidenceImage.pgm')));
    b = nanmedian(AB,3);
	[M N F] = size(X);
    R = (X.^2+Y.^2+Z.^2).^.5;
    r2 = (x.^2+y.^2+z.^2);
%     
     Rs(:,:,1:F,d) = single(R);
%     Bs(:,:,1:F,d) = single(AB);
    xs(:,:,d) = x;
    ys(:,:,d) = y;
    zs(:,:,d) = z;
    D(:,:,d) = r2.^.5;
    Bright(:,:,d) = b;
    Bstd(:,:,d) = nanstd(AB,3);
  else
    load(datafile);
    x = xs(:,:,d);
    y = ys(:,:,d);
    z = zs(:,:,d);
    r2 = D(:,:,d).^2;
    b = Bright(:,:,d);
    [M N] = size(x);
    R = double(Rs(:,:,1:F,d));
  end
  
    % Find the surface norm for this frame
    %   Skip if this has already been done
    if (savefg)
        h = figure;seq(b);
        sr = max(1,round(getrect(h)));
        fg(:,d) = sr;
        close(h);
    else
        sr = fg(:,d);
    end
    sample0 = zeros(size(b));
    for f  =1:floor(size(fg,1)/4)
        sample1 = zeros(size(b));
        sample0rect = sr(f*4-3:f*4);
        sample1(sample0rect(2):sample0rect(2)+sample0rect(4)-1,sample0rect(1):sample0rect(1)+sample0rect(3)-1) = 1;
        sample0 = sample0 | sample1;
    end
        
    nanmap = ~(~sample0|isnan(x)|isnan(y)|isnan(z)|isnan(A(:,:,1)));
    
    XYZ = [x(nanmap)' ; y(nanmap)' ; z(nanmap)' ];
    sample = sort(randsample(length(XYZ),max(round(length(XYZ)/10),min(2200,length(XYZ)))));
    [coefs, P, inliers] = ransacfitplane(XYZ(:,sample), 30, 1);
    norms(:,d) = coefs/norm(coefs(1:3))*sign(coefs(3));
    CosB = abs(sum(A.*repmat(reshape(norms(1:3,d),[1 1 3]),[M N 1]),3));
    Rcalc = -norms(4,d)./CosB;
    
    % Go back and find outlying single frames
    inliermask = abs(R-repmat(Rcalc,[1 1 size(R,3)]))<inlierMM+10 & repmat(sample0,[1 1 size(R,3)]);
    x = sum(X.*inliermask,3)./max(sum(inliermask,3),1);
    y = sum(Y.*inliermask,3)./max(sum(inliermask,3),1);
    z = sum(Z.*inliermask,3)./max(sum(inliermask,3),1);
    mask = sum(inliermask,3)>=minGoodFrames & sample0;
    masks(:,:,d) = mask;
    XYZ = [x(mask)' ; y(mask)' ; z(mask)' ];
  % refit with the better points
    sample = sort(randsample(length(XYZ),round(length(XYZ)/10)));
    [coefs, P, inliers] = ransacfitplane(XYZ(:,sample), 20, 1);
    coefs = coefs/norm(coefs(1:3))*sign(coefs(3));    
    norms(:,d) = coefs/norm(coefs(1:3));
  %  arrow([0 0 0],norms(:,d));
    CosB = abs(sum(A.*repmat(reshape(norms(1:3,d),[1 1 3]),[M N 1]),3));
    Coss(:,:,d) = CosB;
    Rcalc = -norms(4,d)./CosB;
    
    % Usings assumption that B = I*alpha*cos(beta)/D^2
    %Want to see that each pixel has a contant value of I*alpha over frames
    C(:,:,d) = (b.*r2)./CosB;
    d
end
if (savefg)
    save(foregroundfile,'fg');
    save(datafile);
end
    
clear curdir X Y Z AB b R r2 h sr fg savefg mask Rcalc CosB sample0 coefs P inliers XYZ x y z sample inliermask
if (~loaddata)
    save(datafile);
end
%% Test the C values

%figure; hold;

r = (repmat(Cave,[1 1 size(Coss,3)]).*Coss./Brightness).^.5;
diff = r-D.^.5;
aad = nansum(abs(diff).*masks./(D.^.5),3)./max(sum(masks,3),1);
aadmask = sum(masks,3)>1&~isnan(aad)&~isinf(aad);
figure;seq(aad.*aadmask+~aadmask);
figure;hist(aad(aadmask),100);

% compare each tilt separately
for i = 1:size(D,3)
    aad = abs(diff(:,:,i)).*masks(:,:,i)./(D(:,:,i).^.5);
aadmask = sum(masks,3)>1&~isnan(aad)&~isinf(aad);
    
end

%% Figure out pixel vectors from calibration tables
hlut = load('C:\Users\rcrabb\Documents\MATLAB\CaptureScripts\hlut_temp.mat');
N = hlut.cal.header.ncols;
M = hlut.cal.header.nrows;
Xtable = hlut.cal.table.X_Table;
Xtable = reshape(Xtable,[N M])';
Ytable = hlut.cal.table.Y_Table;
Ytable = reshape(Ytable,[N M])';

A = double(Xtable);
A(:,:,2) = Ytable;
A(:,:,3) = 1024;
A = A./repmat((sum(A.^2,3)).^.5,[1 1 3]);
%A = reshape(A,[N*M 3]);
clear Xtable Ytable
%% Draw arrows
figure;hold;
for d = 1:length(norms)
    arrow([0 0 0],norms(1:3,d));
end
clear d

%% make new foreground rectangles
% if exist(foregroundfile)
%     load(foregroundfile);
%     savefg = 0;
% else
%     fg = [];
%     savefg = 1;
% end
% sample0 = [];


sfg = floor(size(fg,1)/4);
for d = 1:length(dirlist)
    
    % Get median value of frames at each pose
    curdir = dirlist(d).name;
    
    Bs = double(loadpgm(fullfile(basedir,curdir,'ConfidenceImage.pgm')));
    B = nanmedian(Bs,3);
    for i = 1:sfg
        sample0rect = fg(:,d);
        B(sample0rect(2):sample0rect(2)+sample0rect(4)-1,sample0rect(1):sample0rect(1)+sample0rect(3)-1) = 0;
    end

    h = figure;seq(B);
    sample0rect = max(1,round(getrect(h)));
    fg(sfg*4+1:(1+sfg)*4,d) = sample0rect;
    close(h);
    %save(foregroundfile,'fg');
    
end
%% Find the cosine of angle between surface normal and pixel vector


