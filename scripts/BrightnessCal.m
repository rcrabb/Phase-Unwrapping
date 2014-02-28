% Intensity Calibration

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
basedir = 'C:\data\BrightnessCal';
subdir = 'set95';
dirlist = dir(fullfile(basedir,strcat(subdir,'_*')));
foregroundfile = fullfile(basedir,strcat('mask_',subdir,'.mat'));
workspacefile = fullfile(basedir,strcat('ws_',subdir,'.mat'));
datafile = fullfile(basedir,strcat('data_',subdir,'.mat'));

%%

clear D Bright C masks cs errs Coss Bstd norms
if exist(foregroundfile)
    load(foregroundfile);
    savefg = 0;
else
    %fg = [];
    savefg = 1;
end
sample0 = [];
snap = 7;
[xg, yg] = meshgrid(1:N,1:M);
inlierMM = 40;
minGoodFrames = 1;
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

    % Keep track of poses
    D(:,:,d) = r2;
    Bs = double(loadpgm(fullfile(basedir,curdir,'ConfidenceImage.pgm')));
    B = nanmedian(Bs,3);
    Bstd(:,:,d) = nanstd(Bs,3);

    % Find the surface norm for this frame
    %   Skip if this has already been done
    if (savefg)
        h = figure;seq(B);
        [pts(:,1) pts(:,2)] = getpts(h);
        pts = pts.*(pts>snap);
        pts(:,1) = pts(:,1).*(pts(:,1)<N-snap)+N*(pts(:,1)>=N-snap);
        pts(:,2) = pts(:,2).*(pts(:,2)<M-snap)+M*(pts(:,2)>=M-snap);
%         pts = floor(pts);
        fg(:,:,d) = pts;
        close(h);
    else
        pts = fg(:,:,d);
    end
    sample0 = inpolygon(xg,yg,pts(:,1),pts(:,2));    
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
    Dcalc(:,:,d) = Rcalc;
    C(:,:,d) = (B.*Rcalc.^2)./CosB;
    %C(:,:,d) = (Bmean.*Dmean.^2)./CosB;
    mask = sum(inliermask,3)>6&sample0;
    masks(:,:,d) = mask;
    sample0 = mask;
    
    
end
if (savefg)
    save(foregroundfile,'fg');
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

save(workspacefile);