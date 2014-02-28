%% Load data for moving captures
function [C masks Coss norms D R AB] = BrightnessExpMoving(A, subdir)
% some parameters
inlierMM = 35;
minsamplesize = 2000;

% Set up files names
basedir = 'C:\data\BrightnessCalSVC';
%subdir = 'aset06';
dirlist = dir(fullfile(basedir,strcat(subdir,'*')));
foregroundfile = fullfile(basedir,strcat('mask_',subdir,'.mat'));
workspacefile = fullfile(basedir,strcat('ws_',subdir,'.mat'));
datafile = fullfile(basedir,strcat('data_',subdir,'.mat'));

clear D Bright C masks cs errs Coss Bstd norms
if exist(foregroundfile)
    load(foregroundfile);
    savefg = 0;
else
    fg = [];
    savefg = 1;
    
% if foreground hasn't been created, then do it
    curdir = dirlist(1).name;
    AB = double(loadpgm(fullfile(basedir,curdir,'ConfidenceImage.pgm')));
    figure;seq(AB);
    pause;
    [M N F] = size(AB);
    d = floor(F/10);
    idx = d+1:d:F;
    if (idx(end)~=F)
        idx(end+1) = F;
    end
    h = figure;seq(AB(:,:,1));
    previd = 1;
    fg(:,previd) = max(1,round(getrect(h)));
    close(h);
    for d = idx
        ab = AB(:,:,d);
        h = figure;seq(ab);
        sr = max(1,round(getrect(h)));
        close(h);
        fg(:,d) = sr;
        if (fg(1,d)==fg(1,previd))
            fg(1,previd:d) = fg(1,d);
        else
            fg(1,previd:d) = round(fg(1,previd):(fg(1,d)-fg(1,previd))/(d-previd):fg(1,d));
        end
        if (fg(2,d)==fg(2,previd))
            fg(2,previd:d) = fg(2,d);
        else
            fg(2,previd:d) = round(fg(2,previd):(fg(2,d)-fg(2,previd))/(d-previd):fg(2,d));
        end
        if (fg(3,d)==fg(3,previd))
            fg(3,previd:d) = fg(3,d);
        else
            fg(3,previd:d) = round(fg(3,previd):(fg(3,d)-fg(3,previd))/(d-previd):fg(3,d));
        end
        if (fg(4,d)==fg(4,previd))
            fg(4,previd:d) = fg(4,d);
        else
            fg(4,previd:d) = round(fg(4,previd):(fg(4,d)-fg(4,previd))/(d-previd):fg(4,d));
        end
        previd = d;
    end
    save(foregroundfile,'fg');
    
end

if (exist(datafile))
    load(datafile)
    loaddata = 1;
else
    loaddata = 0;
end

%figure; hold;
%Rs = zeros([M N F]);


if(~loaddata)
    curdir = dirlist(1).name;
    X = loadpgm(fullfile(basedir,curdir,'XImg.pgm'));
    Y = loadpgm(fullfile(basedir,curdir,'YImg.pgm')); 
    Z = loadpgm(fullfile(basedir,curdir,'ZImg.pgm')); 
    AB = double(loadpgm(fullfile(basedir,curdir,'ConfidenceImage.pgm')));
    [M N F] = size(X);
    R = (X.^2+Y.^2+Z.^2).^.5;
else    
    X = R.*repmat(A(:,:,1),[1 1 size(R,3)]);
    Y = R.*repmat(A(:,:,2),[1 1 size(R,3)]);
    Z = R.*repmat(A(:,:,3),[1 1 size(R,3)]);
    [M N F] = size(X);
end

g = waitbar(0,'Fitting hella planes...');
            
for d = 1:F
    
    x = X(:,:,d);
    y = Y(:,:,d);
    z = Z(:,:,d);
    r =  R(:,:,d);
    r2 = R(:,:,d).^2;
    b = AB(:,:,d);
  
    sample0rect = fg(:,d);    
    sample0 = zeros(size(b));
    sample0(sample0rect(2):sample0rect(2)+sample0rect(4)-1,sample0rect(1):sample0rect(1)+sample0rect(3)-1) = 1;
        
    nanmap = ~(~sample0|isnan(x)|isnan(y)|isnan(z)|isnan(A(:,:,1)));
    
    XYZ = [x(nanmap)' ; y(nanmap)' ; z(nanmap)' ];
    sample = sort(randsample(length(XYZ),min(minsamplesize,length(XYZ)-1)));
    [textnull, coefs, P, inliers] = evalc('ransacfitplane(XYZ(:,sample), 30, 1)');
    norms(:,d) = coefs/norm(coefs(1:3))*sign(coefs(3));
    CosB = abs(sum(A.*repmat(reshape(norms(1:3,d),[1 1 3]),[M N 1]),3));
    Rcalc = -norms(4,d)./CosB;
    
    % Go back and find outlying single frames
    inliermask = abs(r-Rcalc)<inlierMM;    
    mask = inliermask & sample0;
    masks(:,:,d) = mask;
    XYZ = [x(mask)' ; y(mask)' ; z(mask)' ];
  % refit with the better points
    samplesize = min(length(XYZ),max(minsamplesize,round(length(XYZ)/100)));
    sample = sort(randsample(length(XYZ),samplesize));
    [textnull, coefs, P, inliers] = evalc('ransacfitplane(XYZ(:,sample), 20, 1)');
    coefs = coefs/norm(coefs(1:3))*sign(coefs(3));    
    norms(:,d) = coefs/norm(coefs(1:3));
  %  arrow([0 0 0],norms(:,d));
    CosB = abs(sum(A.*repmat(reshape(norms(1:3,d),[1 1 3]),[M N 1]),3));
    Coss(:,:,d) = CosB;
    D(:,:,d) = -norms(4,d)./CosB;
    
    % Usings assumption that B = I*alpha*cos(beta)/D^2
    %Want to see that each pixel has a contant value of I*alpha over frames
    C(:,:,d) = (b.*r2)./CosB;
    waitbar(d/F,g);
%     if mod(d,100) == 0
%         save(datafile,'C','masks','Coss','norms','R','AB','subdir','basedir');
%     end
end
close(g);
if (~loaddata)
   save(datafile,'C','masks','Coss','norms','R','AB','subdir','basedir');
end