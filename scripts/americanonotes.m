%% Americano data
addpath('C:\data\BrightnessProb\data\code');
addpath('C:\Users\rcrabb\Documents\MATLAB\Robust');

basedir = 'C:\data\BrightnessProb\15perc';
fsflist = dir(fullfile(basedir,'*.fsf'));
sample0 = [];

inlierMM = 40;
minGoodFrames = 3;

for d = 1:length(fsflist)
    fsf_file = fullfile(basedir,fsflist(d).name);

    [status fsf] = ReadFSFFile(fsf_file);
    X = fsf.stream(1).data; X(X==0) = NaN; x = nanmedian(X,3);
    Y = fsf.stream(2).data; Y(Y==0) = NaN; y = nanmedian(Y,3);
    Z = fsf.stream(3).data; Z(Z==0) = NaN; z = nanmedian(Z,3);
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
    Bs = fsf.stream(4).data;
    B = nanmedian(Bs,3);

    % Find the surface norm for this frame
    %   Skip if this has already been done
    if (isempty(sample0))
        h = figure;seq(B);
        sample0rect = round(getrect(h));
        sample0 = zeros(size(B));
        sample0(sample0rect(2):sample0rect(2)+sample0rect(4)-1,sample0rect(1):sample0rect(1)+sample0rect(3)-1) = 1;
        close(h);
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
    B = sum(Bs.*inliermask,3)./max(sum(inliermask,3),1);
    Bs(~inliermask) = NaN;
    Brightness(:,:,d) = B;
    Bstd(:,:,d) = nanstd(Bs,3);
  % refit with the better points
    sample = sort(randsample(length(XYZ),max(round(length(XYZ)/10),min(2200,length(XYZ)))));
    [coefs, P, inliers] = ransacfitplane(XYZ(:,sample), 20, 1);
    coefs = coefs/norm(coefs(1:3))*sign(coefs(3));    
    norms(:,d) = coefs/norm(coefs(1:3));
  %  arrow([0 0 0],norms(:,d));
    CosB = abs(sum(A.*repmat(reshape(norms(1:3,d),[1 1 3]),[M N 1]),3));
    Coss(:,:,d) = CosB;
    Rcalc0 = -norms(4,d)./CosB;
    dRcalc(:,:,d) = Rcalc-Rcalc0;
    Rcalc = Rcalc0;
    
    % Usings assumption that B = I*alpha*cos(beta)/D^2
    %Want to see that each pixel has a contant value of I*alpha over frames
    C(:,:,d) = (B.*Rcalc.^2)./CosB;
    inliermask = abs(R-repmat(Rcalc,[1 1 size(R,3)]))<inlierMM & repmat(sample0,[1 1 size(R,3)]);
    mask = sum(inliermask,3)>10;
    masks(:,:,d) = mask;
    sample0 = [];
    
    
end
[M N F] = size(X);

c = C;
c(~masks) = NaN;
stds = nanstd(c,3);
Cave = sum(C.*masks,3)./max(sum(masks,3),1);