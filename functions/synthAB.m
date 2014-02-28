function [AB] = synthAB(depth,albedo,brightnessConst,A,fittingmethod)
% function [AB] = synthAB(depth,albedo,brightnessConst,A,fittingmethod)
%   depth   - depth image measured as radius from lens (instead of Z)
%   albedo  - albedo image (0,1)
%   brightnessConst - array of pixel const values relating depth to brightness
%   A       - pixel direction vector
%   fittingmethod   - how surface normal is computed:
%                       'fitplane', 'triangles', 'cross', 'hack'

% % filter first
% maxdepth = max(depth(:));
% depth = bfilter2(depth/maxdepth,4,[1 .5/maxdepth])*maxdepth;

if (nargin<5)
    fittingmethod = 'cross';
end

[M N] = size(albedo);

if nargin < 4
    A = gethlut();
end
XYZ = repmat(depth,[1 1 3]).*A;

%XYZ(:,:,3) = -XYZ(:,:,3);
% create vectors between adjacent pixels

switch(fittingmethod)
    case 'fitplane'
        w = 5;
        XYZ = padarray(XYZ,[w w],'both');
        albedo = padarray(albedo,[w w],'both');        
        
        % Create waitbar.
        h = waitbar(0,'Fitting normals...');
        set(h,'Name','Surface Normal Progress');
        
        for m = 1+w:M+2*w-w
            for n = 1+w:N+2*w-w
                amask = reshape(albedo(m-w:m+w,n-w:n+w) == albedo(m,n),[(2*w+1)^2 1]);
                pts = reshape(XYZ(m-w:m+w,n-w:n+w,:),[(2*w+1)^2 3]);
                pts = pts(logical(amask),:);
                if size(pts,1) < 3
                    %can't fit a plane
                    mask(m-w,n-w) = 0;
                else
                     pln = FitPlane(pts,size(pts,1));
                     normals(m-w,n-w,1:3) = pln(1:3)*sign(pln(3))/sqrt(sum(pln(1:3).^2));
                     %normals(m-w,n-w,1:3) = pln(1:3)/sqrt(sum(pln(1:3).^2));
                    mask(m-w,n-w) = 1;
                end
            end
            waitbar(m/(M+2*w-w));
        end
        albedo = albedo(1+w:M+2*w-w,1+w:N+2*w-w);
        close(h);
    case 'triangles'

vertshift = padarray(XYZ(2:end,:,:)-XYZ(1:end-1,:,:),[1 0],'both');
horshift = padarray(XYZ(:,2:end,:)-XYZ(:,1:end-1,:),[0 1],'both');

% now recall that for pixel Pij
%   TOPij = vertshift(i,j+1)
%   BOTij = -vertshift(i+1,j+1)
%   LFTij = horshift(i+1,j)
%   RGTij = -horshift(i+1,j+1)

% Upper left normal
normals(:,:,:,1) = cross(vertshift(1:end-1,:,:),horshift(:,1:end-1,:));
% Upper right nomal
normals(:,:,:,2) = cross(vertshift(1:end-1,:,:),-horshift(:,2:end,:));
% Lower right normal
normals(:,:,:,3) = cross(-vertshift(2:end,:,:),-horshift(:,2:end,:));
% Lower left nomal
normals(:,:,:,4) = cross(-vertshift(2:end,:,:),horshift(:,1:end-1,:));

% normalize normals
normals = normals./(repmat(sqrt(sum(normals.^2,3)).*sign(normals(:,:,3,:)),[1 1 3 1]));

% we should only use the normal if pixels are on the same surface (indicated by albedo)
mask = zeros([M N 4]);
mask(2:end,2:end,1) = albedo(2:end,2:end)==albedo(1:end-1,2:end) & albedo(2:end,2:end)==albedo(2:end,1:end-1);
mask(2:end,1:end-1,2) = albedo(2:end,1:end-1)==albedo(1:end-1,1:end-1) & albedo(2:end,1:end-1)==albedo(2:end,2:end);
mask(1:end-1,1:end-1,3) = albedo(1:end-1,1:end-1,1)==albedo(2:end,1:end-1) & albedo(1:end-1,1:end-1,1)==albedo(1:end-1,2:end);
mask(1:end-1,2:end,4) = albedo(1:end-1,2:end,1)==albedo(2:end,2:end) & albedo(1:end-1,2:end)==albedo(1:end-1,1:end-1);

% average normals
normals = sum(normals.*repmat(reshape(mask,[M N 1 4]),[1 1 3 1]),4)./repmat(sum(mask,3),[1 1 3]);

    case 'cross'

    vertshift = padarray(XYZ(3:end,:,:)-XYZ(1:end-2,:,:),[1 0],'both');
    horshift = padarray(-XYZ(:,3:end,:)+XYZ(:,1:end-2,:),[0 1],'both');
    diagshiftul = padarray(XYZ(3:end,1:end-2,:)-XYZ(1:end-2,3:end,:),[1 1],'both');
    diagshiftdl = padarray(XYZ(3:end,3:end,:)-XYZ(1:end-2,1:end-2,:),[1 1],'both');
%     diagshiftul = padarray(XYZ(3:end,3:end,:)-XYZ(1:end-2,1:end-2,:),[1 1],'both');
%     diagshiftdl = padarray(XYZ(1:end-2,3:end,:)-XYZ(3:end,1:end-2,:),[1 1],'both');
    
    % for debugging purposes, normalize these vectors
    vertshift = vertshift./repmat(sqrt(sum(vertshift.^2,3)),[1 1 3]);vertshift(isnan(vertshift))=0;
    horshift = horshift./repmat(sqrt(sum(horshift.^2,3)),[1 1 3]);horshift(isnan(horshift))=0;
    diagshiftul = diagshiftul./repmat(sqrt(sum(diagshiftul.^2,3)),[1 1 3]);diagshiftul(isnan(diagshiftul))=0;
    diagshiftdl = diagshiftdl./repmat(sqrt(sum(diagshiftdl.^2,3)),[1 1 3]);diagshiftdl(isnan(diagshiftdl))=0;

    % up-down normal
    norms(:,:,:,1) = cross(vertshift,horshift);
    % diagonal normal
    norms(:,:,:,2) = cross(diagshiftdl,diagshiftul);

    % normalize normals
    normals = norms./(repmat(sqrt(sum(norms.^2,3)).*sign(norms(:,:,3,:)),[1 1 3 1]));
    %normals = normals.*sign(norms(:,:,3,:));
    mask = zeros([M N 2]);
    mask(2:end-1,2:end-1,1) = albedo(1:end-2,2:end-1)==albedo(3:end,2:end-1) & albedo(2:end-1,1:end-2)==albedo(2:end-1,3:end);
    mask(2:end-1,2:end-1,2) = albedo(1:end-2,1:end-2)==albedo(3:end,3:end) & albedo(3:end,1:end-2)==albedo(1:end-2,3:end);

    % average normals
    normals1 = sum(normals.*repmat(reshape(mask,[M N 1 2]),[1 1 3 1]),4)./repmat(sum(mask,3),[1 1 3]);
    % hack to try to fix bad averaging
    normals(:,:,:,2) = -normals(:,:,:,2);    
    normals1(:,:,:,2) = sum(normals.*repmat(reshape(mask,[M N 1 2]),[1 1 3 1]),4)./repmat(sum(mask,3),[1 1 3]);
    normsize = squeeze(sqrt(sum(normals1.^2,3)));
    [~,normidx] = max(normsize,[],3);
    normzero = normidx==1;normzero(:,:,2) = normidx==2;
    normzero = repmat(reshape(normzero,[M N 1 2]),[1 1 3 1]);
    normals = sum(normals1.*normzero,4);
    
    case 'hack'

    vertshift = padarray(XYZ(3:end,:,:)-XYZ(1:end-2,:,:),[1 0],'both');
    horshift = padarray(-XYZ(:,3:end,:)+XYZ(:,1:end-2,:),[0 1],'both');
    diagshiftul = padarray(XYZ(3:end,1:end-2,:)-XYZ(1:end-2,3:end,:),[1 1],'both');
    diagshiftdl = padarray(XYZ(3:end,3:end,:)-XYZ(1:end-2,1:end-2,:),[1 1],'both');
%     diagshiftul = padarray(XYZ(3:end,3:end,:)-XYZ(1:end-2,1:end-2,:),[1 1],'both');
%     diagshiftdl = padarray(XYZ(1:end-2,3:end,:)-XYZ(3:end,1:end-2,:),[1 1],'both');
    
    % for debugging purposes, normalize these vectors
    vertshift = vertshift./repmat(sqrt(sum(vertshift.^2,3)),[1 1 3]);vertshift(isnan(vertshift))=0;
    horshift = horshift./repmat(sqrt(sum(horshift.^2,3)),[1 1 3]);horshift(isnan(horshift))=0;
    diagshiftul = diagshiftul./repmat(sqrt(sum(diagshiftul.^2,3)),[1 1 3]);diagshiftul(isnan(diagshiftul))=0;
    diagshiftdl = diagshiftdl./repmat(sqrt(sum(diagshiftdl.^2,3)),[1 1 3]);diagshiftdl(isnan(diagshiftdl))=0;

    % up-down normal
    norms(:,:,:,1) = cross(vertshift,horshift);
    % diagonal normal
    norms(:,:,:,2) = cross(diagshiftdl,diagshiftul);

    % normalize normals
    normals = norms./(repmat(sqrt(sum(norms.^2,3)).*sign(norms(:,:,3,:)),[1 1 3 1]));
    %normals = normals.*sign(norms(:,:,3,:));
    mask = zeros([M N 2]);
    mask(2:end-1,2:end-1,1) = albedo(1:end-2,2:end-1)==albedo(3:end,2:end-1) & albedo(2:end-1,1:end-2)==albedo(2:end-1,3:end);
    mask(2:end-1,2:end-1,2) = albedo(1:end-2,1:end-2)==albedo(3:end,3:end) & albedo(3:end,1:end-2)==albedo(1:end-2,3:end);

    % average normals
    normals1 = sum(normals.*repmat(reshape(mask,[M N 1 2]),[1 1 3 1]),4)./repmat(sum(mask,3),[1 1 3]);
    % hack to try to fix bad averaging
    normals(:,:,:,2) = -normals(:,:,:,2);    
    normals1(:,:,:,2) = sum(normals.*repmat(reshape(mask,[M N 1 2]),[1 1 3 1]),4)./repmat(sum(mask,3),[1 1 3]);
    normsize = squeeze(sqrt(sum(normals1.^2,3)));
    [~,normidx] = max(normsize,[],3);
    normzero = normidx==1;normzero(:,:,2) = normidx==2;
    normzero = repmat(reshape(normzero,[M N 1 2]),[1 1 3 1]);
    normals = sum(normals1.*normzero,4);
    
    n1 = normals;
    m1 = sum(mask,3)>0;
    
    % zap out NaN!
    m1(isnan(n1(:,:,1))) = 0;
    n1(isnan(n1)) = 0;
    
    % part 2, use triangles one at a time
    
    vertshift = padarray(XYZ(2:end,:,:)-XYZ(1:end-1,:,:),[1 0],'both');
    horshift = padarray(XYZ(:,2:end,:)-XYZ(:,1:end-1,:),[0 1],'both');

    normals(:,:,:,1) = cross(vertshift(1:end-1,:,:),horshift(:,1:end-1,:));
    normals(:,:,:,2) = cross(vertshift(1:end-1,:,:),-horshift(:,2:end,:));
    normals(:,:,:,3) = cross(-vertshift(2:end,:,:),-horshift(:,2:end,:));
    normals(:,:,:,4) = cross(-vertshift(2:end,:,:),horshift(:,1:end-1,:));

    % normalize normals
    normals = normals./(repmat(sqrt(sum(normals.^2,3)).*sign(normals(:,:,3,:)),[1 1 3 1]));

    % we should only use the normal if pixels are on the same surface (indicated by albedo)
    mask = zeros([M N 4]);
    mask(2:end,2:end,1) = albedo(2:end,2:end)==albedo(1:end-1,2:end) & albedo(2:end,2:end)==albedo(2:end,1:end-1);
    mask(2:end,1:end-1,2) = albedo(2:end,1:end-1)==albedo(1:end-1,1:end-1) & albedo(2:end,1:end-1)==albedo(2:end,2:end);
    mask(1:end-1,1:end-1,3) = albedo(1:end-1,1:end-1,1)==albedo(2:end,1:end-1) & albedo(1:end-1,1:end-1,1)==albedo(1:end-1,2:end);
    mask(1:end-1,2:end,4) = albedo(1:end-1,2:end,1)==albedo(2:end,2:end) & albedo(1:end-1,2:end)==albedo(1:end-1,1:end-1);

    for tri = 1:4
        m2 = repmat( ~m1 & mask(:,:,tri) , [1 1 3]);
        n1 = n1.*(~m2) + normals(:,:,:,tri).*m2;
        m1 = m1 | m2(:,:,1);
    end
    mask = m1;
    normals = n1;

end

normals = normals./(repmat(sqrt(sum(normals.^2,3)),[1 1 3]));                
adotn = abs(dot(A,normals,3));
adotn(sum(mask,3)==0)=NaN;
w = 2;
for m = 1:M
    for n = 1:N
        if sum(mask(m,n,:))==0
            adotn(m,n) = nanmedian(vec(adotn(max(1,m-w):min(M,m+w),max(1,n-w):min(N,n+w))));
        end
    end
end
adotn(isnan(adotn)) = 0;

% Brightness = BrightnessConst * albedo * cos(Beta) / Distance^2

AB = brightnessConst .* albedo .* adotn ./ depth.^2;