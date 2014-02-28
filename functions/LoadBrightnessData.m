function [Clse Dmse corrs stdD] = LoadBrightnessData(basedir,goodframes)

mingoodframes = 200;
% Load in saved data
%[M N] = size(A);
M = 200; N = 320;
%basedir = 'C:\data\BrightnessProb\dojo2\moving';
Clse = zeros(M,N);
Dmse = Clse;
stdD = Clse;
corrs = Clse;

albedo = .15;
maxb = 1820;

for row = 1:M

rname = strcat('row',num2str(row,'%03d'),'.mat');
load(fullfile(basedir,rname));
mask = logical(mask) & b <= maxb;
save(fullfile(basedir,rname),'d','b','cb','mask');
[N F] = size(d);

for n = 1:N 
    if (sum(mask(n,:).*goodframes)>mingoodframes)
    idxs = logical(squeeze(mask(n,:)) & goodframes);
    ns = sum(idxs);
    adotn = cb(n,idxs);
    bright = b(n,idxs);
    dist = d(n,idxs);
    
%     % Y term is Distance
%     X = sqrt(adotn./bright);
%     Clse(row,n) = sum(X.*dist)/sum(X.*X);
%     % Compute mean square error
%     Dmse(row,n) = mean((dist-Clse(row,n)*X).^2);
%     % Compute standard deviation
%     stdD(row,n) = sqrt(Dmse(row,n)*ns/(ns-2));
%     % Find per frame error
%     %Derr(row,n,idxs) = dist-Clse(row,n)*X;
%     correlation = corrcoef(dist,X);
%     corrs(row,n) = correlation(2,1);
    
    % Y term is Brightness
    X = albedo*adotn./dist.^2;
    Clse(row,n) = sum(X.*bright)/sum(X.*X);
    % Compute mean square error
    Dmse(row,n) = mean(dist-sqrt(Clse(row,n)*albedo*adotn./bright));
    % Compute standard deviation
    stdD(row,n) = sqrt(Dmse(row,n)*ns/(ns-2));
    % Find per frame error
  %  Derr(row,n,idxs) = dist-Clse(row,n)*X;
    correlation = corrcoef(bright,X);
    corrs(row,n) = correlation(2,1);
    
    else
    Clse(row,n) = 0;
    Dmse(row,n) = 0;
    corrs(row,n) = 0;
    end
end

end
save(fullfile(basedir,'collected_stats.mat'),'corrs','Clse','stdD','Dmse');
