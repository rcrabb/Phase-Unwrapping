function [R AB] = pgm2mat(directory, force_reload)
% function [RadialDistance ActiveBrightness] = pgm2mat(directory, force_reload)

if nargin < 2
    force_reload = 0;
end
loadfail = 0;

xfile = [directory '\XImg.pgm'];
yfile = [directory '\YImg.pgm'];
zfile = [directory '\ZImg.pgm'];
cfile = [directory '\ConfidenceImage.pgm'];
mfile = [directory '\data.mat'];


if (exist(mfile,'file') && ~force_reload)
    load(mfile);
    if ( exist('AB','var') ~= 1 || exist('R','var') ~= 1)
        disp('failed to load .mat file, loading pgms')
    else
        AB = double(AB);
        return;
    end
end

if ~(exist(xfile,'file')&&exist(yfile,'file')&&exist(zfile,'file')&&exist(cfile,'file'))
    disp('nope')
    return
end


X = loadpgm(xfile);
Y = loadpgm(yfile);
Z = loadpgm(zfile);
AB = uint16(loadpgm(cfile));

R = (X.^2+Y.^2+Z.^2).^.5;
save(mfile,'AB','R');
AB = double(AB);
