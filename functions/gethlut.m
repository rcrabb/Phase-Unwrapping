function [A] = gethlut()

%Figure out pixel vectors from calibration tables
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