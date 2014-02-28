function [pctsChoi0 pctsChoi1 solchoi0 solchoi1 GTs] = comparefulltermschoi(dirchoi0,dirchoi1)
% Compare full terms Ours vs Choi
filenamesChoi0 = dir([dirchoi0 'w*']);
filenamesChoi1 = dir([dirchoi1 'w*']);
dirlistChoi0 = {filenamesChoi0(:).name};
dirlistChoi1 = {filenamesChoi1(:).name};
solours = [];
solchoi = [];
for testnum = 1:length(dirlistChoi0)

    dataname = [dirchoi0 dirlistChoi0{testnum}];
    load(dataname,'s');
    [pctsChoi0(testnum), iChoi0(testnum)] = max([s(:).pctcorrect]);
    solchoi0(:,:,testnum) = s(iChoi0(testnum)).endstate(:,:,end);
    
    dataname = [dirchoi1 dirlistChoi1{testnum}];
    load(dataname,'s','WrapStateGT');
    [pctsChoi1(testnum), iChoi1(testnum)] = max([s(:).pctcorrect]);
    solchoi1(:,:,testnum) = s(iChoi1(testnum)).endstate(:,:,end);
    GTs(:,:,testnum) = WrapStateGT;
end