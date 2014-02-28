function [pctsOurs pctsChoi solours solchoi GTs] = comparefullterms(dirours,dirchoi)
% Compare full terms Ours vs Choi
filenamesOurs = dir([dirours 'w*']);
filenamesChoi = dir([dirchoi 'w*']);
dirlistOurs = {filenamesOurs(:).name};
dirlistChoi = {filenamesChoi(:).name};
solours = [];
solchoi = [];
for testnum = 1:length(dirlistOurs)
    dataname = [dirours dirlistOurs{testnum}];
    load(dataname,'s');
    [pctsOurs(testnum), iOurs(testnum)] = max([s(:).pctcorrect]);
    solours(:,:,testnum) = s(iOurs(testnum)).endstate(:,:,end);
    
    dataname = [dirchoi dirlistChoi{testnum}];
    load(dataname,'s','WrapStateGT');
    [pctsChoi(testnum), iChoi(testnum)] = max([s(:).pctcorrect]);
    solchoi(:,:,testnum) = s(iChoi(testnum)).endstate(:,:,end);
    GTs(:,:,testnum) = WrapStateGT;
end