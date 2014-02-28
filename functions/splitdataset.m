 %% separate NYU dataset
 function [] = splitdataset(fpath,dfields)

 [pathstr,name,EXT] = fileparts(fpath);
 subdir = name;
 subpath = fullfile(pathstr,subdir);
%  if exist(subpath)
%      return
%  end
 %mkdir(subpath);
 
 for p = 1:length(dfields)
     load(fpath,dfields{p});
     eval(['data = ' dfields{p} ';']);
     nd = sum(size(data)>0);
     for i = 1:length(data)
         nameout = fullfile(subpath,[name num2str(i,'%04d'),'.mat']);
         if (nd == 4)         
             eval([dfields{p} '= data(:,:,:,' num2str(i) ');']);
             if ~exist(nameout)
                 save(nameout,dfields{p});
             else
                 save(nameout,dfields{p},'-append');
             end
         elseif (nd == 3)         
             eval([dfields{p} '= data(:,:,' num2str(i) ');']);
             if ~exist(nameout)
                 save(nameout,dfields{p});
             else
                 save(nameout,dfields{p},'-append');
             end
         end
     end
     clear data
 end