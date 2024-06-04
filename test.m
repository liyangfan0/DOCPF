datasets_path='./SatSOT';
results_file='./results';

if ~exist(results_file, 'dir')
    mkdir(results_file);
end
contents = dir(datasets_path);
seq = {};

for i = 1:length(contents)
    if contents(i).isdir && ~strcmp(contents(i).name, '.') && ~strcmp(contents(i).name, '..')
        seq{end+1} = [datasets_path '/' contents(i).name]; 
    end
end


for s=1:length(seq)
   results=runTracker(seq{s},1);
   data = results.res;
   results_file=fullfile(results_file);
   if ~exist(results_file, 'dir')
    mkdir(results_file);
   end
   filePath=fullfile(results_file,[seq{s},'.txt']);
   dlmwrite(filePath, data, 'delimiter', ',');
   disp(['data saved to ', filePath]);
end
