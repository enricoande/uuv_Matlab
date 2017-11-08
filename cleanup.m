%% Root directory of this running .m file:
projectRootDir = fileparts(mfilename('fullpath'));

%% Remove project directories from path:
rmpath(fullfile(projectRootDir,'data'));
rmpath(fullfile(projectRootDir,'dynamics'));
rmpath(fullfile(projectRootDir,'postprocessing'));
rmpath(fullfile(projectRootDir,'preprocessing'));
rmpath(fullfile(projectRootDir,'reinforcement_learning'));
rmpath(fullfile(projectRootDir,'scripts'));

%% leave no trace...
clear projectRootDir;
clear;