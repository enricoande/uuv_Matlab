%% Root directory of this running .m file:
projectRootDir = fileparts(mfilename('fullpath'));

%% Add project directories to path:
addpath(fullfile(projectRootDir,'data'),'-end');
addpath(fullfile(projectRootDir,'dynamics'),'-end');
addpath(fullfile(projectRootDir,'postprocessing'),'-end');
addpath(fullfile(projectRootDir,'preprocessing'),'-end');
addpath(fullfile(projectRootDir,'reinforcement_learning'),'-end');
addpath(fullfile(projectRootDir,'scripts'),'-end');