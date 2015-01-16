function [params] = runScan()

runName = input('Enter the scan name');

params = retCreateDefaultGUIParams;
params.prescanDuration = 0;
params.period = 240;
params.numCycles = 3;
params.interleaves = [];
params.tr = 2;
params.loadMatrix = [];
params.calibration = [];
params.triggerKey = 't';
params.runName = runName;
params.savePath = '/Users/Shared/Matlab_Toolboxes/vistadisp-dev/Applications2/Retinotopy/standard/storedImagesMatrices';
mkdir(fullfile(params.savePath,runName));
saveMatFile = fullfile(['bars-images-', datestr(clock,'yyyymmdd-HHMM-SS'), '.mat']);
params.saveMatrix = fullfile(params.savePath, params.runName, saveMatFile);

% now set rest of the params
params = setRetinotopyParams(params.experiment, params);

% set response device
params = setRetinotopyDevices(params);

% go
doRetinotopyScan(params);