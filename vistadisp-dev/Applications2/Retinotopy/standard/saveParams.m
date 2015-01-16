function [] = saveParams(params,filename)

[scriptdir] = fileparts(mfilename('fullpath'));
paramsPath = fullfile(scriptdir,'savedParams',filename);
save(paramsPath,'params');

end