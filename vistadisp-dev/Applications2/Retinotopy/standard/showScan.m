function [params] = showScan(paramName)

if ~exist('paramName', 'var')
    error('Must supply a parameter file name');
end

[scriptdir] = fileparts(mfilename('fullpath'));

paramsPath = fullfile(scriptdir,'savedParams',paramName);
F = load(paramsPath);
params = F.params;

