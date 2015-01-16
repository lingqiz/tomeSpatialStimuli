function [params] = retPregen(paramName,params)

if ~exist('paramName', 'var')
    error('Must supply a parameter file name');
end
if ~exist('params', 'var'), params = []; end
params = retMenu(params);

% now set rest of the params
params = setRetinotopyParams(params.experiment, params);
[scriptdir] = fileparts(mfilename('fullpath'));
paramsPath = fullfile(scriptdir,'savedParams',paramName);
save(paramsPath,'params');

end