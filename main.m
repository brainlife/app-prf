function [] = main()

if ~isdeployed
	disp('loading paths for IUHPC')
	addpath(genpath('/N/u/brlife/git/jsonlab'))
	addpath(genpath('/N/u/brlife/git/mrTools'))
end

% load my own config.json
config = loadjson('config.json');

% check if stimulus image units were inputted
if isempty(config.stimimageunits)
    display('Stimulus image units were not inputted. App cannot run')
    exit;
else
    stimImageUnits = str2num(config.stimimageunits);
end

% check if frame period was inputted
if isempty(config.frameperiod)
    display('Frame period was not inputted. Default will be used')
    framePeriod = 1.537;
else
    framePeriod = str2num(config.frameperiod);
end

% check if mask manually inputted
if isempty(config.mask)
    display('Mask not selected. All voxels will be analyzed');
    mask = [];
else
    mask = str2num(config.mask);
end

% compute pRF
pRFLife(config.tseries,config.stimimage, framePeriod, stimImageUnits, mask);

end
