function [] = main()

if ~isdeployed
	disp('loading paths for IUHPC')
	addpath(genpath('/N/u/brlife/git/jsonlab'))
	addpath(genpath('/N/u/brlife/git/mrTools'))
	addpath(genpath('/N/u/brlife/git/vistasoft'))
end

% load my own config.json
config = loadjson('config.json');

% check if mask manually inputted
if isempty(config.mask)
    display('Mask not selected. All voxels will be analyzed');
    mask = [];
else
    mask = str2num(config.mask);
end

% compute pRF
pRFLife(config.tseries, config.stimimage, config.frameperiod, config.visual_angle_width, config.visual_angle_height, mask);

end
