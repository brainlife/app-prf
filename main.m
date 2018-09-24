function [] = main()

if ~isdeployed
	disp('loading paths for IUHPC')
	addpath(genpath('/N/u/brlife/git/jsonlab'))
	addpath(genpath('/N/u/brlife/git/mrTools'))
	addpath(genpath('/N/u/brlife/git/vistasoft'))
end

% load my own config.json
config = loadjson('config.json');

% compute pRF
pRFLife(config.tseries, config.stimimage, config.frameperiod, config.visual_angle_width, config.visual_angle_height, config.mask, config.prefitOnly);

end
