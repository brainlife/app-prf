% pRFImport.m
%
%        $Id:$ 
%      usage: pRFImport()
%         by: justin gardner
%       date: 05/27/18
%    purpose: Import completed pRF analysis from BrainLife
%
function retval = pRFImport()

% check arguments
if ~any(nargin == [0])
  help pRFImport
  return
end

% get the path where brain life lives
importPath = mlrGetPathStrDialog('.','Select Brainlife download directory (choose any file)');
if isempty(importPath),return,end

% get just directory
importPath = fileparts(importPath);

% load overlays
overlayFilenames = {'r2.nii','rfWidth.nii','x.nii','y.nii'};
for iOverlay = 1:length(overlayFilenames)
  % get filename
  overlayFilename = fullfile(importPath,overlayFilenames{iOverlay});
  % check for file
  if isfile(fullfile(overlayFilename))
    % load the overlay
    overlayData.(stripext(overlayFilenames{iOverlay})) = mlrImageLoad(overlayFilename);
  else
    disp(sprintf('(pRFImport) Could not find overlay %s in %s',overlayFilenames{iOverlay},importPath));
    return
  end
end

% compute the polar angle and eccentricity
[overlayData.polarAngle overlayData.eccentricity] = cart2pol(overlayData.x,overlayData.y);

% make hot colormap
hotColormap = hot(312);
hotColormap = hotColormap(end-255:end,:);

% display the maps
mrDispOverlay({overlayData.r2 overlayData.polarAngle overlayData.eccentricity overlayData.rfWidth},1,'Concatenation',getMLRView,'overlayNames',{'r2','polarAngle','eccentricity','rfWidth'},'clip',{[0 1] [-pi pi] [0 inf] [0 inf]},'range',{[0 1] [-pi pi] [0 15] [0 15]},'colormapType',{'setRangeToMax','normal','normal','normal'},'cmap',{hotColormap hsv(256) copper(256) copper(256)});
keyboard


