% pRF.m
%
%        $Id:$ 
%      usage: [x y rfWidth r2] = pRFLife(dataFilename, stimImageFilename,framePeriod,stimImageUnits,mask);
%         by: justin gardner
%       date: 05/04/2018
%    purpose: compute pRF analysis given nifti file of data, nifit
%             file of stimulus image and parameters
%             
%             datafilename: nifti file containing functional data from pRF scan (dimensions x,y,z,t)
%             stimImageFilename: stimulus movie (dimensions x, y, t) where x and y are whatever resolution you want for your stim image and t needs to be the same as the datafile
%             stimImageUnits: A 4 vector that give the units for the stimImage as [left bottom deltaX deltaY]. That is, stimImageUnits(1), stimImageUnits(2) are the degrees in visual angle (relative to center of screen) of the left, bottom portion of the stimImage. stimImageUnits(3), stimImageUnits(4) are the difference in degrees between neighboring x and y points of the stimImage (i.e. deltaX, deltaY)
%             mask: vector 3xk of voxels that you want to run pRF for, if empty then will do all voxels
%%
%       e.g.: [polarAngle eccentricity rfWidth r2] = pRFLife('tSeries.nii','stimImage.nii', 1.537, [-16, -11.381 0.405, 0.387],[10 28 16]');
%             with no output arguments, saves files polarAngle.nii, eccentricity.nii, rfWidth.nii and r2.nii
%
function [polarAngle eccentricity rfWidth r2] = pRFLife(dataFilename, stimImageFilename, framePeriod, stimImageUnits, mask)

% load nifti filename 
[d hdr] = mlrImageLoad(dataFilename);

% get scan dims
scanDims = size(d);

% create output images
x = nan(scanDims(1:3));
y = nan(scanDims(1:3));
rfWidth = nan(scanDims(1:3));
r2 = nan(scanDims(1:3));

% open stim image
stimImage.im = mlrImageLoad(stimImageFilename);
stimImage.x = stimImageUnits(1):stimImageUnits(3):((size(stimImage.im,2)-1)*stimImageUnits(3)+stimImageUnits(1));
stimImage.y = stimImageUnits(2):stimImageUnits(4):((size(stimImage.im,1)-1)*stimImageUnits(4)+stimImageUnits(2));
stimImage.t= 0:framePeriod:(scanDims(4)-1)*framePeriod;
[stimImage.x stimImage.y] = meshgrid(stimImage.x,stimImage.y);

% do prefit
prefit.x = min(stimImage.x(:)):10:max(stimImage.x(:));
prefit.y = min(stimImage.y(:)):10:max(stimImage.y(:));
prefit.rfWidth = 1:2:8;
[prefit.x prefit.y prefit.rfWidth] = meshgrid(prefit.x,prefit.y,prefit.rfWidth);
prefit.n = length(prefit.x(:));

% get a canonical HRF
params.lengthInSeconds = 20;
params.timelag = 1;
params.offset = 0;
params.tau = 0.6;
params.exponent = 6;
params.diffOfGamma = 0;
canonicalHRF  = getCanonicalHRF(params,framePeriod);

% init prefitResponse
prefit.response = nan(prefit.n,scanDims(4));

% if the mask is empty, then do all voxels
if isempty(mask)
  % get all coordinates
  [maskX maskY maskZ] = meshgrid(1:scanDims(1),1:scanDims(2),1:scanDims(3));
  % and put them in the mask
  mask(1,:) = maskX(:);
  mask(2,:) = maskY(:);
  mask(3,:) = maskZ(:);
end

% compute prefit models
for iPrefit = 1:prefit.n
  % normalize to 0 mean unit length
  prefit.response(iPrefit,:) = getModelResponse(stimImage,prefit.x(iPrefit),prefit.y(iPrefit),prefit.rfWidth(iPrefit),canonicalHRF,scanDims(4));
end

parpool(str2num(getenv('PROCESS_COUNT')));

% now cycle over voxels in mask and do fit - for now, just do a prefit correlation 
% (i.e. no optimization)
%disppercent(-inf,'(pRFLife) Fitting voxels');
disp('(pRFLife) Fitting voxels');
bestParams_all = nan(4,size(mask,2));

% optimParams
optimParams = optimset('MaxIter',1000,'Display','none');

parfor iVoxel = 1:size(mask,2)
  % get tSeries
  tSeries = squeeze(d(mask(1,iVoxel),mask(2,iVoxel),mask(3,iVoxel),:));
  % normalize tSeries
  tSeries = tSeries-mean(tSeries);
  tSeries = tSeries/sqrt(sum(tSeries.^2));
  % find the best correlation with prefit
  r = prefit.response*tSeries;
  % pick the best model
  [maxR bestModel] = max(r);
  % get initParams
  initParams    = nan(1,3);
  initParams(1) = prefit.x(bestModel);
  initParams(2) = prefit.y(bestModel);
  initParams(3) = prefit.rfWidth(bestModel);

  % optimParams
  optimParams = optimset('MaxIter',1000,'Display','none');

  % do non-linear fit
  [bestParams fval exitflag] = fminsearch(@getModelResidual,initParams,optimParams,tSeries,stimImage,canonicalHRF,scanDims(4));
  % get r for optimal params
  [~,r] = getModelResidual(bestParams,tSeries,stimImage,canonicalHRF,scanDims(4));
  
  bestParams_all(:,iVoxel) = [bestParams(1),bestParams(2),bestParams(3), r^2]';

  %disppercent(iVoxel/size(mask,2));
end
%disppercent(inf);

for iVoxel = 1:size(mask,2)
 % put bestParams into overlays
 x(mask(1,iVoxel),mask(2,iVoxel),mask(3,iVoxel))       = bestParams_all(1,iVoxel);
 y(mask(1,iVoxel),mask(2,iVoxel),mask(3,iVoxel))       = bestParams_all(2,iVoxel);
 rfWidth(mask(1,iVoxel),mask(2,iVoxel),mask(3,iVoxel)) = bestParams_all(3,iVoxel);
 r2(mask(1,iVoxel),mask(2,iVoxel),mask(3,iVoxel))      = bestParams_all(4,iVoxel);   
end

% convert to polarAngle eccentricity
[polarAngle eccentricity] = cart2pol(x,y);

% write out the files
if nargout < 1
  mlrImageSave('polarAngle.nii.gz',polarAngle,hdr);
  mlrImageSave('eccentricity.nii.gz',eccentricity,hdr);
  mlrImageSave('rfWidth.nii.gz',rfWidth,hdr);
  mlrImageSave('r2.nii.gz',r2,hdr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getModelResidual    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [retval r modelResponse] = getModelResidual(params,tSeries,stimImage,canonicalHRF,n)

% get model response
modelResponse = getModelResponse(stimImage,params(1),params(2),params(3),canonicalHRF,n);

% get correlation of model with time series
r = modelResponse * tSeries;

% return 1-r for fminsearch
retval = 1-r;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    modelBoldResponse    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modelBoldResponse = getModelResponse(stimImage,x,y,rfWidth,canonicalHRF,n)

% creat rf
rf = exp(-((stimImage.x-x).^2 + (stimImage.y-y).^2)/(2*rfWidth^2));

% loop as a matrix operation it will generally speed things up significantly
modelNeuralResponse = squeeze(sum(sum(stimImage.im .* repmat(rf,1,1,n))));

%% FIX:For testing with older versions of matlab - use this instead of above
%for i = 1:n
%  modelNeuralResponse(i) = sum(sum(stimImage.im(:,:,i) .* rf));
%end
%modelNeuralResponse = modelNeuralResponse';

% do the convolution
modelBoldResponse = conv(modelNeuralResponse,canonicalHRF.hrf);
% trim length
modelBoldResponse = modelBoldResponse(1:n)';
% normalize
modelBoldResponse = modelBoldResponse-mean(modelBoldResponse);
modelBoldResponse = modelBoldResponse/sqrt(sum(modelBoldResponse.^2));


%%%%%%%%%%%%%%%%%%%%%
%%   getGammaHRF   %%
%%%%%%%%%%%%%%%%%%%%%
function fun = getGammaHRF(time,p)

fun = thisGamma(time,1,p.timelag,p.offset,p.tau,p.exponent)/100;
% add second gamma if this is a difference of gammas fit
if p.diffOfGamma
  fun = fun - thisGamma(time,p.amplitudeRatio,p.timelag2,p.offset2,p.tau2,p.exponent2)/100;
end

%%%%%%%%%%%%%%%%%%%
%%   thisGamma   %%
%%%%%%%%%%%%%%%%%%%
function gammafun = thisGamma(time,amplitude,timelag,offset,tau,exponent)

exponent = round(exponent);
% gamma function
gammafun = (((time-timelag)/tau).^(exponent-1).*exp(-(time-timelag)/tau))./(tau*factorial(exponent-1));

% negative values of time are set to zero,
% so that the function always starts at zero
gammafun(find((time-timelag) < 0)) = 0;

% normalize the amplitude
if (max(gammafun)-min(gammafun))~=0
  gammafun = (gammafun-min(gammafun)) ./ (max(gammafun)-min(gammafun));
end
gammafun = (amplitude*gammafun+offset);


%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getCanonicalHRF   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function hrf = getCanonicalHRF(params,sampleRate)

hrf.time = 0:sampleRate:params.lengthInSeconds;
hrf.hrf = getGammaHRF(hrf.time,params);

% normalize to amplitude of 1
hrf.hrf = hrf.hrf / max(hrf.hrf);

