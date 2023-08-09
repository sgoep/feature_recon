function E = mycanny(Ix,Iy,thresh)
% Matlab-Implementation of Canny
% Ix,Iy is a smooth approximation of the partial derivatives
% THRESH is a two element vector, the first
%   element of which specifies the lower threshold for edge strength, below
%   which all edges are disregarded. The second element specifies the
%   higher threshold, above which all edge pixels are preserved. The range
%   of values allowed is between 0 and 1. If a scalar is specified, it is
%   used as the high threshold, and the low threshold is computed as
%   0.4*high threshold. If you do not specify THRESH, or if THRESH is empty
%   ([]), EDGE chooses low and high values automatically.


% Magic numbers
PercentOfPixelsNotEdges = .7; % Used for selecting thresholds
ThresholdRatio = .4;          % Low thresh is this fraction of the high.

% Calculate Magnitude of Gradient
magGrad = hypot(Ix, Iy);

% Normalize for threshold selection
magmax = max(magGrad(:));
if magmax > 0
    magGrad = magGrad / magmax;
end

% Determine Hysteresis Thresholds
[lowThresh, highThresh] = selectThresholds(thresh, magGrad, PercentOfPixelsNotEdges, ThresholdRatio,[]);

% Perform Non-Maximum Suppression Thining and Hysteresis Thresholding of Edge
% Strength
E = thinAndThreshold(Ix, Iy, magGrad, lowThresh, highThresh);

% thresh = [lowThresh highThresh];

function [lowThresh, highThresh] = selectThresholds(thresh, magGrad, PercentOfPixelsNotEdges, ThresholdRatio, ~)

[m,n] = size(magGrad);

% Select the thresholds
if isempty(thresh)
    counts=imhist(magGrad, 64);
    highThresh = find(cumsum(counts) > PercentOfPixelsNotEdges*m*n,...
        1,'first') / 64;
    lowThresh = ThresholdRatio*highThresh;
elseif length(thresh)==1
    highThresh = thresh;
    if thresh>=1
        error(message('images:edge:singleThresholdOutOfRange'))
    end
    lowThresh = ThresholdRatio*thresh;
elseif length(thresh)==2
    lowThresh = thresh(1);
    highThresh = thresh(2);
    if (lowThresh >= highThresh) || (highThresh >= 1)
        error(message('images:edge:thresholdOutOfRange'))
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : thinAndThreshold
%
function H = thinAndThreshold(dx, dy, magGrad, lowThresh, highThresh)
% Perform Non-Maximum Suppression Thining and Hysteresis Thresholding of
% Edge Strength

% We will accrue indices which specify ON pixels in strong edgemap
% The array e will become the weak edge map.

% E = builtin('_cannyFindLocalMaximaMex', dx,dy,magGrad,lowThresh);
E = images.internal.builtins.cannyFindLocalMaxima(dx,dy,magGrad,lowThresh);

if ~isempty(E)
    [rstrong,cstrong] = find(magGrad>highThresh & E);

    if ~isempty(rstrong) % result is all zeros if idxStrong is empty
        H = bwselect(E, cstrong, rstrong, 8);
    else
        H = false(size(E));
    end
else
    H = false(size(E));
end