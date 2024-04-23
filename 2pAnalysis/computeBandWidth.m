function [Bandwidth] = computeBandWidth( TuningCurve )
% 
% returns the half-maximum, half-bandwidth at the preferred direction of 
% the tuning curve. Assumes a full 360 degree measured tuning curve
%
% Input:
% - TuningCurve: Array with mean response per direction stimulus
%
% Ouput:
% - Bandwidth:      Half bandwidth at half maximum response
% - BandwidthRight: Half bandwidth at the clockwise side of the tuningcurve
% - BandwidthLeft:  Half bandwidth at the counter-clockwise side of the 
%                   tuningcurve
%
%

%% Step 0: Set parameters

HalfMaxPercentage = 1/2;
%    HalfMaxPercentage = 1/sqrt(2);

% set the minimum of the tuning curve to zero
if max(TuningCurve) > 0
    TuningCurve( TuningCurve<0 ) = 0;
else
    warning('TwoPhotonToolbox:TuningCurveBelowZero', ...
        ['Tuning curve does not contain any positive values, ' ...
        'shifted values up by substracting the minimum...']);
    TuningCurve = TuningCurve - min(TuningCurve);
end

%% Step 1: Get the datapoints from the tuning curve
% get number of datapoints on full circle
[~, nDataPoints] = size(TuningCurve);

% calculate width of steps in tuning curve
angleWidth = 360/length(TuningCurve);

%% Step 2: Find the location and value of the peak in the  data
% find array indices with the largest response
PeakIndex = find( TuningCurve == max(TuningCurve) );

% if more than one peak has been found, take the mean of the
% neighboring data points in account
if length(PeakIndex) > 1

    max3peak = zeros(length(PeakIndex),1);
    for d = 1:length(PeakIndex)
        Neighbors = mod( ((PeakIndex(d)-1):1:(PeakIndex(d)+1))-1, length(TuningCurve) ) + 1;
        max3peak(d) = mean( TuningCurve( Neighbors ) );
    end
    MaxMaxIndx = find( max3peak == max(max3peak) );

    % if this doesnt decide, do it with 2 neighbors on each side
    if length(MaxMaxIndx) > 1

        max5peak = zeros(length(PeakIndex),1);
        for d = 1:length(PeakIndex)
            Neighbors = mod( ((PeakIndex(d)-2):1:(PeakIndex(d)+2))-1, length(TuningCurve) ) + 1;
            max5peak(d) = mean( TuningCurve( Neighbors ) );
        end
        MaxMaxIndx = find( max5peak == max(max5peak) );

        % if this doesnt decide, spit out a warning and make a random
        % choice between the maximum responses
        if length(MaxMaxIndx) > 1
            warning('TwoPhotonToolbox:MultiplePreferredDirections', ...
                'Found multiple identical peaks, cannot decide preferred direction, taking random pick...')
            MaxMaxIndx = ceil(rand(1)*length(PeakIndex));
            if MaxMaxIndx == 0
                MaxMaxIndx = 1;
            end
        end                
    end

    PeakIndex = PeakIndex( MaxMaxIndx );
end
peakValue = TuningCurve(PeakIndex);

%% Step 3: Calculate the indexes for the data-points left of the peak and right 
% of the peak
rightSteps = PeakIndex:PeakIndex+(nDataPoints);
leftSteps = PeakIndex:-1:PeakIndex-(nDataPoints);
rightSteps = mod(rightSteps-1,nDataPoints)+1;
leftSteps = mod(leftSteps-1,nDataPoints)+1;

% step through the right half of the tuning curve to find the first 
% datapoint that is 1/sqrt(2) (70.7%) of peak response.
for rs = 1:length(rightSteps)
    pt = TuningCurve( rightSteps(rs) );
    if pt < (peakValue .* HalfMaxPercentage)
        break;
    end
end
rightAngle = (rs-2)*angleWidth;

% interpolate for precise bandwith estimate
intAng = angleWidth - ...
    (( ((peakValue.*HalfMaxPercentage)-TuningCurve(rightSteps(rs))) / ...
    (TuningCurve(rightSteps(rs-1))-TuningCurve(rightSteps(rs))) ) ...
    * angleWidth);

% exact estimate of right angle:
rightAngle = rightAngle + intAng;

if rightAngle < 0
    rightAngle = 0;
end

%% Step 4: Step through the left half of the tuning curve to find the first 
% datapoint that is 1/sqrt(2) (70.7%) of peak response.
for ls = 1:length(leftSteps)
    pt = TuningCurve( leftSteps(ls) );
    if pt < (peakValue .* HalfMaxPercentage)
        break;
    end
end
leftAngle = (ls-2)*angleWidth;

% interpolate for precise bandwith estimate
intAng = angleWidth - ...
    (( ((peakValue.*HalfMaxPercentage)-TuningCurve(leftSteps(ls))) / ...
    (TuningCurve(leftSteps(ls-1))-TuningCurve(leftSteps(ls))) ) ...
    * angleWidth);

% exact estimate of right angle:
leftAngle = leftAngle + intAng;

if leftAngle < 0
    leftAngle = 0;
end

% bandwidth is defined as one-half of the difference between these two angles
Bandwidth = (rightAngle + leftAngle) / 2;
BandwidthRight = rightAngle;
BandwidthLeft = leftAngle;

% if bandwidth is bigger than 180 degrees, give max as answer
if Bandwidth > 180
    Bandwidth = 180;
end
if BandwidthRight > 180
    Bandwidth = 180;
end
if BandwidthLeft > 180
    Bandwidth = 180;
end

