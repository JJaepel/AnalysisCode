function [activeFrameStack,numberOfActiveEvents,eventOnset,eventDuration] = getActiveFrames(imageStack,expParam)
    % Extracts our the active frames in an imaging stack, and returns them. Useful for determining spontaneous 
    %active events.
    if nargin < 2
        expParam = [];
    end
    % Get image dimensions, and then reshape tif stack into a transposed 2d array ==> (x,y,t) to (x*y,t)
    imgStackSize = size(imageStack);
        x  = imgStackSize(1);
        y  = imgStackSize(2);
        t  = imgStackSize(3);
    imageStack = reshape(imageStack,[x*y t]);

    % Parse out variable input arguments and load default parameters where necessary
    if(~isfield(expParam,'ROI')) % Use the entire image as the default ROI, or use the supplied one
        expParam.ROI = true(x*y,1);
    else
        expParam.ROI = reshape(expParam.ROI,[x*y 1]);
        expParam.ROI = logical(expParam.ROI);
    end
    if(~isfield(expParam,'activityThreshold')),              expParam.activityThreshold        = 0.05;        end % how many standard deviations about the mean response that is needed for a pixel to be considered active
    if(~isfield(expParam,'thresholdForActivePixels')),       expParam.thresholdForActivePixels = 0.25;        end % the proportion of pixels that need to be active in an "active" frame
    if(~isfield(expParam,'averageAllFramesInSingleEvents')), expParam.averageAllFramesInSingleEvents = true; end % whether to collapse consecutive active frames into a single event
    if(~isfield(expParam,'thresholdForActivation')) % threshold for a pixel to be considered active. By default is X*standard deviations above the mean. An absolute or percentile threshold could work.
        useSimpleMethod = true;
        if(useSimpleMethod)
            meanImg   = nanmean(imageStack,2);
            stdDevImg =  nanstd(meanImg);
            expParam.thresholdForActivation = meanImg+expParam.activityThreshold*stdDevImg;
        else
            thresholdStack = imageStack;
            thresholdStack(thresholdStack>0) = NaN;
            expParam.thresholdForActivation = expParam.activityThreshold*nanstd(cat(2,-thresholdStack,thresholdStack),[],2);
            clear thresholdStack
        end
    end

    % Get percent active pixels
    percentActivePixels = zeros(t,1);
    for index = 1:t
        percentActivePixels(index) = sum(imageStack(expParam.ROI(:),index)>expParam.thresholdForActivation);
    end 
    percentActivePixels = percentActivePixels ./sum(expParam.ROI(:));
    activeFrames = find(percentActivePixels>expParam.thresholdForActivePixels);
    assert(~isempty(activeFrames),'No active events found. Try lowering activity threshold');

    % Reshape tif stack back into 3d arrays ==> (x*y,t) to (x,y,t)
    imageStack = reshape(imageStack,[x y t]);

    % Either return all event frames or breakdown active frames into contiguous events and average
    if(expParam.averageAllFramesInSingleEvents)
        differenceInFrameLocation = [10; diff(activeFrames)];
        activeEvents              = find(differenceInFrameLocation>1);
        numberOfActiveEvents      = length(activeEvents);
        activeEvents              = [activeEvents; length(activeFrames)+1];
        activeFrames              = [activeFrames; activeFrames(end)+1]; 
        activeFrameStack          = zeros(x,y,numberOfActiveEvents);
        eventOnset                = zeros(numberOfActiveEvents,1);
        eventDuration             = zeros(numberOfActiveEvents,1);
        for event = 1:numberOfActiveEvents
            selectedFrames = find(activeFrames>=activeFrames(activeEvents(event)) & activeFrames<activeFrames(activeEvents(event+1)));
            activeFrameStack(:,:,event) = mean(imageStack(:,:,activeFrames(selectedFrames)),3);
            eventOnset(event)    = activeFrames(selectedFrames(1));
            eventDuration(event) = length(selectedFrames);
        end
    else
        eventOnset    = activeFrames;
        eventDuration = 1+0*activeFrames; 
        numberOfActiveEvents = length(activeFrames);
        activeFrameStack     = imageStack(:,:,activeFrames);
    end
end
