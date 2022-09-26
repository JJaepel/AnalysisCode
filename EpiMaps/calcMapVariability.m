function analysis = calcMapVariability(analysisParams, analysis)

drawingNr = 100;

%% 1.) Overlap of filled Maps between individual trials 
for stim = 1:size(analysis.(analysisParams.field).zScore,3)
   for firstTrial=1:size(analysis.(analysisParams.field).zScore,4)
       firstMap = analysis.(analysisParams.field).indStimTrialMaps(:,:,stim, firstTrial);
       
       %how many of the pixel are active in each area?
       actPx(stim,firstTrial) = nnz(firstMap); %for the whole window
       actPxA19(stim,firstTrial) = nnz(firstMap & analysis.maskA19); %for A19
       actPxV1(stim,firstTrial) = nnz(firstMap & analysis.maskV1); %for V1
       
       for secondTrial = 1:size(analysis.(analysisParams.field).zScore,4) %compare each individual map with each other
           secondMap = analysis.(analysisParams.field).indStimTrialMaps(:,:,stim, secondTrial);
           
           %what is the overlay between the two trials?
           overlapMask = firstMap & secondMap; %for the whole window
           overlapA19 = overlapMask & analysis.maskA19; % for A19
           overlapV1 = overlapMask & analysis.maskV1; %for A19
           
           %calculate the number of pixel for each of them;
           pxlOverlapTrial(firstTrial, secondTrial) = nnz(overlapMask);
           pxlOverlapA19Trial(firstTrial, secondTrial) = nnz(overlapA19);
           pxlOverlapV1Trial(firstTrial, secondTrial) = nnz(overlapV1);
           
           %calculate % of overlapping pixel 
           perOverlapTrial(firstTrial, secondTrial) = pxlOverlapTrial(firstTrial, secondTrial)/actPx(stim,firstTrial);
           perOverlapA19Trial(firstTrial, secondTrial) = pxlOverlapA19Trial(firstTrial, secondTrial)/actPxA19(stim,firstTrial);
           perOverlapV1Trial(firstTrial, secondTrial) = pxlOverlapV1Trial(firstTrial, secondTrial)/actPxV1(stim,firstTrial);
       end
   end
   
   %remove diagonal & save in array
   pxlOverlapTrial(find(eye(size(pxlOverlapTrial)))) = NaN;
   pxlOverlapA19Trial(find(eye(size(pxlOverlapA19Trial)))) = NaN;
   pxlOverlapV1Trial(find(eye(size(pxlOverlapV1Trial)))) = NaN;
   perOverlapTrial(find(eye(size(perOverlapTrial)))) = NaN;
   perOverlapA19Trial(find(eye(size(perOverlapA19Trial)))) = NaN;
   perOverlapV1Trial(find(eye(size(perOverlapV1Trial)))) = NaN;
      
   analysis.(analysisParams.field).pxlOverlay(stim, :) = reshape(pxlOverlapTrial, size(pxlOverlapTrial,1)^2,1);
   analysis.(analysisParams.field).pxlOverlapA19(stim, :) = reshape(pxlOverlapA19Trial, size(pxlOverlapA19Trial,1)^2,1);
   analysis.(analysisParams.field).pxlOverlapV1(stim, :) = reshape(pxlOverlapV1Trial, size(pxlOverlapV1Trial,1)^2,1);
   analysis.(analysisParams.field).perOverlap(stim, :) = reshape(perOverlapTrial, size(perOverlapTrial,1)^2,1);
   analysis.(analysisParams.field).perOverlapA19(stim, :) = reshape(perOverlapA19Trial, size(perOverlapA19Trial,1)^2,1);
   analysis.(analysisParams.field).perOverlapV1(stim, :) = reshape(perOverlapV1Trial, size(perOverlapV1Trial,1)^2,1);
end

%% 2.) Compare to overlap with a randomly drawn map 
for stim = 1:size(analysis.(analysisParams.field).zScore,3)
   for trial=1:size(analysis.(analysisParams.field).zScore,4)
       firstMap = analysis.(analysisParams.field).indStimTrialMaps(:,:,stim, trial);
       for drawing = 1:drawingNr
            randomStim = randi([1 size(analysis.(analysisParams.field).zScore,3)-1]);
            randomTrial= randi([1 size(analysis.(analysisParams.field).zScore,4)]);
            secondMap = analysis.(analysisParams.field).indStimTrialMaps(:,:,randomStim, randomTrial);

            overlapMaskRnd = firstMap & secondMap; %for the whole window
            overlapA19Rnd = overlapMaskRnd & analysis.maskA19; % for A19
            overlapV1Rnd = overlapMaskRnd & analysis.maskV1; %for A19

            %calculate the number of pixel for each of them;
            pxlOverlapTrialRnd(trial, drawing) = nnz(overlapMaskRnd);
            pxlOverlapA19TrialRnd(trial, drawing) = nnz(overlapA19Rnd);
            pxlOverlapV1TrialRnd(trial, drawing) = nnz(overlapV1Rnd);

            %calculate % of overlapping pixel 
            perOverlapTrialRnd(trial, drawing) = pxlOverlapTrialRnd(trial, drawing)/actPx(stim,trial);
            perOverlapA19TrialRnd(trial, drawing) = pxlOverlapA19TrialRnd(trial, drawing)/actPxA19(stim,trial);
            perOverlapV1TrialRnd(trial, drawing) = pxlOverlapV1TrialRnd(trial, drawing)/actPxV1(stim,trial);

        end
   end
   analysis.(analysisParams.field).pxlOverlayRnd(stim, :) = reshape(pxlOverlapTrialRnd, drawingNr *size(analysis.(analysisParams.field).zScore,4),1);
   analysis.(analysisParams.field).pxlOverlapA19Rnd(stim, :) = reshape(pxlOverlapA19TrialRnd, drawingNr *size(analysis.(analysisParams.field).zScore,4),1);
   analysis.(analysisParams.field).pxlOverlapV1Rnd(stim, :) = reshape(pxlOverlapV1TrialRnd, drawingNr *size(analysis.(analysisParams.field).zScore,4),1);
   analysis.(analysisParams.field).perOverlapRnd(stim, :) = reshape(perOverlapTrialRnd, drawingNr *size(analysis.(analysisParams.field).zScore,4),1);
   analysis.(analysisParams.field).perOverlapA19Rnd(stim, :) = reshape(perOverlapA19TrialRnd, drawingNr *size(analysis.(analysisParams.field).zScore,4),1);
   analysis.(analysisParams.field).perOverlapV1Rnd(stim, :) = reshape(perOverlapV1TrialRnd, drawingNr *size(analysis.(analysisParams.field).zScore,4),1);
end

%% 3.) Correlation of activity maps
%restrict area for corr analysis
changeTime = squeeze(std(squeeze(analysis.rawF.roi.stimResponseTrace(:,1,1,:,:)),[],1)); %how much is each pixel changing over time
threshold = (changeTime-min(changeTime(:)))/(max(changeTime(:))-min(changeTime(:)));
threshold(threshold<0.1)=0; threshold(threshold>0.1)=1; %threshold all pixel -> removes inactive pixels
objects = bwlabel(threshold); %find biggest objects
objects(objects == 0) = NaN;
[biggestObject, ~]= mode(objects(:));
threshold(objects ~=biggestObject) = 0;
xVals = find(sum(threshold,2)>5); yVals = find(sum(threshold,1)>5); %find the boundaries of that object

for stim = 1:size(analysis.(analysisParams.field).zScore,3)
   for firstTrial=1:size(analysis.(analysisParams.field).zScore,4)
       firstMap = analysis.(analysisParams.field).indStimTrialMaps(:,:,stim, firstTrial);
       for secondTrial = 1:size(analysis.(analysisParams.field).zScore,4) %compare each individual map with each other
           secondMap = analysis.(analysisParams.field).indStimTrialMaps(:,:,stim, secondTrial);
           
           %how much do both trials correlate? 
           corrMapsTrial(firstTrial, secondTrial) = corr2(firstMap(xVals, yVals), secondMap(xVals, yVals));
           corrMapsA19Trial(firstTrial, secondTrial) = corr2(firstMap(xVals, yVals) & analysis.maskA19(xVals, yVals), secondMap(xVals, yVals) & analysis.maskA19(xVals, yVals));
           corrMapsV1Trial(firstTrial, secondTrial) = corr2(firstMap(xVals, yVals) & analysis.maskV1(xVals, yVals), secondMap(xVals, yVals) & analysis.maskV1(xVals, yVals));
       end
   end
    
    %remove diagonal & save in array
    corrMapsTrial(find(eye(size(corrMapsTrial)))) = NaN;
    corrMapsA19Trial(find(eye(size(corrMapsA19Trial)))) = NaN;
    corrMapsV1Trial(find(eye(size(corrMapsV1Trial)))) = NaN;
    
    analysis.(analysisParams.field).corrMaps(stim, :) = reshape(corrMapsTrial, size(corrMapsTrial,1)^2,1);
    analysis.(analysisParams.field).corrMapsA19(stim, :) = reshape(corrMapsA19Trial, size(corrMapsA19Trial,1)^2,1);
    analysis.(analysisParams.field).corrMapsV1(stim, :) = reshape(corrMapsV1Trial, size(corrMapsV1Trial,1)^2,1);
end

%% 4.) Compare to correlation with a randomly drawn map
for stim = 1:size(analysis.(analysisParams.field).zScore,3)
    for trial = 1:size(analysis.(analysisParams.field).zScore,4)
        firstMap = analysis.(analysisParams.field).indStimTrialMaps(:,:,stim, trial);
        for drawing = 1:drawingNr
            randomStim = randi([1 size(analysis.(analysisParams.field).zScore,3)-1]);
            randomTrial= randi([1 size(analysis.(analysisParams.field).zScore,4)]);
            secondMap = analysis.(analysisParams.field).indStimTrialMaps(:,:,randomStim, randomTrial);
            
            %how much do both trials correlate? 
            corrMapsTrialRnd(drawing, trial) = corr2(firstMap(xVals, yVals), secondMap(xVals, yVals));
            corrMapsA19TrialRnd(drawing, trial) = corr2(firstMap(xVals, yVals) & analysis.maskA19(xVals, yVals), secondMap(xVals, yVals) & analysis.maskA19(xVals, yVals));
            corrMapsV1TrialRnd(drawing, trial) = corr2(firstMap(xVals, yVals) & analysis.maskV1(xVals, yVals), secondMap(xVals, yVals) & analysis.maskV1(xVals, yVals));
        end
    end
    analysis.(analysisParams.field).corrMapsRnd(stim, :) = reshape(corrMapsTrialRnd, drawingNr *size(analysis.(analysisParams.field).zScore,4),1);
    analysis.(analysisParams.field).corrMapsA19Rnd(stim, :) = reshape(corrMapsA19TrialRnd, drawingNr *size(analysis.(analysisParams.field).zScore,4),1);
    analysis.(analysisParams.field).corrMapsV1Rnd(stim, :) = reshape(corrMapsV1TrialRnd, drawingNr *size(analysis.(analysisParams.field).zScore,4),1);
end

%% 5.) If you add all individual trials on top of each other, how much of the possibly used area is used/how much are the pixel converging?
overlapAll=zeros(size(analysis.(analysisParams.field).zScore,1), size(analysis.(analysisParams.field).zScore,2),size(analysis.(analysisParams.field).zScore,3));
boundaryAreaA19=zeros(size(analysis.(analysisParams.field).zScore,1), size(analysis.(analysisParams.field).zScore,2),size(analysis.(analysisParams.field).zScore,3));
boundaryAreaV1=zeros(size(analysis.(analysisParams.field).zScore,1), size(analysis.(analysisParams.field).zScore,2),size(analysis.(analysisParams.field).zScore,3));

for stim = 1:size(analysis.(analysisParams.field).zScore,3)
   for trial=1:size(analysis.(analysisParams.field).zScore,4)
       Map = analysis.(analysisParams.field).indStimTrialMaps(:,:,stim, trial);
       overlapAll(:,:,stim) = overlapAll(:,:,stim)+Map;
   end
   boundaryAreaA19(:,:,stim) = overlapAll(:,:,stim).*analysis.maskA19;
   boundaryAreaV1(:,:,stim) = overlapAll(:,:,stim).*analysis.maskV1;
end
%to count pixel, make sure that all of them are 0 (not used at all) or 1
%(used at least once)
boundaryAreaV1(boundaryAreaV1 > 0) = 1;
boundaryAreaA19(boundaryAreaA19 > 0) = 1;

SumBoundAreaA19 = squeeze(nansum(boundaryAreaA19,[1 2]));
SumBoundAreaV1 = squeeze(nansum(boundaryAreaV1,[1 2]));

%simulate random distribution of active pixel
A19Pixels = find(analysis.maskA19 == 1); %what are the indicators of the A19/V1 masks
V1Pixels = find(analysis.maskV1 == 1);
for drawing = 1:drawingNr %go up to 1000
   overlapAllRnd=zeros(size(analysis.(analysisParams.field).zScore,1), size(analysis.(analysisParams.field).zScore,2),size(analysis.(analysisParams.field).zScore,3));
   for stim = 1:size(analysis.(analysisParams.field).zScore,3)
       for trial=1:size(analysis.(analysisParams.field).zScore,4)
          %randomly drawing pixels from the corresponding area corresponding 
          %to how many pixels were active in that area at that trial 
          randA19pxlNumbers = randperm(numel(A19Pixels), actPxA19(stim,trial));
          randA19pxls = A19Pixels(randA19pxlNumbers);
          randV1pxlNumbers = randperm(numel(V1Pixels), actPxV1(stim,trial));
          randV1pxls = V1Pixels(randV1pxlNumbers);
          randAllpxls = [randA19pxls' randV1pxls'];
          
          %put the pixel together to a map
          randMap = zeros(size(analysis.(analysisParams.field).zScore,1)*size(analysis.(analysisParams.field).zScore,2),1);
          randMap(randAllpxls) = 1;
          randMap = reshape(randMap,size(analysis.(analysisParams.field).zScore,1), size(analysis.(analysisParams.field).zScore,2));
          
          %add the individual maps on top of each other
          overlapAllRnd(:,:,stim) = overlapAllRnd(:,:,stim)+randMap;
       end
       %divide them back into areas
       boundaryAreaA19Rnd(:,:,stim) = overlapAllRnd(:,:,stim).*analysis.maskA19;
       boundaryAreaV1Rnd(:,:,stim) = overlapAllRnd(:,:,stim).*analysis.maskV1;
   end
   %to count pixel, make sure that all of them are 0 (not used at all) or 1
   %(used at least once)
   boundaryAreaV1Rnd(boundaryAreaV1Rnd > 0) = 1;
   boundaryAreaA19Rnd(boundaryAreaA19Rnd > 0) = 1;
   
   %how many pixels are being used in total?
   SumBoundAreaA19RndAll(drawing,:) = squeeze(nansum(boundaryAreaA19Rnd,[1 2]));
   SumBoundAreaV1RndAll(drawing,:) = squeeze(nansum(boundaryAreaV1Rnd,[1 2]));
end

%get the average
SumBoundAreaA19Rnd=median(SumBoundAreaA19RndAll,1)';
SumBoundAreaV1Rnd=median(SumBoundAreaV1RndAll,1)';

%convergence Factor: how much do the active pixel converge together?
% = 1 - Percentage of actual pixels used of possible pixel used
% the higher, the more the areas are converging
analysis.(analysisParams.field).convFactA19 = 1-(SumBoundAreaA19./SumBoundAreaA19Rnd);
analysis.(analysisParams.field).convFactV1 = 1-(SumBoundAreaV1./SumBoundAreaV1Rnd);

%should you also take into account the size of the area and how much of it
%is being used?
%convFactV1 = ((SumBoundAreaV1./SumBoundAreaV1Rnd)).*SumBoundAreaV1Rnd/length(V1Pixels);


