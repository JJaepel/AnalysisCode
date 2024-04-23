function [localDeltaOriRandom, localDeltaOriSomaRandom, localDispValuesRandom] = randomizeInputPositionForlocalMeasurements(allSpines, allBranches, cellDetails, inputType)



%% Step 1: Find the ori-selective good inputs
%switch depending on the inputType
switch inputType
    case 'A19'
        inputROIs = find([allSpines.input] == 1);
    case 'V1'
        inputROIs = find([allSpines.input] == 2);
end

%make selectors for good and ori-selective spines
%which ones are 2p spines
TwoPROIsNR = find([allSpines.TwoPMatch] ==1);
%which of those are responsive
goodROIs = find([allSpines.good] == 1);
%which of those are selective
oriTwoPROIs = find([allSpines.OSI] > 0.1); %all oriselectROIs, size = allfuncSpines
oriSelect = TwoPROIsNR(oriTwoPROIs);
oriGood = intersect(oriSelect, goodROIs);

%intersect with inputs
oriGoodInputs = intersect(oriGood, inputROIs);

%% Step 2: Prelocate variables to save the data in - each input is resampled 1000 times
localDeltaOriRandom = NaN(1000, length(oriGoodInputs));
localDeltaOriSomaRandom = NaN(1000, length(oriGoodInputs));
localDispValuesRandom = NaN(1000, length(oriGoodInputs));

%% Step 3: Go through all inputs after each other, find the branch, randomly assign a position and exchange its property with the ROI that was on there, calculate localDeltaOri etc. and save it

for i = 1:length(oriGoodInputs)
    dendriteOnCell = find([allBranches.CellNr] == allSpines(oriGoodInputs(i)).CellNr);
    sameDendriteNr = find([allBranches.dendNr] == allSpines(oriGoodInputs(i)).Dendrite);
    dendriteID = intersect(dendriteOnCell, sameDendriteNr);
    
    if size(allBranches(dendriteID).ROIs,2) > 2 %should be at least 3 spines on there
       dendriteROIs = allBranches(dendriteID).ROIs;
       %get the input ID and save its data separately
       inputID = find([dendriteROIs.ROINrOnBranch] == allSpines(oriGoodInputs(i)).ROINrOnBranch);
       
       if ~isempty(inputID)
           inputIDData = dendriteROIs(inputID);

           %what are the good spines on the soma?
           goodSpinesOnSameDend = intersect(find([dendriteROIs.OSI] > 0.1), find([dendriteROIs.good] == 1));

           %make sure there are at least 3 in total
           if length(goodSpinesOnSameDend) > 2
               %get soma preference
               SomaPref = cellDetails(allSpines(oriGoodInputs(i)).CellNr).oriPref;
               
               %now randomize its position 1000times
               for r = 1:1000
                   temp = dendriteROIs;
                   %randomize the position, using one that is also a good input             
                   newInputPosOrder = randperm((length(goodSpinesOnSameDend)));
                   newInputPos = goodSpinesOnSameDend(newInputPosOrder(1));

                   %change the distance to Soma between the new input and the old
                   %input
                   try
                      temp(inputID).distToSoma = temp(newInputPos).distToSoma;
                   catch
                      disp(['i = ' num2str(i) ' , new Pos =  ' num2str(newInputPos)])
                   end
                   temp(newInputPos).distToSoma = inputIDData.distToSoma;

                   %now find the ones that are within 2.5/5 um of its position
                   allDist = abs(temp(newInputPos).distToSoma - [temp(goodSpinesOnSameDend).distToSoma]);
                   withinRadius = find(allDist < 5);

                   %make sure that it is at last 3 in total
                   if length(withinRadius) > 2
                        %1) Local dispersion 
                        prefOri = [dendriteROIs(goodSpinesOnSameDend(withinRadius)).prefOri];
                        Theta = (prefOri/360)*2*pi;
                        %get the mean theta of the local environment
                        meanTheta = circ_mean(2*Theta')/2;
                        if meanTheta < 0
                            meanTheta = meanTheta + pi;
                        end
                        %calculate mean distance of all spines from that angle
                        localDispersion = (1/(length(2*Theta)))*sum(pi-abs(pi-abs(2*Theta-2*meanTheta)))/2;
                        localDispValuesRandom(r,i) = rad2deg(localDispersion);


                        % 2) deltaOri to Soma
                        meanOri = rad2deg(meanTheta);
                        localDeltaOriSomaRandom(r,i) = abs(meanOri-SomaPref);
                        if localDeltaOriSomaRandom(r,i) > 90
                           localDeltaOriSomaRandom(r,i) = 180-localDeltaOriSomaRandom(r,i);
                        end

                        % 3) deltaOri to spine
                        %what is the deltaOri to the spine
                        try
                        deltaOris = abs(prefOri-temp(newInputPos).prefOri);
                        catch
                            disp(['i = ' num2str(i) ' , new Pos =  ' num2str(newInputPos)])
                        end
                        deltaOris(deltaOris>90) = 180-deltaOris(deltaOris>90);
                        %what is the mean deltaOri?
                        localDeltaOriRandom(r,i) = mean(deltaOris);
                   end
               end
           end
       end
    end
end

%% Step 4: Reformat the vector and remove the nans

%reshape
localDeltaOriRandom = reshape(localDeltaOriRandom, 1,1000*length(oriGoodInputs));
localDeltaOriSomaRandom = reshape(localDeltaOriSomaRandom, 1,1000*length(oriGoodInputs));
localDispValuesRandom = reshape(localDispValuesRandom, 1,1000*length(oriGoodInputs));

%remove
localDeltaOriRandom(isnan(localDeltaOriRandom)) = [];
localDeltaOriSomaRandom(isnan(localDeltaOriSomaRandom)) = [];
localDispValuesRandom(isnan(localDispValuesRandom)) = [];
