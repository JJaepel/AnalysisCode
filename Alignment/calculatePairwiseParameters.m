function pwMeasures = calculatePairwiseParameters(Spines)

% Calculates pairwise measurements for spines on the same dendritic segment
%
% Input:
% - Spines: all ROIs with their properties

% Steps:
% 1.) Define structure
% 2.) Go through each dendrite 
% 3.) Do all pairwise measurements
% 
% Output: 
% - pwMeasures: a structure containing measurement differences between
% pairs for differenct parameters, each row is one pair, measuring:
%    - if both spines are responsive
%    - if both spines are selective (OSI/DSIvect > 0.1)]
%    - distance between them
%    - deltaOri
%    - deltaDir
%    - deltaOSI
%    - deltaDSI
%    - deltaDSIvect
%    - tuningCurve correation
%    - same field of view?
%    - trial-to-trial correlation

% Written by Juliane Jaepel
% Max Planck Florida Institute for Neuroscience
% Version 1.0: February, 2024

%restructure as following: go dendrite by dendrite, find all ROIs on their,
%for all the pairs calculate distances, deltaOri etc, then multiply their
%responsiveness good x good, as well as their logical osi (> 0.1 yes/no) as
%well as their logical dsi to get ori-selective pairs as well as
%dsi-selective pairs

%% 1.) Define the structure
pwMeasures = struct;

%% 2.) Find all dendrites with ROIs and loop through them
%get all dendrites that have at least one ROI
dendNrs = unique([Spines.Dendrite]);

%start the counter for the row
pairCounter = 1;

%go through all those dendrites
for d = 1:length(dendNrs)
    %find all spines on that dendrite
    spineOnDend = find([Spines.Dendrite] == dendNrs(d));
    if length(spineOnDend) > 1 %if there is more than one spine
        %nested loop to generate the pairs
        for i = 1:numel(spineOnDend)
            for j = i+1:numel(spineOnDend)  
                %only add them if there are 2p matched
                if Spines(spineOnDend(i)).TwoPMatch && Spines(spineOnDend(j)).TwoPMatch
                    %general info              
                    pwMeasures(pairCounter).Branch =Spines(spineOnDend(i)).Branch;
                    pwMeasures(pairCounter).Dendrite = dendNrs(d);
                    pwMeasures(pairCounter).BranchOrder = Spines(spineOnDend(i)).BranchOrder;
                    pwMeasures(pairCounter).type = Spines(spineOnDend(i)).type;
                    pwMeasures(pairCounter).spineA = Spines(spineOnDend(i)).Nr;
                    pwMeasures(pairCounter).spineB = Spines(spineOnDend(j)).Nr;

                    %distance
                    pwMeasures(pairCounter).distance = abs([Spines(spineOnDend(i)).distToSoma]-[Spines(spineOnDend(j)).distToSoma]');
                    %responsive
                    pwMeasures(pairCounter).goodPair = Spines(spineOnDend(i)).good * Spines(spineOnDend(j)).good;
                    
                    %selective for ori and dir (good pair as well as
                    %OSI/DSI for both higher than 0.1
                    pwMeasures(pairCounter).oriSelect = (Spines(spineOnDend(i)).OSI > 0.1 && Spines(spineOnDend(j)).OSI)* pwMeasures(pairCounter).goodPair;
                    pwMeasures(pairCounter).dirSelect = (Spines(spineOnDend(i)).DSIvect > 0.1 && Spines(spineOnDend(j)).DSIvect)*pwMeasures(pairCounter).goodPair;
                    
                    %deltaOri
                    pwMeasures(pairCounter).deltaOri = abs([Spines(spineOnDend(i)).prefOri] - [Spines(spineOnDend(j)).prefOri]');
                    if pwMeasures(pairCounter).deltaOri>90
                       pwMeasures(pairCounter).deltaOri = 180-pwMeasures(pairCounter).deltaOri;
                    end
                    %deltaDir
                    pwMeasures(pairCounter).deltaDir = abs([Spines(spineOnDend(i)).prefDir] - [Spines(spineOnDend(j)).prefDir]');
                    if pwMeasures(pairCounter).deltaDir>180
                       pwMeasures(pairCounter).deltaDir = 360-pwMeasures(pairCounter).deltaDir;
                    end

                    %deltaOSI
                    pwMeasures(pairCounter).deltaOSI = abs([Spines(spineOnDend(i)).OSI] - [Spines(spineOnDend(j)).OSI]');

                    %deltaDSI
                    pwMeasures(pairCounter).deltaDSI = abs([Spines(spineOnDend(i)).DSI] - [Spines(spineOnDend(j)).DSI]');

                    %deltaDSIvect
                    pwMeasures(pairCounter).deltaDSIvect = abs([Spines(spineOnDend(i)).DSIvect] - [Spines(spineOnDend(j)).DSIvect]');
                    
                    %tuningCurve correlation
                    pwMeasures(pairCounter).curveCorr = corr(Spines(spineOnDend(i)).meanResp', Spines(spineOnDend(j)).meanResp');
                    
                    %check if they are in the same field of view
                    if Spines(spineOnDend(i)).expNr == Spines(spineOnDend(i)).expNr
                        pwMeasures(pairCounter).sameExp = 1;
                        
                        %also do trial to trial correlation
                        pwMeasures(pairCounter).trialCorr = trialToTrialCorrelation(Spines(spineOnDend(i)).cycRes, Spines(spineOnDend(j)).cycRes);
                    else
                        pwMeasures(pairCounter).sameExp = 0;
                        pwMeasures(pairCounter).trialCorr = [];
                    end
                    
                    pairCounter = pairCounter + 1;
                end
                
            end
        end
    end
    
end
