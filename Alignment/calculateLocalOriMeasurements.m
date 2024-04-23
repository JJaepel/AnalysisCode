function Spines = calculateLocalOriMeasurements(Spines, Soma)

% Calculates local orientation measurements within different radiuses:
% - local dispersion: measure of how homeogeneous the local environment is
% -> how many degree difference in terms of std
% - homeogeneity index: similar measurement, going from 0 to 1
% - localDeltaOri: what is the deltaOri to other spines depend on the
% radius?
% - localDeltaOriSoma: how different is the local environment from the
% soma?
% - nearestSimilarPrefOri: what is the distance to the nearest spine with a
% similar ori pref (deltaOri <=22.5 deg) - NaN if there is none
% - how selective the environment is -> median OSI
%
% Input:
% - Spines: all ROIs with their properties
% - Soma: Soma with its properties
% 
% Output: 
% - Spines: all ROIs with the new measurements
%
% Criteria for including Spines
% - spines are responsive
% - within a certain radius or distance
% - spines are selective (OSI > 0.1) -> not for selectivity measurement

%which dist and radius do you want to look at?
radius = [2.5, 5, 7.5, 10]; %within 2.5/5/7.5/10
dist = [0, 5, 10, 15]; %between 0-5/5-10/10-15/more than 15

% 0) Add new columns to the Spines
Spines = arrayfun(@(x) setfield(x, 'localDispersion', NaN(1,4)), Spines);
Spines = arrayfun(@(x) setfield(x, 'HI', NaN(1,4)), Spines); 
Spines = arrayfun(@(x) setfield(x, 'localDeltaOri', NaN(1,4)), Spines);
Spines = arrayfun(@(x) setfield(x, 'localDeltaOriSoma', NaN(1,4)), Spines);
Spines = arrayfun(@(x) setfield(x, 'nearestSimilarPrefOri', NaN(1,1)), Spines);
Spines = arrayfun(@(x) setfield(x, 'localOSI', NaN(4,1)), Spines);

%make sure we are only looking at the good spines that are selective
goodSpines = find([Spines.good] == 1);
selectiveSpines = find([Spines.OSI] > 0.1);
TwoPNr = find([Spines.TwoPMatch] == 1);
oriSelectSpines = TwoPNr(selectiveSpines); 
goodSelectSpines = intersect(goodSpines,oriSelectSpines);

%% I. Local dispersion, homeogenetiy & deltaOri vs. distance
%now go through all selective spines for homeogeneity index, and local
%dispersion
for gSS = 1:length(goodSelectSpines)
    % 1) Select spines on the same dendrite
    s = goodSelectSpines(gSS);
    dendrNr = Spines(s).Dendrite; %get the dendrite Nr of the spine
    spineOnSameDend = find([Spines.Dendrite] == dendrNr); %find all spines on that dend
    goodSpinesOnSameDend = intersect(goodSelectSpines, spineOnSameDend); %make sure they are also good & selective
    
    % 2) Get dist to all Spines
    allDist = abs(Spines(s).distToSoma - [Spines(goodSpinesOnSameDend).distToSoma]);
    
    % 3) For each radius, calculate the measurements
    for r = 1:length(radius)
        withinRadius = find(allDist < radius(r));
        if length(withinRadius) > 2
            
            % 3a) Local dispersion 
            prefOri = [Spines(goodSpinesOnSameDend(withinRadius)).prefOri];
            Theta = (prefOri/360)*2*pi;
            
            %get the mean theta of the local environment
            meanTheta = circ_mean(2*Theta')/2;
            if meanTheta < 0
                meanTheta = meanTheta + pi;
            end
            
            %calculate mean distance of all spines from that angle
            localDispersion = (1/(length(2*Theta)))*sum(pi-abs(pi-abs(2*Theta-2*meanTheta)))/2;
            localDisp(r) = rad2deg(localDispersion);
            
            % 3b) Homeogeneity index
            thresholdOrientations=exp(2*1i*prefOri);
            HI(r) = 1-(abs(1/length(withinRadius) * sum(thresholdOrientations)));
            if HI(r) == 0
                HI(r) = NaN;
            end
            
            % 3c) deltaOri to Soma
            meanOri = rad2deg(meanTheta);
            deltaOriSoma(r) = abs(meanOri-Soma.prefOri);
            if deltaOriSoma(r) > 90
                deltaOriSoma(r) = 180-deltaOriSoma(r);
            end
        else
            localDisp(r) = NaN;
            HI(r) = NaN;
            deltaOriSoma(r) = NaN;
        end
    end
    
    % 4) For each distance (0-5, 5-10, 10-15, >15, calculate the delta Ori
    spineOri = Spines(s).prefOri;
    for d = 1:length(dist)
        if d < length(dist)
           %which ones are between the two distances?
           withinRadius = find(allDist >= dist(d));
           withinRadius2 = find(allDist < dist(d+1));
           allWithin = intersect(withinRadius, withinRadius2);
           
           if length(allWithin) > 2
               %what is their prefOri?
               allOris = [Spines(goodSpinesOnSameDend(allWithin)).prefOri];
               %what is the deltaOri to the spine
               deltaOris = abs(allOris-spineOri);
               deltaOris(deltaOris>90) = 180-deltaOris(deltaOris>90);
               %what is the mean deltaOri?
               localDeltaOri(d) = mean(deltaOris);
           else
               localDeltaOri(d) = NaN;
           end
        else
            %which ones are greater than the distance
            farAway = find(allDist > dist(d));
            if length(farAway) > 2
               %what is their prefOri?
               allOris = [Spines(goodSpinesOnSameDend(farAway)).prefOri];
               %what is the deltaOri to the spine
               deltaOris = abs(allOris-spineOri);
               deltaOris(deltaOris>90) = 180-deltaOris(deltaOris>90);
               %what is the mean deltaOri?
               localDeltaOri(d) = mean(deltaOris);
           else
               localDeltaOri(d) = NaN;
           end
            
        end
    end
    
    % 5) Find the nearest neighbor with a deltaOri < 25
    %get all Oris on the segment
    allOrisOnDend = [Spines(goodSpinesOnSameDend).prefOri];
    %get the difference to the spine
    allDeltaOri = abs(allOrisOnDend-spineOri);
    allDeltaOri(allDeltaOri>90) = 180-allDeltaOri(allDeltaOri>90);
    %set the spine itself to NaN, otherwise distance will be 0
    allDeltaOri(find(goodSpinesOnSameDend == s)) = NaN;
    %find the ones that are below 22.5
    similarSpines = find(allDeltaOri <=22.5);
    if ~isempty(similarSpines)
        distSimSpines = allDist(similarSpines);
        Spines(s).nearestSimilarPrefOri = min(distSimSpines);
    else
        Spines(s).nearestSimilarPrefOri = NaN;
    end
    
    % 6) Save to structure
    Spines(s).localDispersion = localDisp;
    Spines(s).HI = HI;
    Spines(s).localDeltaOriSoma = deltaOriSoma;
    Spines(s).localDeltaOri = localDeltaOri;
end

%% II. Mean OSI of environment
%now go through all good spines for median OSI
for gs = 1:length(goodSpines)
    % 1) Select spines on the same dendrite
    sp = goodSpines(gs);
    dendrNr = Spines(sp).Dendrite; %get the dendrite Nr of the spine
    spineOnSameDend = find([Spines.Dendrite] == dendrNr); %find all spines on that dend
    goodSpinesOnSameDend = intersect(goodSpines, spineOnSameDend); %make sure they are also good
    
     % 2) Get dist to all Spines
    allDist = abs(Spines(sp).distToSoma - [Spines(goodSpinesOnSameDend).distToSoma]);
    
    % 3) For each radius, calculate the measurements
    for r = 1:length(radius)
        withinRadius = find(allDist < radius(r));
        if length(withinRadius) > 2
            allOSI = [Spines(goodSpinesOnSameDend(withinRadius)).OSI];
            localOSI(r) = nanmean(allOSI);
        else
            localOSI(r) = NaN;
        end
    end
    
     % 4) Save to structure
     Spines(sp).localOSI = localOSI;
end
