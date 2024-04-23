function Spines = calculateLocalDirMeasurements(Spines, Soma)

% Calculates local direction measurements within different radiuses:
% - local dispersion: measure of how homeogeneous the local environment is
% -> how many degree difference in terms of std
% - homeogeneity index: similar measurement, going from 0 to 1
% - localDeltaDir: what is the deltaDir to other spines depend on the
% radius?
% - localDeltaDirSoma: how different is the local environment from the
% soma?
% - nearestSimilarPrefDir: what is the distance to the nearest spine with a
% similar ori pref (deltaDir <=45 deg) - NaN if there is none
% - how selective the environment is -> median DSIvect
%
% Input:
% - Spines: all ROIs with their properties
% 
% Output: 
% - Spines: all ROIs with the new measurements
%
% Criteria for including Spines
% - spines are responsive
% - within a certain radius or distance
% - spines are selective (DSI > 0.1) -> not for selectivity measurement

%which dist and radius do you want to look at?
radius = [2.5, 5, 7.5, 10]; %within 2.5/5/7.5/10
dist = [0, 5, 10, 15]; %between 0-5/5-10/10-15/more than 15

% 0) Add new columns to the Spines
Spines = arrayfun(@(x) setfield(x, 'localDirDispersion', NaN(1,4)), Spines);
Spines = arrayfun(@(x) setfield(x, 'HIDir', NaN(1,4)), Spines); 
Spines = arrayfun(@(x) setfield(x, 'localDeltaDir', NaN(1,4)), Spines);
Spines = arrayfun(@(x) setfield(x, 'localDeltaDirSoma', NaN(1,4)), Spines);
Spines = arrayfun(@(x) setfield(x, 'nearestSimilarPrefDir', NaN(1,1)), Spines);
Spines = arrayfun(@(x) setfield(x, 'localDSI', NaN(4,1)), Spines);

%make sure we are only looking at the good spines that are selective
goodSpines = find([Spines.good] == 1);
selectiveSpines = find([Spines.DSI] > 0.1);
TwoPNr = find([Spines.TwoPMatch] == 1);
dirSelectSpines = TwoPNr(selectiveSpines); 
goodSelectSpines = intersect(goodSpines,dirSelectSpines);

%% I. Local dispersion, homeogenetiy & deltaDir vs. distance
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
            prefDirs = [Spines(goodSpinesOnSameDend(withinRadius)).prefDir];
            Theta = (prefDirs/360)*2*pi;
            
            %get the mean angle of the local environment
            meanTheta = circ_mean(Theta');
            
            %calculate mean distance of all spines from that angle
            localDispersion = (1/(length(Theta)))*sum(pi-abs(pi-abs(Theta-meanTheta)));
            localDisp(r) = rad2deg(localDispersion);
            
            % 3b) Homeogeneity index
            thresholdDirs=exp(1i*prefDirs);
            HI(r) = 1-(abs(1/length(withinRadius) * sum(thresholdDirs)));
            if HI(r) == 0
                HI(r) = NaN;
            end
            
             % 3c) deltaDir to Soma
            meanDir = rad2deg(meanTheta);
            deltaDirSoma(r) = abs(meanDir-Soma.prefDir);
            if deltaDirSoma(r) > 180
                deltaDirSoma(r) = 360-deltaDirSoma(r);
            end
        else
            localDisp(r) = NaN;
            HI(r) = NaN;
            deltaDirSoma(r) = NaN;
        end
    end
    
    % 4) For each distance (0-5, 5-10, 10-15, >15, calculate the delta Ori
    spineDir = Spines(s).prefDir;
    for d = 1:length(dist)
        if d < length(dist)
           %which ones are between the two distances?
           withinRadius = find(allDist >= dist(d));
           withinRadius2 = find(allDist < dist(d+1));
           allWithin = intersect(withinRadius, withinRadius2);
           
           if length(allWithin) > 2
               %what is their prefOri?
               allDirs = [Spines(goodSpinesOnSameDend(allWithin)).prefDir];
               %what is the deltaOri to the spine
               deltaDirs = abs(allDirs-spineDir);
               deltaDirs(deltaDirs>180) = 360-deltaDirs(deltaDirs>180);
               %what is the mean deltaOri?
               localDeltaDir(d) = mean(deltaDirs);
           else
               localDeltaDir(d) = NaN;
           end
        else
            %which ones are greater than the distance
            farAway = find(allDist > dist(d));
            if length(farAway) > 2
               %what is their prefOri?
               allDirs = [Spines(goodSpinesOnSameDend(farAway)).prefDir];
               %what is the deltaOri to the spine
               deltaDirs = abs(allDirs-spineDir);
               deltaDirs(deltaDirs>180) = 360-deltaDirs(deltaDirs>180);
               %what is the mean deltaOri?
               localDeltaDir(d) = mean(deltaDirs);
           else
               localDeltaDir(d) = NaN;
           end
            
        end
    end
    
    % 5) Find the nearest neighbor with a deltaDir < 45
    %get all Oris on the segment
    allDirsOnDend = [Spines(goodSpinesOnSameDend).prefDir];
    %get the difference to the spine
    allDeltaDir = abs(allDirsOnDend-spineDir);
    allDeltaDir(allDeltaDir>180) = 360-allDeltaDir(allDeltaDir>180);
    %set the spine itself to NaN, otherwise distance will be 0
    allDeltaDir(find(goodSpinesOnSameDend == s)) = NaN;
    %find the ones that are below 45
    similarSpines = find(allDeltaDir <=45);
    if ~isempty(similarSpines)
        distSimSpines = allDist(similarSpines);
        Spines(s).nearestSimilarPrefDir = min(distSimSpines);
    else
        Spines(s).nearestSimilarPrefDir = NaN;
    end
    
    % 6) Save to structure
    Spines(s).localDirDispersion = localDisp;
    Spines(s).HIDir = HI;
    Spines(s).localDeltaDirSoma = deltaDirSoma;
    Spines(s).localDeltaDir = localDeltaDir;
end

%% II. Mean DSI of environment
%now go through all good spines for median DSI
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
            allDSI = [Spines(goodSpinesOnSameDend(withinRadius)).DSIvect];
            localDSI(r) = mean(allDSI);
        else
            localDSI(r) = NaN;
        end
    end
    
     % 4) Save to structure
     Spines(sp).localDSI = localDSI;
end
