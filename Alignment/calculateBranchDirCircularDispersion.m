function Dendrites = calculateBranchMeasurements(Dendrites, Spines)

% Caluculates dir circular dispersion for branches as a measure of how
% homeogeneous a branch is -> this is in direction space

% Input:
% - Spines: all ROIs with their properties
% - Dendrites: all branches with their properties
% 
% Output: 
% - Dendrites: all branches with the new property circular dispersion
%
% Criteria for including Spines/Dendrites
% - spines are responsive
% - spines are selective (DSI > 0.1)
% - at least 3 per branch


%add a new column to the Dendrites
Dendrites = arrayfun(@(x) setfield(x, 'dirCircDispersion', NaN(1,1)), Dendrites); 
Dendrites = arrayfun(@(x) setfield(x, 'meanDir', NaN(1,1)), Dendrites); 

%make sure we are only looking at the good spines
goodSpines = find([Spines.good] == 1);
selectiveSpines = find([Spines.DSI] > 0.1);

%now go through all dendrites
for d = 1:length(Dendrites)
    %select spines on the dendrite
    spineOnDend = find([Spines.Dendrite] == d); %find all spines on that dend
    goodSpinesOnDend = intersect(goodSpines, spineOnDend); %make sure they are responsive
    goodSelecSpinesOnDend = intersect(goodSpinesOnDend, selectiveSpines); %make sure they are selective
    
    if length(goodSelecSpinesOnDend)>2
        %get the dir pref of all the spines
        prefDir =[Spines(goodSpinesOnDend).prefDir];
        %get the theta
        Theta = (prefDir/360)*2*pi;
        
        %get the direction of the dendrite
        meanTheta = circ_mean(Theta');
        if meanTheta < 0
            meanTheta = meanTheta + 2*pi;
        end

        %calculate mean distance of all spines from that angle
        dirCircDispersion = (1/(length(Theta)))*sum(pi-abs(pi-abs(Theta-meanTheta)));
        dirCircDispersionDeg = rad2deg(dirCircDispersion);

        %write it to the array
        Dendrites(d).dirCircDispersion = dirCircDispersionDeg;
        Dendrites(d).meanDir = rad2deg(meanTheta);
    else
        Dendrites(d).dirCircDispersion = NaN;
        Dendrites(d).meanDir =NaN;
    end
end

