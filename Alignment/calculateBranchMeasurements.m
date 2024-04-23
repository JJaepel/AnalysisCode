function Dendrites = calculateBranchMeasurements(Dendrites, Spines, Soma)

% Calculates various measurements across the dendritic segment
%
% Input:
% - Spines: all ROIs with their properties
% - Dendrites: all branches with their properties
% 
% Output: 
% - Dendrites: all branches with the new property circular dispersion &
% mean orientation
%
% Criteria for including Spines/Dendrites
% - spines are responsive
% -[ spines are selective (OSI/DSIvect > 0.1)]
% - at least 3 per branch

%add a new column to the Dendrites
Dendrites = arrayfun(@(x) setfield(x, 'circDispersion', NaN(1,1)), Dendrites); 
Dendrites = arrayfun(@(x) setfield(x, 'meanOri', NaN(1,1)), Dendrites); 
Dendrites = arrayfun(@(x) setfield(x, 'deltaOri', NaN(1,1)), Dendrites);
Dendrites = arrayfun(@(x) setfield(x, 'dirCircDispersion', NaN(1,1)), Dendrites); 
Dendrites = arrayfun(@(x) setfield(x, 'meanDir', NaN(1,1)), Dendrites); 
Dendrites = arrayfun(@(x) setfield(x, 'deltaDir', NaN(1,1)),Dendrites);
Dendrites = arrayfun(@(x) setfield(x, 'medianOSI', NaN(1,1)), Dendrites); 
Dendrites = arrayfun(@(x) setfield(x, 'medianDSI', NaN(1,1)), Dendrites);
Dendrites = arrayfun(@(x) setfield(x, 'medianDSIvect', NaN(1,1)), Dendrites);
% Dendrites = arrayfun(@(x) setfield(x, 'pwDistOriSelect', NaN(1,1)), Dendrites); 
% Dendrites = arrayfun(@(x) setfield(x, 'pwDistDirSelect', NaN(1,1)), Dendrites);
% Dendrites = arrayfun(@(x) setfield(x, 'pwDistAll', NaN(1,1)), Dendrites);
% Dendrites = arrayfun(@(x) setfield(x, 'pwDeltaOri', NaN(1,1)), Dendrites); 
% Dendrites = arrayfun(@(x) setfield(x, 'pwDeltaDir', NaN(1,1)), Dendrites);
% Dendrites = arrayfun(@(x) setfield(x, 'pwDeltaOSI', NaN(1,1)), Dendrites);
% Dendrites = arrayfun(@(x) setfield(x, 'pwDeltaDSI', NaN(1,1)), Dendrites);
% Dendrites = arrayfun(@(x) setfield(x, 'pwDeltaDSIvect', NaN(1,1)), Dendrites);

%make sure we are only looking at the good spines
goodSpines = find([Spines.good] == 1);
selectiveSpines = find([Spines.OSI] > 0.1);
selectiveDirSpines = find([Spines.DSIvect] > 0.1);

%now go through all dendrites
for d = 1:length(Dendrites)
    %select spines on the dendrite
    spineOnDend = find([Spines.Dendrite] == d); %find all spines on that dend
    goodSpinesOnDend = intersect(goodSpines, spineOnDend); %make sure they are responsive
    goodOriSelecSpinesOnDend = intersect(goodSpinesOnDend, selectiveSpines); %make sure they are selective
    goodDirSelecSpinesOnDend = intersect(goodSpinesOnDend, selectiveDirSpines); %make sure they are selective
    
    %circularDispersion & meanOri
    if length(goodOriSelecSpinesOnDend)>2
        %get the ori pref of all the spines
        prefOri =[Spines(goodOriSelecSpinesOnDend).prefOri];
        %get the theta
        Theta = (prefOri/360)*2*pi;

        %get the orientation of the dendrite
        meanTheta = circ_mean(2*Theta')/2;
        if meanTheta < 0
            meanTheta = meanTheta+pi;
        end

        %calculate mean distance of all spines from that angle
        circDispersion = (1/(length(2*Theta)))*sum(pi-abs(pi-abs(2*Theta-2*meanTheta)))/2;
        circDispersionDeg = rad2deg(circDispersion);

        %write it to the array
        Dendrites(d).circDispersion = circDispersionDeg;
        Dendrites(d).meanOri =  rad2deg(meanTheta);
        
%         %get the pairwise distance
%         pwDist = abs([Spines(goodOriSelecSpinesOnDend).distToSoma]-[Spines(goodOriSelecSpinesOnDend).distToSoma]');
%         m = triu(true(size(pwDist)),1); %get only the values of the upper triangle without the diagonal
%         Dendrites(d).pwDistOriSelect = pwDist(m);
%         
%         %get the pairwise difference in prefOri
%         pwDeltaOri = abs([Spines(goodOriSelecSpinesOnDend).prefOri] - [Spines(goodOriSelecSpinesOnDend).prefOri]');
%         pwDeltaOri(pwDeltaOri > 90) = 180-pwDeltaOri(pwDeltaOri>90);
%         Dendrites(d).pwDeltaOri = pwDeltaOri(m);
        
    else
        Dendrites(d).circDispersion = NaN;
        Dendrites(d).meanOri =NaN;
    end
    
    %deltaOri to Soma
    Dendrites(d).deltaOri = abs(Dendrites(d).meanOri-Soma.prefOri);
    if Dendrites(d).deltaOri > 90
        Dendrites(d).deltaOri = 180-Dendrites(d).deltaOri;
    end
    
    %dirCircDispersion, meanDir, pairwise distance prefDir
    if length(goodDirSelecSpinesOnDend)>2
        %get the dir pref of all the spines
        prefDir =[Spines(goodDirSelecSpinesOnDend).prefDir];
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
        
%         %get the pairwise distance
%         pwDist = abs([Spines(goodDirSelecSpinesOnDend).distToSoma]-[Spines(goodDirSelecSpinesOnDend).distToSoma]');
%         m = triu(true(size(pwDist)),1); %get only the values of the upper triangle without the diagonal
%         Dendrites(d).pwDistDirSelect = pwDist(m);
%         
%         %get the pairwise difference in prefDir
%         pwDeltaDir = abs([Spines(goodDirSelecSpinesOnDend).prefDir] - [Spines(goodDirSelecSpinesOnDend).prefDir]');
%         pwDeltaDir(pwDeltaDir > 180) = 360-pwDeltaDir(pwDeltaDir>180);
%         Dendrites(d).pwDeltaDir = pwDeltaDir(m);
    else
        Dendrites(d).dirCircDispersion = NaN;
        Dendrites(d).meanDir =NaN;
    end
    
    %deltaDir to Soma
    Dendrites(d).deltaDir = abs(Dendrites(d).meanDir-Soma.prefDir);
    if Dendrites(d).deltaDir > 180
        Dendrites(d).deltaDir = 360-Dendrites(d).deltaDir;
    end
    
    %Branch selectivity
    if length(goodSpinesOnDend)>2
        %median selectivity
        Dendrites(d).medianOSI = median([Spines(goodSpinesOnDend).OSI]);
        Dendrites(d).medianDSI = median([Spines(goodSpinesOnDend).DSI]);
        Dendrites(d).medianDSIvect = median([Spines(goodSpinesOnDend).DSIvect]);
        
%         %get the pairwise distance
%         pwDist = abs([Spines(goodSpinesOnDend).distToSoma]-[Spines(goodSpinesOnDend).distToSoma]');
%         m = triu(true(size(pwDist)),1); %get only the values of the upper triangle without the diagonal
%         Dendrites(d).pwDistAll = pwDist(m);
%         
%         %get the pairwise difference in selectivity markers
%         pwDeltaOSI = abs([Spines(goodDirSelecSpinesOnDend).OSI] - [Spines(goodDirSelecSpinesOnDend).OSI]');
%         Dendrites(d).pwDeltaOSI = pwDeltaOSI(m);
%         pwDeltaDSI = abs([Spines(goodDirSelecSpinesOnDend).DSI] - [Spines(goodDirSelecSpinesOnDend).DSI]');
%         Dendrites(d).pwDeltaDSI = pwDeltaDSI(m);
%         pwDeltaDSIvect = abs([Spines(goodDirSelecSpinesOnDend).DSIvect] - [Spines(goodDirSelecSpinesOnDend).DSIvect]');
%         Dendrites(d).pwDeltaDSIvect = pwDeltaDSIvect(m);
    end
    
    
end

