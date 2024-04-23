function Dendrites = calculateBranchSelectivity(Dendrites, Spines)

% Calculates median selectivity indices for branches

% Input:
% - Spines: all ROIs with their properties
% - Dendrites: all branches with their properties
% 
% Output: 
% - Dendrites: all branches with the new properties medianOSI, medianDSI &
% medianDSIvect
%
% Criteria for including Spines/Dendrites
% - spines are responsive
% - at least 3 per branch

%add a new column to the Dendrites
Dendrites = arrayfun(@(x) setfield(x, 'medianOSI', NaN(1,1)), Dendrites); 
Dendrites = arrayfun(@(x) setfield(x, 'medianDSI', NaN(1,1)), Dendrites);
Dendrites = arrayfun(@(x) setfield(x, 'medianDSIvect', NaN(1,1)), Dendrites);

%make sure we are only looking at the good spines
goodSpines = find([Spines.good] == 1);

for d = 1:length(Dendrites)
    %select spines on the dendrite
    spineOnDend = find([Spines.Dendrite] == d); %find all spines on that dend
    goodSpinesOnDend = intersect(goodSpines, spineOnDend); %make sure they are responsive
    
    if length(goodSpinesOnDend)>2
        Dendrites(d).medianOSI = median([Spines(goodSpinesOnDend).OSI]);
        Dendrites(d).medianDSI = median([Spines(goodSpinesOnDend).DSI]);
        Dendrites(d).medianDSIvect = median([Spines(goodSpinesOnDend).DSIvect]);
    end
end