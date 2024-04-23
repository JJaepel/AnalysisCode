function Spines = calculateSpineHI(Spines)
% Caluculates homeogeneity index for spines on the same branch for a
% variety or radiuses

% Input:
% - Spines: all ROIs with their properties
% - dist: a 1xnumber of dist array
% 
% Output: 
% - Spines: a 1x num dist array stating homeogeneity index for different
% distances added as a new column -> the closer it is to 1, the more
% homeogeneous it its

%which dist do you want to look at?
dist = [2.5, 5, 10];

%add a new column to the Spines
Spines = arrayfun(@(x) setfield(x, 'HI', NaN(1,1)), Spines); 

%make sure we are only looking at the good spines
goodSpines = find([Spines.good] == 1);

for gs = 1:length(goodSpines)
    %select spines on the same dendrite
    s = goodSpines(gs);
    dendrNr = Spines(s).Dendrite; %get the dendrite Nr of the spine
    spineOnSameDend = find([Spines.Dendrite] == dendrNr); %find all spines on that dend
    goodSpinesOnSameDend = intersect(goodSpines, spineOnSameDend); %make sure they are responsive

    %get dist to all Spines
    allDist = abs(Spines(s).distToSoma - [Spines(goodSpinesOnSameDend).distToSoma]);
    
    for d = 1:length(dist)
        withinRadius = find(allDist < dist(d));
        preferences = zeros(length(withinRadius),1);
        for wrs = 1:length(withinRadius)
            preferences(wrs) = Spines(goodSpinesOnSameDend(withinRadius(wrs))).funcData.prefOri;
        end
        thresholdOrientations=exp(2*1i*preferences);
        HI(d) = 1-(abs(1/length(withinRadius) * sum(thresholdOrientations)));
        if HI(d) == 0
            HI(d) = NaN;
        end
    end
    Spines(s).HI = HI;
end

