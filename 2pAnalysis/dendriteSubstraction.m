function ce = dendriteSubstraction(ce, info)

% This function removes the dendritic component from the spine ROIs
%
% Inputs:
% - ce: structure containing in each row one ROI with anatomical
% information, raw trace, dff trace and stimulus chopped trace = cyc
% - info: structure containing all the information about this experiment
%
% Ouput:
% - ce: updated with cycRes, containing the spineOnly signal trace, and
% corr as a measurement of how much the spine is still related to the
% dendrite
% -ROIChopped.mat: ROI file containing the updated ce

%% Step 1: Go through all ROIs
for cc = 1:length(ce)
    %decided whether it is a spine
    if ce(cc).spine
        %find dendrite
        dendPoint = cc;
        while ~ce(dendPoint).dendrite %go forward until you hit the dendrite
            dendPoint = dendPoint+1;
        end
        
        %get the spine signal, multiply it with 0.9, remove inf
        Spdff = ce(cc).dff;
        Spdff = Spdff(1:round(length(Spdff)*.9));
        Spdff(isinf(Spdff)) = 0;
        
        %calculate noise by finding trace below median
        Spdff_sub = Spdff(Spdff < nanmedian(Spdff)+abs(min(Spdff)));
        %noiseM = nanmedian(Spdff_sub); % should be near 0
        noiseSD = nanstd(Spdff_sub);
        
        %fit the slope between dendrite and spine
        slope = robustfit(ce(dendPoint).cyc(:),ce(cc).cyc(:));

        %dendritic scalar applied (from robust fit) to cyc and subtraction
        ce(cc).slope = slope(2);
        ce(cc).cycRes = ce(cc).cyc - slope(2).*ce(dendPoint).cyc;
        ce(cc).cycRes(ce(cc).cycRes < 0) = 0;
        
        %slight variation to get the rawres
        rSp =  ce(cc).dff;
        rDn =  ce(dendPoint).dff;
        rSp = rSp - slope(2).*rDn;
        rSp(rSp < -noiseSD) = -noiseSD;
        rSp(isinf(rSp)) = 0;
        rDn(isinf(rDn)) = 0;
        ce(cc).rawRes = rSp;
        rSp(rSp <= 0) = nan;
        rDn(rDn <= 0) = nan;
        
        %how much are spine and dendrite still correlated?
        [r,~] = corrcoef(rSp,rDn,'rows','pairwise');
        ce(cc).corr = r(2);

    else
        %if it is a dendrite, just set these empty
        ce(cc).cycRes = [];
        ce(cc).slope = [];
        ce(cc).corr = -1;
    end
end

%% Step 2: Save it
save([info.saveDir 'ROIsChopped.mat'], 'ce', '-mat') 

