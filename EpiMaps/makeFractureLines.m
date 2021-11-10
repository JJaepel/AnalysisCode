function fracture_filtered = makeFractureLines(ROI, activeFrames, corrTable, padValue, analysisParams)

ROILabeled = zeros(size(ROI));
ROILabeled(ROI) = 1:sum(ROI(:));

%adds zeros on th edges of the ROI
roiPad = padarray(ROILabeled,padValue,0,'both');
roiSidePad = zeros(size(roiPad,1),padValue);
roiPad = [roiSidePad roiPad roiSidePad];

%active frame with mask, side edges are NaN
mask = squeeze(activeFrames(:,:,1));
mask(mask==0) = NaN;

maskPad = padarray(mask,padValue,NaN,'both');
maskSidePad = zeros(size(maskPad,1),padValue);
maskPad = [maskSidePad maskPad maskSidePad];
maskPad(maskPad==0) = NaN;

%compute all correlation maps and store in cell array
corrMaps = {};
RoiVector = reshape(ROILabeled,[],1);

for i = 1:length(RoiVector)
    corrMapTemp =  recomputeImage(real(corrTable(i,:)),ROI);
    corrMapTemp(isnan(corrMapTemp))=0;
    corrMaps{i,1} = corrMapTemp;
end

%define storage for  the fracture map
fractureMap = zeros(size(maskPad,1),size(maskPad,2));
 
for x =(1+padValue):(length(maskPad) - padValue)
    tic
    for y = (1+padValue):(length(maskPad) - padValue) %length without pad zeros
        %if the index is NaN on the mask, set the fracture strenght to NaN
        if isnan(maskPad(x,y))
            fractureMap(x,y) = NaN;
        else
        
            %find the index for the central seed
            seedIndex = roiPad(x,y);
            
            %compute the correlation map for the initial seed
            corrInitial = corrMaps{seedIndex,1};
            
            %find all neighbors within a certain range
            indexRangeNeighbors = roiPad(x-padValue:x+padValue,y-padValue:y+padValue);
            indexRangeVec = reshape(indexRangeNeighbors,[],1);
            %index_range_m = index_range_mat .* index_range_m; % mask for values only in x,y directions
                      
            %for all neighbors, calculate the weighted correlation
            weightedCorrelation = zeros(size(indexRangeVec,1),1);
            
            
            for i = 1:length(indexRangeVec)
                neighborIndex = indexRangeVec(i);
                %find x and y coordinates
                if neighborIndex == 0 % if the neighbor is a padding value
                    weightedCorrelation(i) = NaN;
                elseif i == ceil(length(indexRangeVec)/2) %if it is the seed point, set to zero
                    weightedCorrelation(i) = NaN;           
                else
                    [row,col] = find(roiPad==neighborIndex); %find the position on the initial padded ROI
                    
                    if isnan(maskPad(row,col))
                        weightedCorrelation(i) = NaN; %if that is masked, set it to NaN
                    else
                        %otherwise compute corrleation map
                        corrFracture = corrMaps{neighborIndex,1};
                        %calc correaltion between both maps
                        corrAll = corr2(corrInitial,corrFracture);
                        %calc distance between initial and neighbour
                        distance = sqrt((x - row)^2+(y  - col)^2);
                        %weight correlation by distance
                        weightedCorrAll = (1-corrAll)/distance;
                        if isinf(weightedCorrAll)
                            weightedCorrelation(i) = NaN;
                        else
                            weightedCorrelation(i) = weightedCorrAll;
                        end
                    end
                end
                %average over all neighbours
                fractureMap(x,y) = nanmean(weightedCorrelation(:));
            end
        end
    end
    toc
    if mod(x,10) == 0
        disp(['Current row: ' num2str(x)])
    end
 end  

fracture_filtered = medfilt2(fractureMap);
fr = fractureMap; 
fr = fr - min(fr(:)) ; %normalisation
fr = fr / max(fr(:)) ;
    
figure; imagesc(fr)
colormap (jet)
colorbar;
saveas(gcf, fullfile(analysisParams.saveDirectory, ['FractureMap_Size' num2str(padValue) '.png']))