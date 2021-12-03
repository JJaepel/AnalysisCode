function analysis = makeMasks(data, analysis, analysisParams, verbose)

%use PCA to make masks, 1st dimension is window, 2nd is areas
[~, mixedfilters, ~] = SVDsimple(data.rawF,1);
data.PCAsRaw = mixedfilters;

%make maks for whole window
figure
imagesc(data.PCAsRaw(:,:,1))
saveas(gcf, [analysisParams.saveDirectory, 'Mask1stPCA_Raw.png'])
close gcf

mask = data.PCAsRaw(:,:,1);
if mask(1,1) > 0
    mask(mask > 0) = 0;
    mask(mask < 0) = 1;
else
    mask(mask > 0) = 1;
    mask(mask < 0) = 0;
end
analysis.mask = logical(mask);

%remove BV
analysis.maskBV = removeBVfromMask(analysis, analysisParams.changeThreshold);

%make mask for subregions
tempMask = data.PCAsRaw(:,:,2);
tempMask(~analysis.maskBV)=0;

posMask = tempMask;
posMask(posMask > 0) = 1;
posMask(posMask < 0) = 0;

negMask = tempMask;
negMask(negMask > 0) = 0;
negMask(negMask < 0) = 1;

figure
imagesc(tempMask)
colorbar
if verbose
    posMask19=isequal(input('Are positive correlations in area 19? (Y/N): ','s'),'Y');
    if posMask19
        analysis.maskV1 = negMask;
        analysis.maskA19 = posMask;
    else
        analysis.maskV1 = posMask;
        analysis.maskA19 = negMask;
    end
else
    disp('Positive correlations = A19 - Please change if needed!')
    analysis.maskV1 = negMask;
    analysis.maskA19 = posMask;
end
saveas(gcf, [analysisParams.saveDirectory, 'Mask2ndPCA_Raw.png'])
close gcf