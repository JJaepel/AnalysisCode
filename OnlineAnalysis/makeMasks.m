function analysis = makeMasks(rawF, analysis, saveDirectory, changeThreshold,verbose)

if nargin < 5
    verbose = 1;
end
if nargin < 4
    changeThreshold = 0;
end

%use PCA to make masks, 1st dimension is window, 2nd is areas
[~, mixedfilters, ~] = SVDsimple(rawF,1);
PCAsRaw = mixedfilters;

%make maks for whole window
figure
imagesc(PCAsRaw(:,:,1))
saveas(gcf, [saveDirectory, 'Mask1stPCA_Raw.png'])
close gcf

mask = PCAsRaw(:,:,1);
if mask(1,1) > 0
    mask(mask > 0) = 0;
    mask(mask < 0) = 1;
else
    mask(mask > 0) = 1;
    mask(mask < 0) = 0;
end
analysis.mask = logical(mask);

%remove BV
analysis.maskBV = removeBVfromMask(analysis, changeThreshold);

%make mask for subregions
tempMask = PCAsRaw(:,:,2);
tempMask(~analysis.maskBV)=0;

posMask = tempMask;
posMask(posMask > 0) = 1;
posMask(posMask < 0) = 0;

negMask = tempMask;
negMask(negMask > 0) = 0;
negMask(negMask < 0) = 1;

figure
imagesc(tempMask)
saveas(gcf, [saveDirectory, 'Mask2ndPCA_Raw.png'])
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
close gcf