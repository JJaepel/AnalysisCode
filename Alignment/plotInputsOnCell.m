function plotInputsOnCell(plotVeriInputs, plotAllInputs)

%switchboard
if nargin < 1
    plotVeriInputs = 1;
elseif nargin < 2
    plotAllInputs = 1;
    plotVeriInputs = 1;
end

%load the data
[File, DataPath] = uigetfile('*.mat', 'Please select file for loading cell data','Y:\Data_JJ_Ferret_project\JJ_FerretProject_Organized');
a = load([DataPath File]);

%% make an image with all dendrites
types = {a.STEDData.type};
typeInfo = cellfun(@(x) find(contains(x, 'Dendrite')),types,'UniformOutput',false);
notDendrite = cellfun(@isempty,typeInfo); %all other cells are empty
Dendrites = find(~notDendrite); 

newImage = zeros(size(a.CellImage,1), size(a.CellImage,2), length(Dendrites));
for d = 1:length(Dendrites)
    Left = a.STEDData(Dendrites(d)).Left; Up = a.STEDData(Dendrites(d)).Up;
    tempImage = a.STEDData(Dendrites(d)).Image;
    newImage(Up:Up+size(tempImage,1)-1,Left:Left+size(tempImage,2)-1,d) = tempImage;
end
newImageMean= mean(newImage,3);    
newImageMean = mat2gray(newImageMean);

figure
imagesc(newImageMean)
colormap('gray')
axis off
hold on

%% now plot all inputs on top
% 1.) Find spines
spineInfo = cellfun(@(x) find(contains(x, 'Spine')),types,'UniformOutput',false);
notSpines = cellfun(@isempty,spineInfo); %all other cells are empty
SpineIDs = find(~notSpines);

% 2.) Find inputs
inputCell = {a.STEDData.Input};
inputInfo = cellfun(@(x) find(contains(x, 'yes')),inputCell,'UniformOutput',false);
notInputs = cellfun(@isempty,inputInfo); %all other cells are empty
InputIDs = find(~notInputs);

% 3.) For those inputs, find the center of the image and plot it on  top

if plotAllInputs
    for i = 1:length(SpineIDs)
        spineNr = SpineIDs(i);
        PosLeft = a.STEDData(spineNr).Left + size(a.STEDData(spineNr).Image,2)/2;
        PosUp = a.STEDData(spineNr).Up + size(a.STEDData(spineNr).Image,1)/2;
        plot(PosLeft, PosUp, 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'Green');
    end
end

if plotVeriInputs
    for i = 1:length(InputIDs)
        spineNr = InputIDs(i);
        PosLeft = a.STEDData(spineNr).Left + size(a.STEDData(spineNr).Image,2)/2;
        PosUp = a.STEDData(spineNr).Up + size(a.STEDData(spineNr).Image,1)/2;
        plot(PosLeft, PosUp, 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'Red');
    end
end

inputData = struct;
for i=1:length(InputIDs)
    spineNr = InputIDs(i);
    inputData(i).Name = a.STEDData(spineNr).Name;
    inputData(i).Image = a.STEDData(spineNr).Image;
    inputData(i).PosLeft = a.STEDData(spineNr).Left;
    inputData(i).PosUp = a.STEDData(spineNr).Up;
end

%% save the image
todayDate = today('datetime');
F = getframe(gcf);
[X, ~] = frame2im(F);
if plotAllInputs
    if plotVeriInputs
        imwrite(X, [DataPath filesep 'Inputs_Checked_Verified' char(todayDate) '.tif'], 'tiff'); %all checked and verified inputs
    else
        imwrite(X, [DataPath filesep 'Inputs_Checked' char(todayDate) '.tif'], 'tiff'); %all checked inputs
    end
else
    if plotVeriInputs
        imwrite(X, [DataPath filesep 'Inputs_Verified' char(todayDate) '.tif'], 'tiff'); %all verified inputs
    else
        imwrite(X, [DataPath filesep 'Cell_Overview' char(todayDate) '.tif'], 'tiff'); %only cell overview
    end
end

