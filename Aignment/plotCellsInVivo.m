%switch variables
animal = 'F2688';
TwoPhontondir = 'Z:\Juliane\Data\2P_Data\';
bufferInUm = 30;

%NOTES: 1 unit of movement in the axis are 0.81 um! Try to first align D5
%and D6 using these measurements

%% 1. Find all dendrites
%read in the exp info for all 
filePath = 'F:\Organization\Animals\';
file = '2pExpByStimulusSpines.xlsx';
[~, xls_txt, xls_all]=xlsread([filePath file], 'driftingGrating');
exp_info = findExpInfo(xls_txt, xls_all);
animalInfo = cellfun(@(x) find(contains(x, animal)), exp_info.animal, 'UniformOutput', false); %which ones contain the animal
notAnimal = cellfun(@isempty,animalInfo); %all other cells are empty
indAnimal = find(~notAnimal); %all the dendrites
flags = find(exp_info.flag); %remove flagged experiments
indAnimal = indAnimal(indAnimal~=flags);

%% 2. For each dendrite recording, get information about size, position, the image and ROI positions
imagePos = struct;

for ind = 1:length(indAnimal)
    expID = char(exp_info.exp_id{indAnimal(ind)});
    animalName = char(exp_info.animal{indAnimal(ind)});
    disp(['Reading in ' char(expID)])
    
    imagePos(ind).expID= expID;
    imagePos(ind).type = char(exp_info.region{indAnimal(ind)});
    
    %read in header
    reader = ScanImageTiffReader([TwoPhontondir animalName filesep expID filesep expID '_00001_00001.tif']);
    meta=reader.metadata();

    %find position of sample in header
    startAxes = strfind(meta,'SI.hMotors.axesPosition = [');
    splitPos = strfind(meta,']');
    endAxes = splitPos(find(splitPos > startAxes, 1));
    axesPosLabel = meta(startAxes:endAxes);

    %split position into variables
    axesPosLabelCell = strsplit(axesPosLabel);
    imagePos(ind).xPos = str2double(axesPosLabelCell{3}(2:end));
    imagePos(ind).yPos = str2double(axesPosLabelCell{4});
    imagePos(ind).zPos = str2double(axesPosLabelCell{5}(1:end-1));

    %find size of image
    Zoom = regexp(meta, '(?<=scanZoomFactor = )\d+\.?\d*', 'match');
    sizeImage = 823.0500 / str2double(Zoom{1});
    imagePos(ind).sizeImage = sizeImage;
    imagePos(ind).zoom = Zoom;

    %add image of the region
    filePath = [TwoPhontondir animalName filesep expID filesep 'Registered\slice1\Projection\Projection.tif'];
    imagePos(ind).image = imread(filePath);
    
    %SCALE IMAGE IN TERMS OF MIN AND MAX TO MAKE THEM BETTER COMPARABLE BETWEEEN REGIONS!!!! 
    
    %add scaled image of 10 pixel per um
    scale = imagePos(ind).sizeImage*10/size(imagePos(ind).image,2);
    imagePos(ind).imageScaled = imresize(imagePos(ind).image,scale);
    
    %add ROI data
    ROIFile = load([TwoPhontondir animalName filesep expID filesep 'Registered\slice1\Projection\ROIs.mat']);
    imagePos(ind).ROIs = ROIFile.data.roi;
end

%% 3. Get parameters for plotting in x & y
%define region for plotting
minX = min([imagePos.xPos]);
maxX = max([imagePos.xPos]);
minY = min([imagePos.yPos]);
maxY = max([imagePos.yPos]);

startX = maxX; %let's set 0 to minX/minY - buffer for the extend of the imaging region
startY = maxY;

%define size of image in total - one z per image minus cell image
sizeY = size(imagePos(d).imageScaled,1);%(abs(maxX-minX)+bufferInUm*2)*10; %size in um -> lets have 10 pixels per um
sizeX = size(imagePos(d).imageScaled,1);%(abs(maxY-minY)+bufferInUm*2)*10;

%write X and Y position in relation to startX and startY
for d = 1:length(imagePos)
    imagePos(d).relXPos = abs(-imagePos(d).xPos + startX);
    imagePos(d).relYPos = abs(-imagePos(d).yPos + startY);
    imagePos(d).pixelXPos = abs(-imagePos(d).xPos + startX)*10;
    imagePos(d).pixelYPos = abs(-imagePos(d).yPos + startY)*10;
    
    imagePos(d).relXSomaPos = imagePos(1).xPos - imagePos(d).xPos;
    imagePos(d).relYSomaPos = imagePos(1).yPos - imagePos(d).yPos;
    imagePos(d).pixelXPos = abs(-imagePos(d).xPos + startX)*10;
    imagePos(d).pixelYPos = abs(-imagePos(d).yPos + startY)*10;
end

%sort the dendrites by their z position
[~, sortedDend] = sort([imagePos.zPos]);

%plot all dendrites in a different z plane
allDendrites = zeros(sizeX, sizeY, length(imagePos)-1);
z = 1;
for image = 1:length(imagePos)
    d = sortedDend(image);
    if length(imagePos(d).type) < 7  % 'cells' is less characters than 'dendrite'
        continue
    end
    
    startImageX = imagePos(d).pixelXPos-size(imagePos(d).imageScaled,1)/2;
    stopImageX = startImageX + size(imagePos(d).imageScaled,1);
    
    startImageY = imagePos(d).pixelYPos-size(imagePos(d).imageScaled,1)/2;
    stopImageY = startImageY + size(imagePos(d).imageScaled,2);
    
    allDendrites(startImageY:stopImageY-1,startImageX:stopImageX-1,z) = imagePos(d).imageScaled;
    z = z+1;
end

%make a projection of it 
a = mean(allDendrites,3);
figure
imagesc(a)

% figure
% for i = 1:6
%     subplot(2,3,i)
%     imagesc(allDendrites(:,:,i))
% end