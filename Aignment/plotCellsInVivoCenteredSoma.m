%switch variables
animal = 'F2688';
TwoPhontondir = 'Z:\Juliane\Data\2P_Data\';
anaDir = 'F:\Data\ImageAnalysis\';
bufferInUm = 30;
posToUmFactor = 0.82;
responses = 1;

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
for f = 1:length(flags)
    indAnimal = indAnimal(indAnimal~=flags(f));
end

%% 2. For each dendrite recording, get information about size, position, the image and ROI positions
imagePos = struct;

for ind = 1:length(indAnimal)
    expID = char(exp_info.exp_id{indAnimal(ind)});
    animalName = char(exp_info.animal{indAnimal(ind)});
    disp(['Reading in ' char(TwoPhontondir) char(animalName) '/' char(expID) '/' char(expID) '_00001_00001.tif'])
    
    imagePos(ind).expID= expID;
    imagePos(ind).type = char(exp_info.region{indAnimal(ind)});
    
    %read in header
    reader = ScanImageTiffReader([TwoPhontondir animalName filesep expID filesep expID '_00001_00001.tif']);
    meta=reader.metadata();

    %find position of sample in header
    startAxes = strfind(meta,'SI.hMotors.axesPosition = [');
    %startSample = strfind(meta,'SI.hMotors.samplePosition = [');
    splitPos = strfind(meta,']');
    endAxes = splitPos(find(splitPos > startAxes, 1));
    %endSample = splitPos(find(splitPos > startSample, 1));
    axesPosLabel = meta(startAxes:endAxes);
    %SamplePosLabel = meta(startSample:endSample);

    %split position into variables
    axesPosLabelCell = strsplit(axesPosLabel);
    imagePos(ind).xPos = str2double(axesPosLabelCell{3}(2:end));
    imagePos(ind).yPos = str2double(axesPosLabelCell{4});
    imagePos(ind).zPos = str2double(axesPosLabelCell{5}(1:end-1));

%     SamplePosLabelCell = strsplit(SamplePosLabel);
%     samplePos(1) = str2double(SamplePosLabelCell{3}(2:end));
%     samplePos(2) = str2double(SamplePosLabelCell{4});
%     samplePos(3) = str2double(SamplePosLabelCell{5}(1:end-1));
    %imagePos(ind).samplePos = samplePos;

    %find size of image
    Zoom = regexp(meta, '(?<=scanZoomFactor = )\d+\.?\d*', 'match');
    sizeImage = 823.0500 / str2double(Zoom{1});
    imagePos(ind).sizeImage = sizeImage;
    imagePos(ind).zoom = Zoom;

    %add image of the region
    filePath = [TwoPhontondir animalName filesep expID filesep 'Registered\slice1\Projection\Projection.tif'];
    rawImage = imread(filePath);
    imagePos(ind).image = mat2gray(rawImage);
    
    %SCALE IMAGE IN TERMS OF MIN AND MAX TO MAKE THEM BETTER COMPARABLE BETWEEEN REGIONS!!!! 
    
    %add scaled image of 10 pixel per um
    scale = imagePos(ind).sizeImage*10/size(imagePos(ind).image,2);
    imagePos(ind).imageScaled = imresize(imagePos(ind).image,scale);
    
    %add ROI data
    ROIFile = load([TwoPhontondir animalName filesep expID filesep 'Registered\slice1\Projection\ROIs.mat']);
    imagePos(ind).ROIs = ROIFile.data.roi;
    
    %add functional information
    if responses
        anaData = load([anaDir animalName filesep expID filesep 'AnaData.mat']);
        anaROI = anaData.analysis.rawRes.roi;
        for r = 1:length(imagePos(ind).ROIs)
            imagePos(ind).ROIs(r).prefOri = anaROI(r).preferredOrientation;
            imagePos(ind).ROIs(r).OSIFit = anaROI(r).OSIFit;
            imagePos(ind).ROIs(r).isResp = anaROI(r).isResponseSignificant;
        end
    end
end
disp('Lets take a break here')

%% 3. Get parameters for plotting in x & y
%define region for plotting
minX = min([imagePos.xPos]);
maxX = max([imagePos.xPos]);
minY = min([imagePos.yPos]);
maxY = max([imagePos.yPos]);

startX = maxX+bufferInUm*2; %let's set 0 to minX/minY - buffer for the extend of the imaging region
startY = maxY+bufferInUm*2;

%define size of image in total - one z per image minus cell image
sizeY = (abs(maxX-minX)+bufferInUm*2)*10; %size in um -> lets have 10 pixels per um
sizeX = (abs(maxY-minY)+bufferInUm*2)*10;

%write X and Y position in relation to startX and startY
for d = 1:length(imagePos)
    imagePos(d).relXPos = abs(-imagePos(d).xPos + startX);
    imagePos(d).relYPos = abs(-imagePos(d).yPos + startY);
    imagePos(d).pixelXPos = abs(-imagePos(d).xPos + startX)*10*posToUmFactor;
    imagePos(d).pixelYPos = abs(-imagePos(d).yPos + startY)*10*posToUmFactor;
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
    
    startImageX = ceil(imagePos(d).pixelXPos)-size(imagePos(d).imageScaled,1)/2;
    imagePos(d).startImageX = startImageX;
    stopImageX = ceil(startImageX + size(imagePos(d).imageScaled,1));
    
    startImageY = ceil(imagePos(d).pixelYPos)-size(imagePos(d).imageScaled,1)/2;
    imagePos(d).startImageY = startImageY;
    stopImageY = ceil(startImageY + size(imagePos(d).imageScaled,2));
    
    allDendrites(startImageY:stopImageY-1,startImageX:stopImageX-1,z) = imagePos(d).imageScaled;
    z = z+1;
end

%make a projection of it 
a = mean(allDendrites,3);
figure
imagesc(a)
colormap('gray')
axis image
hold on

%plot ROI positions on top of it
disp('Now plot ROIs')
LUT = hsv(180);
for image = 1:length(imagePos)
    if length(imagePos(image).type) < 7  % 'cells' is less characters than 'dendrite'
        continue
    end
    for r = 1:length(imagePos(image).ROIs)
        xpos = imagePos(image).ROIs(r).xPos +  imagePos(image).startImageX;
        ypos = imagePos(image).ROIs(r).yPos +  imagePos(image).startImageY;
        
        if responses
            if imagePos(image).ROIs(r).isResp == 1
               plot(xpos,ypos,'ok','MarkerSize',5','MarkerFaceColor',LUT(1+floor(imagePos(image).ROIs(r).prefOri),:));
                hold on
            end
        else
            plot(xpos,ypos,'ok','MarkerSize',10','MarkerFaceColor','red');
            hold on
        end
    end
end