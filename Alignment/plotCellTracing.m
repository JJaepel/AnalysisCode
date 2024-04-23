confRecon = 'Z:\Saachi\Tracing\F2688LSL05Tracing\ConfocalTraces';
confReconFiles = dir([confRecon '\*.swc']); %tracing data
imgFile = dir([confRecon '\*.jpg']); %projection
ROIFiles = dir([confRecon filesep 'Spine Tracings' filesep '\*.zip']); %spine or ROI data
imgSizeConf = [4060, 4060]; %number of pixels
umSize = [302.1508, 302.1508]; %size in microns 

plotImage = 1;

%% load data
%load image
template = imread([confRecon filesep imgFile.name]);
template = imadjust(template);

%load tracing Data
xTrace = [];
yTrace = [];
for f = 1:length(confReconFiles)
    file = [confRecon filesep confReconFiles(f).name];
    b = loadfilelist(file);
    a = b(4:end);
    c = cellfun(@str2num, a, 'UniformOutput', false);
    traceInfo = cell2mat(c');
    xTrace = [xTrace; traceInfo(:,3)];
    yTrace = [yTrace; traceInfo(:,4)];
end

%load ROI data
cvsROIs = [];
for r = 1:length(ROIFiles)
    file = [confRecon filesep 'Spine Tracings' filesep ROIFiles(r).name];
    [ROI] = ReadImageJROI(file);
    cvsROIs = [cvsROIs ROI];
end

%% process data
%calc umfactor to change from micrometer to pixelposition
umfactor = imgSizeConf(1)/umSize(1);

%separate tracing data into dendrites by finding jumps
dffXTrace = diff(xTrace); dffYTrace = diff(yTrace);
jumpsMins = find(dffXTrace < -1.5); jumpsMinsY = find(dffYTrace < -1.5);
jumpMax = find(dffXTrace > 1.5); jumpMaxY = find(dffYTrace > 1.5);
jumps = [jumpsMins; jumpMax; jumpsMinsY; jumpMaxY];
jumps = sort(jumps);
jumps = unique(jumps);

x = 1;
for j = 1:length(jumps)
    dendrite{j}.xTrace = xTrace(x:jumps(j))*umfactor;
    dendrite{j}.yTrace = yTrace(x:jumps(j))*umfactor;
    x = jumps(j)+1;
end

%process ROIs
for i = 1:length(cvsROIs)
    tempR = poly2mask(cvsROIs{i}.mnCoordinates(:,1), cvsROIs{i}.mnCoordinates(:,2),imgSizeConf(1),imgSizeConf(2));
    %CHECK FOR DOUBLE REGIONPROPS!!! USE LARGEST ONE
    STATS = regionprops(tempR,'Centroid', 'Area', 'PixelList', 'ConvexHull' );
    [~, idx] = max([STATS.Area]);
    
    ROIs(i).perimeter = cvsROIs{i}.mnCoordinates+1;
    ROIs(i).body = STATS.PixelList;
    if length(ROIs(i).perimeter) > 100
        ROIs(i).type = 'soma';
    else
        ROIs(i).type = 'spine';
    end
    ROIs(i).size = length(ROIs(i).body(:,2));
    ROIs(i).xPos = STATS(idx).Centroid(:,1)+1;
    ROIs(i).yPos = STATS(idx).Centroid(:,2)+1;
end

%% plot data
figure
%imgData
if plotImage
    imagesc(template)
    axis image
    colormap('gray')
else
    set(gca, 'YDir','reverse')
end
hold on

%traces
for dend = 1:length(dendrite)
    plot(dendrite{dend}.xTrace ,dendrite{dend}.yTrace, 'LineWidth',1.5, 'color', 'red');
    hold on
end

%spines
numROI = 0;
for i = 1:length(ROIs)
    type = ROIs(i).type;
    switch type
        case 'spine'
            xpos= ROIs(i).xPos;
            ypos= ROIs(i).yPos;
            plot(xpos,ypos,'ok','MarkerSize',2,'MarkerFaceColor', 'red');
            hold on
            numROI = numROI+1;
        case 'soma'
            plot(ROIs(i).body(:,1), ROIs(i).body(:,2),'color', 'red');
            hold on
    end
end
disp(num2str(numROI))
axis off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(confRecon, ['CellTracing' num2str(plotImage) '.png']))

