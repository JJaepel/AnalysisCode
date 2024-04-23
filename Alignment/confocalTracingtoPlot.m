function confocalTracingtoPlot(cellInfo, cell)

confDir = ['Z:\Juliane\InputAnalysis\' char(cellInfo.animal{cell}) filesep char(cellInfo.cellName{cell}) '\B - Confocal\' char(cellInfo.slice{cell}) filesep];
saveDir = [confDir 'Dendrite reconstruction'];
if ~exist(saveDir, 'dir')
    % make new file directory
    mkdir(saveDir); 
end

confReconFiles = dir([confDir '\Dendrites\*.swc']); %tracing data
imgFile = dir([confDir '\*.jpg']); %projection
confFile = dir([confDir '\*.tif']);
ROIFiles = dir([confDir filesep 'Spines' filesep '\*.zip']); %spine or ROI data
%number of pixels
plotImage = 1;

%% load image
%load image
template = imread([confDir filesep imgFile.name]);
template = imadjust(template);
imgSizeConf = size(template); 

%get the micrometer to pixel factor from the confocal image
a = imfinfo([confDir filesep confFile.name]);
umfactor = a(1).XResolution; %x & y resolution
ImageDesc = a(1).ImageDescription;

%% load tracing Data -> for each dendrite after each other
ROIcounter = 1;
denNr = 1;
ROIs = struct;
Branches = struct;
allDendrites = [];
dendNumApi = cellInfo.apicalDend{cell};
dendriteNrs = cellInfo.dendriteNr{cell};
if length(dendriteNrs) > 1
    temp = strsplit(dendriteNrs, ','); %split them by the comma
    dendriteNrs = cellfun(@ (x) str2double(x), temp); %change to numbers
end
if dendriteNrs(1) > 1
    denNr = 101;
end
for e = 1:length(dendriteNrs)
    d = dendriteNrs(e);
    
    %%DENDRITES
    %first load the dendrite tracings for that dendrite and find how many
    %segments it is
    dendrite = struct;
    dendriteCounter = 1;
    segmentNr = 1;
    
    %apical or basal?
    apInd = find(dendNumApi == d);
    if isempty(apInd)
        type = 'basal';
    else
        type = 'apical';
    end
    
    filesDend = cellfun(@(x) contains(x, ['Dendrite' sprintf('%02d',d)]),{confReconFiles.name}); %which are the files for that dendrite?
    segment = struct;
    for f = 1:length(filesDend)
        if filesDend(f) == 1
            file = [confDir filesep 'Dendrites' filesep confReconFiles(f).name];
            b = loadfilelist(file);
            a = b(4:end);
            c = cellfun(@str2num, a, 'UniformOutput', false);
            traceInfo = cell2mat(c');
            Trace = traceInfo(:,3:5);
            
            %find if there are branch points in that file
            dffXTrace = diff(Trace(:,1)); dffYTrace = diff(Trace(:,2));
            jumpsMins = find(dffXTrace < -1.5); jumpsMinsY = find(dffYTrace < -1.5);
            jumpMax = find(dffXTrace > 1.5); jumpMaxY = find(dffYTrace > 1.5);
            jumps = [jumpsMins; jumpMax; jumpsMinsY; jumpMaxY];
            jumps = sort(unique(jumps));            
            %save all segments
            if ~isempty(jumps)
                x = 1;
                for j = 1:length(jumps)
                    segment(segmentNr).Trace = Trace(x:jumps(j),:);
                    x = jumps(j)+1;
                    segmentNr = segmentNr+1;
                end
                %add the last segment
                segment(segmentNr).Trace = Trace(x:end,:);
                segmentNr = segmentNr+1;
            else
                segment(segmentNr).Trace = Trace;
                segmentNr = segmentNr+1;
            end
            
        end
    end
    startPoint =  segment(1).Trace(1,:);
    c = 1;
    %now find connection points between segments
    branchPointC = 1;
    branchPoint = zeros(1,4); %coordinates and dist to soma
    branchLocation = zeros(1,4); %nearestSegment, Pos on Segment, Origin Segment
    for s=1:segmentNr-1
        %whats the minimal distance between the start of this segment and all other
        %segments?
        distSeg = zeros(segmentNr-1,1); Pos =distSeg;
        for ss=1:segmentNr-1
            t=ones(length(segment(ss).Trace),1);
            dist = sqrt((segment(s).Trace(1,1)*t - segment(ss).Trace(:,1)).^2 + (segment(s).Trace(1,2)*t - segment(ss).Trace(:,2)).^2);
            [distSeg(ss),Pos(ss)] = min(dist);
        end
        Pos(s) = nan; distSeg(s) = inf;
        [minDist,nearSeg] = min(distSeg);
        if Pos(nearSeg) > 1
            branchPoint(branchPointC,1:3)= segment(nearSeg).Trace(Pos(nearSeg),:);
            branchPoint(branchPointC,4) = pdist([branchPoint(branchPointC,1:3);startPoint]);
            branchLocation(branchPointC,:) = [s, 1, nearSeg,Pos(nearSeg)];
            branchPointC = branchPointC+1;
        end
    end
    %get all the unique branchpoints and the segments
    diffBranchPoints = abs(diff(branchPoint(:,1)))+abs(diff(branchPoint(:,2))); %find branchpoints that are very similar in position (within 1 um in each direction)
    doublet = find(diffBranchPoints < 2);
    branchPoint(doublet,:) = [];
    branchLocation(doublet,:) = [];
    [~, order] = sort(branchPoint(:,4)); %order them by distance from the starting point
    branchPoints = zeros(length(order),8);
    for i = 1:length(order)
        branchPoints(i,1) = branchPoint(order(i),1); %xCoord
        branchPoints(i,2) = branchPoint(order(i),2); %yCoord
        branchPoints(i,3) = branchPoint(order(i),3); %zCoord
        branchPoints(i,4) = branchLocation(order(i),1); %first Segment
        branchPoints(i,5) = branchLocation(order(i),2); %first Segment Pos
        branchPoints(i,6)= branchLocation(order(i),3); %second Segment
        branchPoints(i,7) = branchLocation(order(i),4); %second Segment Pos
        branchPoints(i,8) = branchPoint(order(i),4); %dist from soma
    end
    
    %let's assume that the first segment starts from the cell to the first
    %branch point of that segment
    SPseg1 = find(branchPoints(:,4) == 1);
    SPseg2 = find(branchPoints(:,6) == 1);
    SPseg = sort([branchPoints(SPseg1,5); branchPoints(SPseg2,7)]);
    
    dendrite(dendriteCounter).Branch= d; %what branch are we starting at?
    dendrite(dendriteCounter).dendNr= denNr;
    dendrite(dendriteCounter).startPoint = 0; %starts at soma = 0
    dendrite(dendriteCounter).endPoint = 1; %first branchpoint = 1
    if segmentNr > 2
        dendrite(dendriteCounter).coord = segment(1).Trace(1:SPseg(1),:);
        dendrite(dendriteCounter).pixelCoord = segment(1).Trace(1:SPseg(1),1:2)*umfactor;
    else
        dendrite(dendriteCounter).coord = segment(1).Trace;
        dendrite(dendriteCounter).pixelCoord = segment(1).Trace(:,1:2)*umfactor;
    end
    %generally, the distance between coord is 0.0744 -> length of trace is
    %number of coord -1 *0.0744
    dendrite(dendriteCounter).length = (size(dendrite(dendriteCounter).coord,1)-1) *0.0744;
    dendrite(dendriteCounter).denOnBranch = dendriteCounter;
    switch type
        case 'apical' 
            dendrite(dendriteCounter).BranchOrder = 2;
            dendrite(dendriteCounter).distToSoma = 0; %%Cave: change that to reflect the difference in z
        case 'basal'
            dendrite(dendriteCounter).BranchOrder = 1;
            dendrite(dendriteCounter).distToSoma = 0;
    end
    dendrite(dendriteCounter).type = type;
    dendriteCounter = dendriteCounter +1;
    denNr = denNr+1;
    
    %now lets go through all branchpoints
    if segmentNr > 2
        for bP = 1:size(branchPoints,1)
            % 1) Lets start with seg1 which starts at this position
            dendrite(dendriteCounter).Branch= d;
            dendrite(dendriteCounter).dendNr= denNr;
            dendrite(dendriteCounter).startPoint = bP;
            SPseg1 = branchPoints(bP,4); %what segment is the branch on? should start at 1
            %let's see if that branch is having a branchpoint later on
            bEnd = find(branchPoints(:,6) == SPseg1);
            if isempty(bEnd) %if not, then the whole dendrite can be put into this
                dendrite(dendriteCounter).endPoint = NaN;
                dendrite(dendriteCounter).coord = segment(SPseg1).Trace;
                dendrite(dendriteCounter).pixelCoord = segment(SPseg1).Trace(:,1:2)*umfactor;
            else %otherwise, look for the next branching point off this segment
                endPos = sort(branchPoints(bEnd,7));
                dendrite(dendriteCounter).coord = segment(SPseg1).Trace(1:endPos(1),:);
                dendrite(dendriteCounter).pixelCoord = segment(SPseg1).Trace(1:endPos(1),1:2)*umfactor;
                %what is the nearest branchpoint?
                t=ones(size(branchPoints,1),1);
                dist = sqrt((dendrite(dendriteCounter).coord(end,1)*t - branchPoints(:,1)).^2 + (dendrite(dendriteCounter).coord(end,2)*t - branchPoints(:,2)).^2);
                [~,dendrite(dendriteCounter).endPoint] = min(dist);

            end
            dendrite(dendriteCounter).length = (size(dendrite(dendriteCounter).coord,1)-1) *0.0744;
            %now let's look how far away it is from the soma by searching for
            %the previous dendrite
            previousDend = find([dendrite.endPoint] == dendrite(dendriteCounter).startPoint);
            dendrite(dendriteCounter).distToSoma = dendrite(previousDend).length + dendrite(previousDend).distToSoma;
            dendrite(dendriteCounter).denOnBranch = dendriteCounter;
            dendrite(dendriteCounter).BranchOrder = dendrite(previousDend).BranchOrder+1;
            dendrite(dendriteCounter).type = type;
            dendriteCounter = dendriteCounter +1;
            denNr = denNr+1;

            % 2) Now look at the seg2 from which it branches off
            dendrite(dendriteCounter).Branch= d;
            dendrite(dendriteCounter).dendNr= denNr;
            dendrite(dendriteCounter).startPoint = bP;
            SPseg2 = branchPoints(bP,6); %what segment is the branch on? 
            SPseg2Pos = branchPoints(bP,7);% what position does it start on?
            %let's see if that branch is having further branchpoints later on;
            bSeg2End = find(branchPoints(:,6) == SPseg2);
            bSeg2EndPos = sort(branchPoints(bSeg2End,7));
            bEnd = find(bSeg2EndPos > SPseg2Pos);
            if isempty(bEnd) %if not, then the whole dendrite can be put into this
                dendrite(dendriteCounter).endPoint = NaN;
                dendrite(dendriteCounter).coord = segment(SPseg2).Trace(SPseg2Pos:end,:);
                dendrite(dendriteCounter).pixelCoord = segment(SPseg2).Trace(SPseg2Pos:end,1:2)*umfactor;
            else %otherwise, use this until you get to the next branching point off this segment
                dendrite(dendriteCounter).coord = segment(SPseg2).Trace(SPseg2Pos:bSeg2EndPos(bEnd(1)),:);
                dendrite(dendriteCounter).pixelCoord = segment(SPseg2).Trace(SPseg2Pos:bSeg2EndPos(bEnd(1)),1:2)*umfactor;
                %what is the nearest branchpoint?
                t=ones(size(branchPoints,1),1);
                dist = sqrt((dendrite(dendriteCounter).coord(end,1)*t - branchPoints(:,1)).^2 + (dendrite(dendriteCounter).coord(end,2)*t - branchPoints(:,2)).^2);
                [~,dendrite(dendriteCounter).endPoint] = min(dist);
            end
            dendrite(dendriteCounter).length = (size(dendrite(dendriteCounter).coord,1)-1) *0.0744;
            %now let's look how far away it is from the soma by searching for
            %the previous dendrite
            previousDend = find([dendrite.endPoint] == dendrite(dendriteCounter).startPoint);
            dendrite(dendriteCounter).distToSoma = dendrite(previousDend).length + dendrite(previousDend).distToSoma;        
            dendrite(dendriteCounter).denOnBranch = dendriteCounter;
            dendrite(dendriteCounter).BranchOrder = dendrite(previousDend).BranchOrder+1;
            dendrite(dendriteCounter).type = type;
            dendriteCounter = dendriteCounter +1;
            denNr = denNr+1;
        end
    end
    
    %%ROIs
    %Load the spine data for this dendrite
    filesSpines = cellfun(@(x) contains(x, ['Dendrite' sprintf('%02d',d)]),{ROIFiles.name}); %which are the files for that dendrite?
    file = [confDir filesep 'Spines' filesep ROIFiles(filesSpines).name];
    [cvsROIs] = ReadImageJROI(file);
    
    %process ROIs
    for r = 1:size(cvsROIs,2)
        tempR = poly2mask(cvsROIs{r}.mnCoordinates(:,1), cvsROIs{r}.mnCoordinates(:,2),imgSizeConf(1),imgSizeConf(2));
        %CHECK FOR DOUBLE REGIONPROPS!!! USE LARGEST ONE
        STATS = regionprops(tempR,'Centroid', 'Area', 'PixelList', 'ConvexHull' );
        [~, idx] = max([STATS.Area]);
        
        if length(cvsROIs{r}.mnCoordinates) > 200 %then it is the soma
            continue
        end
        
        ROIs(ROIcounter).Nr = sprintf('%04d',ROIcounter);
        ROIs(ROIcounter).perimeter = cvsROIs{r}.mnCoordinates+1;
        try
            ROIs(ROIcounter).body = STATS.PixelList;
            ROIs(ROIcounter).size = length(ROIs(ROIcounter).body(:,2));
            ROIs(ROIcounter).xPos = STATS(idx).Centroid(:,1)+1;
            ROIs(ROIcounter).yPos = STATS(idx).Centroid(:,2)+1;
        catch
            ROIs(ROIcounter).body = cvsROIs{r}.mnCoordinates(1,:);
            ROIs(ROIcounter).size = 1;
            ROIs(ROIcounter).xPos = cvsROIs{r}.mnCoordinates(1,1);
            ROIs(ROIcounter).yPos = cvsROIs{r}.mnCoordinates(1,2);
        end
        ROIs(ROIcounter).slice = str2double(cvsROIs{r}.strName(1:4));
        
        %find the right z position by looking for the closest dendrite 
        %whats the minimal distance between the spine and dendrites 2D
        minDenDist = zeros(1,length(dendrite), 1);
        denPos = minDenDist;
        for den=1:size(dendrite,2)
            SpineX = ROIs(ROIcounter).xPos;
            SpineY = ROIs(ROIcounter).yPos;
            t = ones(length(dendrite(den).coord),1);
            denDist = sqrt((SpineX*t - dendrite(den).pixelCoord(:,1)).^2 + (SpineY*t - dendrite(den).pixelCoord(:,2)).^2);
            [minDenDist(den),denPos(den)] = min(denDist);
        end
        [~,nearDend]= min(minDenDist);
        ROIs(ROIcounter).zPos = dendrite(nearDend).coord(denPos(nearDend),3);
        ROIs(ROIcounter).Branch = d;
        ROIs(ROIcounter).ROINrOnBranch = r;
        ROIs(ROIcounter).Dendrite = dendrite(nearDend).dendNr;
        ROIs(ROIcounter).denOnBranch = dendrite(nearDend).denOnBranch;
        ROIs(ROIcounter).distToSoma = dendrite(nearDend).distToSoma+denPos(nearDend)*0.0744;
        ROIs(ROIcounter).BranchOrder = dendrite(nearDend).BranchOrder; 
        ROIs(ROIcounter).type = type;
        ROIcounter= ROIcounter+1;
    end
    
    denColors = cbrewer('div', 'Spectral', length(dendrite));
    denColors(denColors>1) = 1;
    dendOnBranch = [dendrite.dendNr];
    
    %lets plot to see how this is doing for different branches
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
        plot(dendrite(dend).pixelCoord(:,1) ,dendrite(dend).pixelCoord(:,2), 'LineWidth',1.5, 'color', denColors(dend,:));
        hold on
        %spines
        SpinesOnDen = find([ROIs.Dendrite] == dendOnBranch(dend));
        if ~isempty(SpinesOnDen)
            for i = 1:length(SpinesOnDen)
                xpos= ROIs(SpinesOnDen(i)).xPos;
                ypos= ROIs(SpinesOnDen(i)).yPos;

                plot(xpos,ypos,'ok','MarkerSize',2,'MarkerFaceColor', denColors(dend,:));
                hold on
            end
        end
    end
    
    %save information about branch
    Branches(e).type = type;
    Branches(e).dendrites = dendrite;
    ROIsOnBranch = find([ROIs.Branch] == d);
    Branches(e).ROIs = ROIs(ROIsOnBranch);
    
    %coordinates for loading only that branch data in Matching
    pixelsAll = {dendrite.pixelCoord};
    Branches(e).xMin = floor(min(cell2mat(cellfun(@(x) min(x(:,1)), pixelsAll, 'UniformOutput',false)))); 
    Branches(e).xMax = ceil(max(cell2mat(cellfun(@(x) max(x(:,1)), pixelsAll, 'UniformOutput',false))));
    Branches(e).yMin = floor(min(cell2mat(cellfun(@(x) min(x(:,2)), pixelsAll, 'UniformOutput',false)))); 
    Branches(e).yMax = ceil(max(cell2mat(cellfun(@(x) max(x(:,2)), pixelsAll, 'UniformOutput',false)))); 
    
    slices = {dendrite.coord};
    slicesROIs = cell2mat({Branches(e).ROIs.slice});
       
    Branches(e).zMin = min(floor(min(cell2mat(cellfun(@(x) min(x(:,3)), slices, 'UniformOutput',false)))), min(slicesROIs)); 
    Branches(e).zMax = max(ceil(max(cell2mat(cellfun(@(x) max(x(:,3)), slices, 'UniformOutput',false)))), max(slicesROIs));

    allDendrites = [allDendrites, dendrite];
    clear dendrite
    
    axis off
    set(gcf, 'color', 'w');
    %save tracing for each dendrite
    saveas(gcf, fullfile(saveDir, ['CellTracing_Dendrite_' sprintf('%02d',d) '.png']))
end

%% add the soma
somaFile = dir([confDir filesep 'Spines' filesep '*Soma*' ]);
file = [confDir filesep 'Spines' filesep somaFile(1).name];
[cvsROI] = ReadImageJROI(file);
tempR = poly2mask(cvsROI.mnCoordinates(:,1), cvsROI.mnCoordinates(:,2),imgSizeConf(1),imgSizeConf(2));
STATS = regionprops(tempR,'Centroid', 'Area', 'PixelList', 'ConvexHull' );
Soma.xPos = STATS.Centroid(:,1)+1;
Soma.yPos = STATS.Centroid(:,2)+1;
Soma.perimeter = cvsROI.mnCoordinates;

%% plot all tracing
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

denColors = cbrewer('div', 'Spectral', length(allDendrites));
denColors(denColors>1) = 1;
for allD = 1:length(allDendrites)
    plot(allDendrites(allD).pixelCoord(:,1) ,allDendrites(allD).pixelCoord(:,2), 'LineWidth',1.5, 'color', denColors(allD,:));
    hold on
    
    %spines
    SpinesOnDen = find([ROIs.Dendrite] == allD);
    if ~isempty(SpinesOnDen)
        for i = 1:length(SpinesOnDen)
            xpos= ROIs(SpinesOnDen(i)).xPos;
            ypos= ROIs(SpinesOnDen(i)).yPos;

            plot(xpos,ypos,'ok','MarkerSize',2,'MarkerFaceColor', denColors(allD,:));
            hold on
        end
    end
end

axis off
set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDir, 'CellTracing_AllDendrites.png'))

%% save data
save([confDir filesep 'cellReconstruction.mat'], 'allDendrites','ROIs', 'Branches','Soma','-mat') 


