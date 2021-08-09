
    
function resizedTifStack = read_Tiffs2019(filePath,imgScaling,updateFrequency,useWaitBar,dataType,codeVersion)
% tifStack = read_Tiffs(filePath,updateFrequency)
%
% A fast method to read a tif stack: calls the tiff library directly.
% Used the same method described in this web-link: 
% http://www.matlabtips.com/how-to-load-tiff-stacks-fast-really-fast/
% 
% filePath - fileDirectory and name of image stack
% imgScaling - Scales size of images by specified scaling factor (default is no scaling).
%              0.5 drops image size by half, while 2x doubles image size
% updateFrequency - Based on the number of images, how often in percentage
% of files read to inform user of current progress (0-100%)
% useWaitBar - Displays the progress of reading in the current stack.
% dataType - data type of the imaging data (default is 'uint16', but can be 'single' or 'float')
% codeVersion - Can specify whether to use the 'tifflib' function to read
%               images, the newer 'Tiff' MATLAB object to read images (default method), 
%               or the 'ScanImage' tiff reader (fastest).
%
% Written by David Whitney (10/3/2013)
% Updated on 12/20/2018
% David.Whitney@mpfi.org
% Max Planck Florida Institude

if(nargin<2), imgScaling      = 0.5;      end
if(nargin<3), updateFrequency = 5;        end
if(nargin<4), useWaitBar      = false;    end
if(nargin<5), dataType        = 'uint16'; end
if(nargin<6), codeVersion     = 'ScanImage';   end
t0=tic();
disp(['Reading Image Stack - ' filePath]);

% Read TIF Header and Setup Information
if(~strcmp(codeVersion,'ScanImage'))
    useImfinfo = true;
    if(useImfinfo)
        InfoImage=imfinfo(filePath);
        xImage=InfoImage(1).Width;
        yImage=InfoImage(1).Height;
        NumberOfImages=length(InfoImage);
    else
        hTif = Tiff(filePath);
        firstImage = hTif.read();
        xImage=size(firstImage,2);
        yImage=size(firstImage,1);
        NumberOfImages=1; 
        while(~hTif.lastDirectory)
            hTif.nextDirectory();
            NumberOfImages=NumberOfImages+1; 
        end
    end

    % Initialize MATLAB array to contain tif stack
    originalTifStack = zeros(yImage,xImage,NumberOfImages,dataType);

    disp(['Finished Reading Image Header - ' num2str(toc(t0)) ' seconds Elapsed']);
else
    NumberOfImages = 1;
end
tHeader = tic();

% setup wait bar
if(useWaitBar)
    h = waitbar(0,'Opening Tif image...', 'Name', 'Open TIF Image', 'Pointer', 'watch');
    currentPosition = get(h,'Position');
    offset = 100;
    set(h,'Position',[currentPosition(1)-offset/2 currentPosition(2) currentPosition(3)+offset currentPosition(4)]);
    currentPosition = get(get(h,'Children'),'Position');
    set(get(get(h,'Children')),'Position',[currentPosition(1)-offset/2 currentPosition(2) currentPosition(3)+offset currentPosition(4)]);
else
    updateFrequency = round((updateFrequency/100)*NumberOfImages);
end

% Get imaging data
switch codeVersion
    case 'ScanImage'
        tiffReader = ScanImageTiffReader(filePath);
        originalTifStack = tiffReader.data();
        originalTifStack = permute(originalTifStack,[2 1 3]);
        tiffReader.close();
    case 'tifflib'
        % uses the tifflib function to read images fast
        currentImage     = zeros(yImage,xImage);
        FileID = tifflib('open',filePath,'r');
        rps    = tifflib('getField',FileID,Tiff.TagID.RowsPerStrip);
        for i=1:(NumberOfImages)
            if(useWaitBar)
                waitbar(i/NumberOfImages,h, ['Image ' num2str(i) ' of ' num2str(NumberOfImages) ' - ' num2str(toc) 's Elapsed - ' num2str((NumberOfImages-i)*toc/i) 's Left']);
            else
                if(mod(i+1,updateFrequency)==0)
                    disp([num2str(round(100*i/NumberOfImages)) '% Done Reading Image Stack - ' num2str(toc) ' seconds Elapsed']);
                end
            end

            tifflib('setDirectory',FileID,i-1);
            warning('OFF','MATLAB:imagesci:tiffmexutils:libtiffWarning');
            % Go through each strip of data.
            rps = min(rps,yImage);
            for r = 1:rps:yImage
              row_inds = r:min(yImage,r+rps-1);
              stripNum = tifflib('computeStrip',FileID,r)-1;
              currentImage(row_inds,:) = tifflib('readEncodedStrip',FileID,stripNum);
            end
            originalTifStack(:,:,i) = currentImage;
        end
        tifflib('close',FileID);
        warning('ON','MATLAB:imagesci:tiffmexutils:libtiffWarning')
    case 'Tiff'
        % Setup TIF object and Read-In Basic Information
        hTif = Tiff(filePath);
          
        warning('OFF','MATLAB:imagesci:tiffmexutils:libtiffWarning');
        for i=1:(NumberOfImages)
            if(useWaitBar)
                waitbar(i/NumberOfImages,h, ['Image ' num2str(i) ' of ' num2str(NumberOfImages) ' - ' num2str(toc) 's Elapsed - ' num2str((NumberOfImages-i)*toc/i) 's Left']);
            else
                if(mod(i+1,updateFrequency)==0)
                    disp([num2str(round(100*i/NumberOfImages)) '% Done Reading Image Stack - ' num2str(toc(t0)) ' seconds Elapsed']);
                end
            end

            originalTifStack(:,:,i) = hTif.read();
            if(i == NumberOfImages)
                hTif.close();
                warning('ON','MATLAB:imagesci:tiffmexutils:libtiffWarning')
            else
                hTif.nextDirectory();                
            end
        end
end
if(strcmp(codeVersion,'ScanImage'))
    yImage         = size(originalTifStack,1);
    xImage         = size(originalTifStack,2);
    NumberOfImages = size(originalTifStack,3);
end
disp(['Finished Reading Image Stack - ' num2str(toc(tHeader)) ' seconds Elapsed']);

% If necessary rescale images
if(imgScaling ~= 1 && imgScaling>0)
    scaledX = round(xImage*imgScaling);
    scaledY = round(yImage*imgScaling);
    resizedTifStack  = zeros(scaledY,scaledX,NumberOfImages,dataType);
    for(i=1:NumberOfImages)
        resizedTifStack(:,:,i) = imresize(originalTifStack(:,:,i),[scaledY scaledX]); % Scales image size
    end
else
    resizedTifStack = originalTifStack;
end      
if(useWaitBar),close(h); end

