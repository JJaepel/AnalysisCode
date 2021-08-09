function tifStack = readTiffsScaled(filePath,imgScaling,updateFrequency,useWaitBar)
    warning off;
    if(nargin<2), imgScaling = 0.5;    end
    if(nargin<3), updateFrequency = 5; end
    if(nargin<4), useWaitBar = false;  end
    tic;

    disp(['Reading Image Stack - ' filePath]);

    % Read TIF Header and Setup Information
    InfoImage=imfinfo(filePath);
        xImage=InfoImage(1).Width;
        yImage=InfoImage(1).Height;
        NumberOfImages=length(InfoImage);
        disp(['Finished Reading Image Header - ' num2str(toc) ' seconds Elapsed']);

    % use wait bar
    if(useWaitBar)
        h = waitbar(0,'Opening Tif image...', 'Name', 'Open TIF Image', 'Pointer', 'watch');
    %     currentPosition = get(h,'Position');
    %     offset = 100;
    %     set(h,'Position',[currentPosition(1)-offset/2 currentPosition(2) currentPosition(3)+offset currentPosition(4)]);
    %     currentPosition = get(get(h,'Children'),'Position');
    %     set(get(get(h,'Children')),'Position',[currentPosition(1)-offset/2 currentPosition(2) currentPosition(3)+offset currentPosition(4)]);
    else
        updateFrequency = round((updateFrequency/100)*NumberOfImages);
    end

    % Initialize MATLAB array to contain tif stack
    scaledX = round(xImage*imgScaling);
    scaledY = round(yImage*imgScaling);
    tifStack     = zeros(scaledY,scaledX,NumberOfImages,'uint16');

    codeVersion = 'alternativeMethod'; % both methods seem pretty similar in performance, but original method is a bit faster
    switch codeVersion
        case 'originalMethod'
            % uses the tifflib function to read images fast
            FileID = tifflib('open',filePath,'r');
            rps    = tifflib('getField',FileID,Tiff.TagID.RowsPerStrip);
            rps    = min(rps,yImage);
            for i=1:(NumberOfImages)
                % display read progress
                if(useWaitBar)
                    waitbar(i/NumberOfImages,h, ['Image ' num2str(i) ' of ' num2str(NumberOfImages) ' - ' num2str(toc) 's Elapsed - ' num2str((NumberOfImages-i)*toc/i) 's Left']);
                else
                    if(mod(i+1,updateFrequency)==0)
                        disp([num2str(round(100*i/NumberOfImages)) '% Done Reading Image Stack - ' num2str(toc) ' seconds Elapsed']);
                    end
                end

                % turn off warnings
                warning('OFF','MATLAB:imagesci:tiffmexutils:libtiffWarning');

                % Go through each strip of data.
                currentImage = zeros(yImage,xImage);
                tifflib('setDirectory',FileID,i-1);
                for r = 1:rps:yImage
                  row_inds = r:min(yImage,r+rps-1);
                  stripNum = tifflib('computeStrip',FileID,r)-1;
                  currentImage(row_inds,:) = tifflib('readEncodedStrip',FileID,stripNum);
                end

                % Rescale data
                if(imgScaling ~= 1 && imgScaling>0)
                    tifStack(:,:,i) = imresize(currentImage,[scaledY scaledX]); % Scales image size
                else
                    tifStack(:,:,i) = currentImage;
                end
            end
            tifflib('close',FileID);
            disp(['Finished Reading Image Stack - ' num2str(toc) ' seconds Elapsed']);
            warning('ON','MATLAB:imagesci:tiffmexutils:libtiffWarning')
        case 'alternativeMethod'
            % Setup TIF object and Read-In Basic Information
            hTif = Tiff(filePath);

            warning('OFF','MATLAB:imagesci:tiffmexutils:libtiffWarning');
            for i=1:NumberOfImages
                if(useWaitBar)
                    waitbar(i/NumberOfImages,h, ['Image ' num2str(i) ' of ' num2str(NumberOfImages) ' - ' num2str(toc) 's Elapsed - ' num2str((NumberOfImages-i)*toc/i) 's Left']);
                else
                    if(mod(i+1,updateFrequency)==0)
                        disp([num2str(round(100*i/NumberOfImages)) '% Done Reading Image Stack - ' num2str(toc) ' seconds Elapsed']);
                    end
                end

                if(imgScaling ~= 1 && imgScaling>0)
                    tifStack(:,:,i) = imresize(hTif.read(),[scaledY scaledX]); % Scales image size
                else
                    tifStack(:,:,i) = hTif.read();
                end
                if(i == NumberOfImages)
                    hTif.close();
                    warning('ON','MATLAB:imagesci:tiffmexutils:libtiffWarning')
                    disp(['Finished Reading Image Stack - ' num2str(toc) ' seconds Elapsed']);
                else
                    hTif.nextDirectory();                
                end
            end
    end

    if(useWaitBar),close(h); warning on; end
end