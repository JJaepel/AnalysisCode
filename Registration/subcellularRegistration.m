function subcellularRegistration(server, animal, name, level)

%This function registers the raw data and writes the new tif files in the 
%folder. Current version will be based on SI2016 format and BigTiff reader
%options to use downsample technique or nonrigid method or to use red
%channel to do registration.
%
%Input:
% -server: run it on RAID (prefered) or server?
% -animal: name of animal to define folder
% -name: exp name to define folder
% -level: is it a volume imaging?
%
%Output:
%- starts the registration, which will produce the registered tiffs, as
%well as the projection in the given folder


%% Step 0: Set some parameters
%type
regtype = 'downsample'; %Alternatives: 'dft'  
 
%downsample opts
downsampleRates = [1/8, 1/4, 1/2, 1];
maxMovement = 1/4;

%1/2 spatial resampling
imagSpatSamp = 0; 

%% Step 1: Data location and tfile folders
if server == 0
    drive = 'F:\';
else 
    drive = 'Z:\Juliane\';
end

TwoPhontondir = [drive 'Data\2P_Data\'];
tifDirectory= [TwoPhontondir animal filesep name];
baseTifDir = [tifDirectory '\Registered'];

%if it is volume images, you need to go through all the slices individually
if level
    numSlices = 5;
else
    numSlices = 1;
end

%go through one slice at the time
for slice = 1:numSlices
    
    %go to directory and make folders
    if ~exist(baseTifDir, 'dir')
        % make new file directory
        mkdir(baseTifDir); 
    end
    outputDir = [baseTifDir '\slice' num2str(slice) '\'];
    mkdir(outputDir)
    cd(tifDirectory);
    
    %% Step 2: Make the template
    disp(['reading in data from ' tifDirectory ' and grabbing templates']);
    template = [];
    
    %get all the tif files in the directory
    files  = dir('*.tif');
    
    tic;
    %go through the first files
    for fileNum = 1:4 %fix or do all? 
        
        %read in the data using ScanImageTiffReader
        imgStack = ScanImageTiffReader([cd,'/',files(fileNum).name]).data;
        imgStack = squeeze(imgStack);
        dat = squeeze(squeeze(sum(sum(imgStack,1),2)));
        
        %which image has the highest brightness?
        [~,id] = max(dat);
        
        %find the surrounding 40 images and make the mean image 
        if id < 15
            template(:,:,fileNum) = mean(imgStack(:,:,id:id+35),3);
        elseif id > size(imgStack,3)-25
            template(:,:,fileNum) = mean(imgStack(:,:,id-15:end),3);
        else
            template(:,:,fileNum) = mean(imgStack(:,:,id-14:id+25),3);
        end
    end
    
    %which template from the different files is the brightest
    [~,id] = max(sum(sum(template,1),2));
    template = squeeze(template(:,:,id)); %make this as your template
    
    %rotate & flip the template for saving
    templateRot = flip(template);
    templateRot = imrotate(templateRot,-90);
    
    %save it
    disp('Saving template');
    filename = [outputDir 'template.tif'];
    imwrite(uint16(templateRot),filename,'tif','writemode','overwrite');
    toc;
    
    %resize it if downsampling
    if imagSpatSamp
        template = imresize(template,.5);
    end

    %% Step 3: Read & Downsample files if selected
    % go through all files, one after each other
    for fileNum = 1:length(files)
        
        %read om the file
        imgStack = ScanImageTiffReader([cd,'/',files(fileNum).name]).data;
        imgStack = squeeze(imgStack);

        %2x spatial downsampling if selected
        if imagSpatSamp
            resamppx = size(imgStack,1)/2;
            for frnum = 1:size(imgStack,3)
                im = squeeze(imgStack(:,:,frnum));
                im = imresize(im,.5);
                im(im<0) = 0;
                imgStack(1:resamppx,1:resamppx,frnum) = im;
            end
            imgStack = imgStack(1:resamppx,1:resamppx,:);
            disp(['spat resamp and thresh done'])
        end

        %split channels if needed
        if level
            imgStack = imgStack(:,:,slice:numSlices:end);
        else
            imgStack = imgStack(:,:,1:end);
        end
        [height,width,depth] = size(imgStack);


        %% Step 4: Do the registration
        %Type A: Downsample: downsample in different rates (spatially) and find the
        %best registration for these areas
        if strcmp(regtype,'downsample')
            tic;
            regOffsetsX = zeros(depth,1);
            regOffsetsY = zeros(depth,2);
            
            %go through all downsample Rates to find the optimal Offsets
            for r=1:length(downsampleRates)
                sampRate = downsampleRates(r);
                if r>1
                    prevsampRate = downsampleRates(r-1);
                end
                disp(['registering images, iteration ' num2str(r)]);
                downHeight = height*sampRate;
                downWidth = width*sampRate;
                
                %resize the template acccordingly 
                templateImg = imresize(template, [downHeight,downWidth],'bilinear');
                for d=1:depth
                    regImg = imresize(imgStack(:,:,d), [downHeight,downWidth],'bilinear');
                    if r==1
                        %initial offset
                        minOffsetY = -round(maxMovement*downHeight/2);
                        maxOffsetY = round(maxMovement*downHeight/2);

                        minOffsetX = -round(maxMovement*downWidth/2);
                        maxOffsetX = round(maxMovement*downWidth/2);
                    else
                        %we are refining an earlier offset
                        minOffsetY = regOffsetsY(d)*sampRate - sampRate/prevsampRate/2;
                        maxOffsetY = regOffsetsY(d)*sampRate + sampRate/prevsampRate/2;

                        minOffsetX = regOffsetsX(d)*sampRate - sampRate/prevsampRate/2;
                        maxOffsetX = regOffsetsX(d)*sampRate + sampRate/prevsampRate/2;
                    end
                    bestCorrValue = -1;
                    bestCorrX = 0;
                    bestCorrY = 0;
                    for y=minOffsetY:maxOffsetY
                        for x=minOffsetX:maxOffsetX
                            %determine the offsets in X and Y for which the overlap
                            %between the images correlates best
                            subTemplateY1 = 1+max(y,0);
                            subTemplateY2 = downHeight+min(y,0);
                            subTemplateX1 = 1+max(x,0);
                            subTemplateX2 = downWidth+min(x,0);

                            subRegY1 = 1+max(-y,0);
                            subRegY2 = downHeight+min(-y,0);
                            subRegX1 = 1+max(-x,0);
                            subRegX2 = downHeight+min(-x,0);

                            subTemplateImg = templateImg(subTemplateY1:subTemplateY2,...
                                subTemplateX1:subTemplateX2);
                            subRegImg = regImg(subRegY1:subRegY2,subRegX1:subRegX2);

                            corrValue = corr2(subRegImg, subTemplateImg);
                            if corrValue > bestCorrValue
                                bestCorrX = x;
                                bestCorrY = y;
                                bestCorrValue = corrValue;
                            end
                        end
                    end
                    regOffsetsY(d) = bestCorrY*1/sampRate;
                    regOffsetsX(d) = bestCorrX*1/sampRate;
                end
            end

            disp('Registration offsets in (Y,X) format:');
            %now do the registration
            for d=1:depth
                img = imgStack(:,:,d);
                shiftedY1 = 1+max(regOffsetsY(d),0);
                shiftedY2 = height+min(regOffsetsY(d),0);
                shiftedX1 = 1+max(regOffsetsX(d),0);
                shiftedX2 = width+min(regOffsetsX(d),0);

                subRegY1 = 1+max(-regOffsetsY(d),0);
                subRegY2 = height+min(-regOffsetsY(d),0);
                subRegX1 = 1+max(-regOffsetsX(d),0);
                subRegX2 = width+min(-regOffsetsX(d),0);

                shiftedImg = ones(height,width)*double(min(min(img)));
                shiftedImg(shiftedY1:shiftedY2,shiftedX1:shiftedX2) = img(subRegY1:subRegY2,...
                    subRegX1:subRegX2);
                imgStack(:,:,d) = shiftedImg;
            end
            toc; 


        %Type B: Non-rigid
        elseif strcmp(regtype,'nonrigid')
            %set the options
            options_rigid = NoRMCorreSetParms('d1',size(imgStack,1),'d2',size(imgStack,2),...
                'grid_size',[128,128],'overlap_pre',64,'bin_width',50,'max_shift',25,'us_fac',50,'iter',1);
            options_nonrigid.use_parallel = 1;

            tic; 
            %correlate the stack to the template
            [M_final,~,~] = normcorre(imgStack,options_rigid,template); 
            imgStack = M_final;
            toc;
            
        %Type C: Dft, not implemented!!!
        elseif strcmp(regtype,'dft')
            disp('Lets try something new')
        end


        %% Step 5: Save the stacks as new files
        disp('Saving tiff stack');
        %get fo;e ma,e
        filename = [outputDir '\stack_c1_' sprintf('%02i',fileNum) '.tif'];
        for d = 1:depth
            %if it its the first frame, use overwrite, else just append it
            %with the other frames
            if d==1
                imwrite(uint16(imgStack(:,:,d)'),filename,'tif','writemode','overwrite');
            else
                imwrite(uint16(imgStack(:,:,d)'),filename,'tif','writemode','append');
            end
        end

        disp(['Done, Saved in  ', outputDir]);
    end
end


