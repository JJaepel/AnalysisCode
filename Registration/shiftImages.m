function [registeredImage] = shiftImages(sourceImage,row_shift, col_shift, shiftInDiscreteValues)

if(nargin < 4), shiftInDiscreteValues = true; end

% Compute registered version of sourceImage
if(shiftInDiscreteValues)
    % shift image discretely by whole pixel values
    row_shift = round(row_shift);
    col_shift = round(col_shift);
    registeredImage = circshift(sourceImage,[row_shift col_shift]);
    
    % make rows and cols zero where we shifted
    if(row_shift>0), registeredImage(1:row_shift,:,:) = 0; elseif(row_shift<0), registeredImage((end+row_shift+1):end,:,:) = 0; end
    if(col_shift>0), registeredImage(:,1:row_shift,:) = 0; elseif(col_shift<0), registeredImage(:,(end+row_shift+1):end,:) = 0; end         
else
    for currentChannel = 1:size(sourceImage,3) 
        fftedSrcImg = fft2(padImage(sourceImage(:,:,currentChannel)));

        if(usfac > 0)
            [nr,nc]=size(fftedSrcImg);
            Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
            Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
            [Nc,Nr] = meshgrid(Nc,Nr);
            registeredChannel = fftedSrcImg.*exp(1i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
            registeredImagePadded(:,:,currentChannel) = registeredChannel*exp(1i*diffphase);
        elseif(usfac == 0)
            registeredImagePadded(:,:,currentChannel) = fftedSrcImg*exp(1i*diffphase);
        end
        registeredImagePadded(:,:,currentChannel) = abs(ifft2(registeredImagePadded(:,:,currentChannel)));
        registeredImage(:,:,currentChannel) = removePadding(registeredImagePadded(:,:,currentChannel),srcImagePadding);
    end
end

return

function [paddedImage,paddingInfo] = padImage(img)
% pad images and put on square array of size 2^N (zero padding)
paddingInfo.imgSize     = size(img);
paddingInfo.ext         = round(2.^(ceil(log(max(paddingInfo.imgSize))/log(2)))); %(This will be the size of the filters/arrays)
paddingInfo.YLocation   = round((paddingInfo.ext-paddingInfo.imgSize(1))/2)+[1:paddingInfo.imgSize(1)]; 
paddingInfo.XLocation   = round((paddingInfo.ext-paddingInfo.imgSize(2))/2)+[1:paddingInfo.imgSize(2)];
paddedImage = zeros(paddingInfo.ext);  
paddedImage(paddingInfo.YLocation,paddingInfo.XLocation) = img;
return
function img = removePadding(paddedImage,paddingInfo)
% removing padding
img = paddedImage(paddingInfo.YLocation,paddingInfo.XLocation);
return