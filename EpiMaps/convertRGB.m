function imrgb = convertRGB(mtx)
    % Red-White-Blue polarization map for ON-OFF
    stimnum= size(mtx, 3);
    cmap = lbmap(stimnum, 'RedBlue');   
    mtx = (mtx - min(mtx(:)))/ (max(mtx(:))-min(mtx(:))); 
    imrgb=zeros(size(mtx,1), size(mtx,2), 3);
   
    for i = 1:2
        imrgb(:,:,1) = imrgb(:,:,1)+ cmap(i,1) * mtx(:,:,i);
        imrgb(:,:,2) = imrgb(:,:,2)+cmap(i,2) * mtx(:,:,i);
        imrgb(:,:,3) = imrgb(:,:,3)+cmap(i,3) * mtx(:,:,i);
    end
end