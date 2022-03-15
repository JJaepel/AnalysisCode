function seedBasedCorrDelay(Movie, index, i)

%interactive version of seedBasedCorr

% Movie=imresize(Movie,0.5);
sz=size(Movie);
Movie=reshape(Movie,sz(1)*sz(2),sz(3));
Amean = mean(Movie,2); %avg at each pixel location in the image over time
CorMovie=reshape(Amean,[sz(1),sz(2)]);

figure
while 1
    imagesc(CorMovie);     % display image
    colormap(jet)
    caxis([0 1])
    [xi, yi, but] = ginput(1);      % get a point
    if ~isequal(but, 1)             % stop if not button 1
        break
    end
    xi=ceil(xi); yi=ceil(yi);
    ind=sub2ind([sz(1) sz(2)], yi, xi); 

    pix=squeeze(Movie(ind,1:end-i));
    corM=corr(Movie(index,1+i:end)',pix');
    CorMovie=reshape(corM,sz(1),sz(2));
    CorMovie(yi-1:yi+1,xi-1:xi+1)=1;
end

% iList=index(5:500:length(index));
% for i=1:length(iList)
%     pix=squeeze(Movie(iList(i),:));
%     corM=corr(Movie(index,:)',pix');
%     frameCor=zeros(sz(1)*sz(2),1); frameCor(index)=corM;
%     frameCor=reshape(frameCor,sz(1),sz(2));
%     
%     
%     [r,c]=ind2sub(sz,iList(i));
%     frameCor(r-3:r+3,c-3:c+3)=1;
%     CorMovie(:,:,i)=frameCor;
%     i
% end
% beep

% B=CorMovie;
% B(B<0)=0;
% B(B>.25)=0.25;
% B=rescale(B);
% B=floor(B*255);


% v = VideoWriter('CorMovie_filtered10.avi', 'Indexed AVI');
% v.FrameRate=2;
% v.Colormap=jet(256);
% open(v);
% for k = 1:size(B,3) 
%    writeVideo(v,(B(:,:,k)));
% end
% 
% close(v);