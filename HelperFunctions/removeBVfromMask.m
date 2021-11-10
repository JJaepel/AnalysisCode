function MaskBV = removeBVfromMask(analysis, changeThreshold)

Image = analysis.baseImg;
Image=Image./imgaussfilt(Image,size(Image,1)./60);
Image(~analysis.mask)=0;

NoFiltImage=analysis.baseImg;
NoFiltImage(~analysis.mask)=0;

Perc=10;
Thresh=percentile2D(Image(analysis.mask),Perc,1);

figure
subplot(1,2,1)
imagesc(NoFiltImage)
axis off
axis equal
subplot(1,2,2)
imagesc(Image)
axis off
axis equal
caxis([Thresh Thresh+.001])
if changeThreshold
    Done = 0;
else
    Done = 1;
end

while not(Done)
Good=isequal(input(['Percentile is ' num2str(Perc) '. Looks good? (Y/N): '],'s'),'Y');
if Good
    Done=1;
else
    Perc=input('Input new percentile: ');
    Thresh=percentile2D(Image(analysis.mask),Perc,1);
    caxis([Thresh Thresh+.001])
end
end

MaskBV=(analysis.mask&(Image>Thresh));
