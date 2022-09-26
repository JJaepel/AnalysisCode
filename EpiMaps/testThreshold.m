function threshold = testThreshold(image)

h=figure;
set(h, 'units','normalized','outerposition',[0 0 1 1]);

subplot(1,2,1)
contour(image);
set(gca,'ydir','reverse')
axis off
axis equal
colorbar

threshold = 1.5;
temp = zeros(size(image,1), size(image,2));
temp(image<threshold) = 0;
temp(image>threshold) = 1;

subplot(1,2,2)
imagesc(temp)
axis off
axis equal
set(gcf, 'color', 'w');

Done = 0;
while not(Done)
    Good=isequal(input(['Threshold is ' num2str(threshold) '. Looks good? (Y/N): '],'s'),'Y');
    if Good
        Done=1;
    else
        threshold=input('Input new threshold: ');
        temp = zeros(size(image,1), size(image,2));
        temp(image<threshold) = 0;
        temp(image>threshold) = 1;
        imagesc(temp)
        axis off
        axis equal
    end
end
close gcf


