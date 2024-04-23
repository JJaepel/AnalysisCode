function correlation_domains(corrTable, ROI)

%this function makes mask from the correlation domains. It creates masks with positive
%and negative correlations and sepatares them into single domains


seed_filtered = (recomputeImage(real(analysis.corrTable(26667,:)),analysis.ROI)); 
seed_filtered(isnan(seed_filtered))=0;

seed_mask = seed_filtered;
seed_mask(seed_mask > 0.1) = 1;
seed_mask(seed_mask < 0.1) = 0;
seed_mask = logical(seed_mask);
figure; imagesc(seed_mask)

domains_positive = bwpropfilt(seed_mask,'perimeter',20);
figure; imagesc(domains_positive)

seed_mask_negative = seed_filtered;

seed_mask_negative(seed_mask_negative > -0.05) = 0;
seed_mask_negative(seed_mask_negative < -0.05) = 1;
seed_mask_negative = logical(seed_mask_negative);
figure; imagesc(seed_mask_negative)

domains_negative = bwpropfilt(seed_mask_negative,'perimeter',20);
figure; imagesc(domains_negative)

domains = bwpropfilt(seed_mask,'perimeter',30);
figure; imagesc(domains)


[B,L] = bwboundaries(domains,'noholes'); %L - matrix with labeled domains



single_positive = L;
single_positive(single_positive > 7) = 0;
single_positive(single_positive < 7) = 0;


[B,L] = bwboundaries(domains_negative,'noholes'); %L - matrix with labeled domains



single_negative = L;
single_negative(single_negative > 5) = 0;
single_negative(single_negative < 5) = 0;



FramesDomains = zeros(size(data.activeFrames));

for i = 1:size(data.activeFrames,3)
    frame_mask = data.activeFrames(:,:,i);
    frame_mask(isnan(frame_mask))=0;
    mean_f = mean(mean(frame_mask));
    std_f = std(std(frame_mask));
    frame_mask(frame_mask<(mean_f + 2.5* std_f)) = 0;
    frame_mask(frame_mask>(mean_f + 2.5* std_f)) = 1;
    frame_mask = logical(frame_mask);
    domains_frame = bwpropfilt(frame_mask,'perimeter',40);
     figure;imagesc(domains_frame)

    [B,L] = bwboundaries(domains,'holes');
FramesDomains(:,:,i) = L; 
end
    figure;
    imshow(label2rgb(L, @gray, [.5 .5 .5]))
    hold on
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1)
    end
end

end