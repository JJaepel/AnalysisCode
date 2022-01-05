function [outputArg1,outputArg2] = onset_V1_A19(HighF,maskBV, maskV1,maskA19,dff)
%This function takes hight-pass spatial filter from dF/F signal for the whole ROI.
% It applies mask to remove BV and separates masks from PC2 to separate V1 and A19.
%Then the sum over space is calculated. Next, it takes derivatives for 
%detection of the onset of events

%dff = data.dff;
% 0. take high-pass filter from raw df/f frames
 analysis.ROIActive =true( [size(dff,1),size(dff,2)]);
 [LowHighData, HighF ] = LowHighNormalize(double(dff), analysis.ROIActive,maskBV);

% maskBV = analysis.maskBV; maskV1 = analysis.maskV1; HighF = data.HighF; maskA19 = analysis.maskA19; 
% 1. apply blood vessels mask and V1 PC2 mask to the high-pass df/f filter
v_one = HighF .* maskBV .* maskV1;
v_once_size = size(v_one);
        x  = v_once_size(1);
        y  = v_once_size(2);
        t  = v_once_size(3);
v_one_elements = nnz(v_one(:,:,1));    
v_one_sum = reshape(v_one, [x*y t]);
v_one_sum = sum(v_one_sum,1); 
v_one_sum(isnan(v_one_sum))=0; % sum of events over space for V1 events
v_one_sum = zscore(v_one_sum);
% 2. make X vector, translating frames into sec
frame_time = importdata('C:\Users\rulah\Desktop\frametrigger.txt');

frame_time_one = mean(diff(frame_time));
x_vec = zeros(size(v_one_sum,2));
all_time = frame_time(end) - frame_time(1);
    
    for i = 2 : size(v_one_sum,2)
        x_vec(1) =  frame_time_one;
        x_vec(i) =  x_vec(i-1) + frame_time_one;
    end
    x_vec = x_vec(:,1)';
%3. take derivatives and threshhold them by std, remove neighbouring values
%V1

dydx = gradient(v_one_sum) ./ gradient(x_vec);
derivatives_one = zscore(dydx);
onset_der_one = gradient(derivatives_one);
peaks_der_one = find(onset_der_one>std(onset_der_one)); % indices for frames with onset

for i = 1:length(peaks_der_one)-1
 if peaks_der_one(i+1) - peaks_der_one(i) < 3
    peaks_der_one(i+1) = NaN;
end 
end
peaks_der_one = peaks_der_one(~isnan(peaks_der_one));

onset_one = NaN(1,length(v_one_sum));

for i = 1: length(peaks_der_one)
onset_one(peaks_der_one(i)) = v_one_sum(peaks_der_one(i));
end
onset_one(onset_one>0.3) = NaN; % to remove onset values on peaks

%4. same for A19
a_nin = HighF .* maskBV .* maskA19;
a_nin_size = size(a_nin);
        x  = a_nin_size(1);
        y  = a_nin_size(2);
        t  = a_nin_size(3);
a_nin_elements = nnz(a_nin(:,:,1));    
a_nin_sum = reshape(a_nin, [x*y t]);
a_nin_sum = sum(a_nin_sum,1);
a_nin_sum(isnan(a_nin_sum))=0; % sum of events over space for V1 events
a_nin_sum = zscore(a_nin_sum);

dydx_n = gradient(a_nin_sum) ./ gradient(x_vec);
derivatives_nin = zscore(dydx_n);
onset_der_nin = gradient(derivatives_nin);
peaks_der_nin = find(onset_der_nin>std(onset_der_nin)); % indices for frames with onset

for i = 1:length(peaks_der_nin)-1
 if peaks_der_nin(i+1) - peaks_der_nin(i) < 3
    peaks_der_nin(i+1) = NaN;
end 
end
peaks_der_nin = peaks_der_nin(~isnan(peaks_der_nin));

onset_nin = NaN(1,length(a_nin_sum));

for i = 1: length(peaks_der_nin)
onset_nin(peaks_der_nin(i)) = a_nin_sum(peaks_der_nin(i));
end
onset_nin(onset_nin>0.3) = NaN; % to remove onset values on peaks

%5. 

% cross corrrelation for 2 vectors with peaks derivatives onsets
[c,lags] = xcorr(peaks_der_one(1,2:end), peaks_der_nin(:));
[~,lag_value] = max(c);
lag_value = lags(1,lag_value); %shows the value of mismatch in the number of frames
%figure; stem(lags,c)

bool_der = length(peaks_der_one) > length(peaks_der_nin); %defines length of the next 'for loop' as length of bigger vector
if bool_der == 0
    len_derivative = length(peaks_der_nin);
else 
    len_derivative = length(peaks_der_one);
end

peak_first = []; %takes first element as first peak 
if peaks_der_one(1,1) - peaks_der_nin(1,1) > 0 || peaks_der_one(1,1) - peaks_der_nin(1,1) == 0
    peak_first = 2; %peaks_der_nin(1,1);
else
    peak_first = 1; %peaks_der_one(1,1);
end

%substract each element in the second vector from the first vector and
%compare, if it's > then 5, then it's a single spike 

peak_order = [];
for i = 1:5
for ii = 1:len_derivative
    peak_order(i) = find(abs(peaks_der_one(1,1) - peaks_der_nin(1,1))< 5);
end
end

if isempty(peak_order) == 1



%%version 2
%1. separate events into individual and simultaneous

%1 = individual peak,V1. 2 = individual peak,A19. 3 = V1 first. 4 = A19
%first. 5 = simultaneous onset at the same frame
peak_order_one = zeros(1, length(peaks_der_one));

for loop_vone = 1:length(peaks_der_one)
    for loop_vnin = 1:length(peaks_der_nin)
        simultaneous = find(abs(peaks_der_one(1,loop_vone) - peaks_der_nin(1,loop_vnin)) < 5);
        if isempty(simultaneous) == 0
            if peaks_der_one(1,loop_vone) - peaks_der_nin(1,loop_vnin) <0
                peak_order_one(1,loop_vone) = 3; % V1 first   
            elseif peaks_der_one(1,loop_vone) - peaks_der_nin(1,loop_vnin) >0
                peak_order_one(1,loop_vone) = 4; % A19 first
            else
                peak_order_one(1,loop_vone) = 5; 
            end
            break
        else
            continue
        end
    end

    if isempty(simultaneous)
       peak_order_one(1,loop_vone) = 1;
    end
end

%for A19
%1 = individual peak,V1. 2 = individual peak,A19. 3 = V1 first. 4 = A19
%first. 5 = simultaneous onset at the same frame
peak_order_nineteen = zeros(1, length(peaks_der_nin));

for loop_vnin = 1:length(peaks_der_nin)
    for loop_vone = 1:length(peaks_der_one)
        simultaneous = find(abs(peaks_der_nin(1,loop_vnin) - peaks_der_one(1,loop_vone)) < 5);
        if isempty(simultaneous) == 0
            if peaks_der_nin(1,loop_vnin) - peaks_der_one(1,loop_vone)  < 0
                peak_order_nineteen(1,loop_vnin) = 4; % A19 first   
            elseif peaks_der_nin(1,loop_vnin) - peaks_der_one(1,loop_vone) >0
                peak_order_nineteen(1,loop_vnin) = 3; % V1 first
            else
                peak_order_nineteen(1,loop_vnin) = 5; 
            end
            break
        else
            continue
        end
    end

    if isempty(simultaneous)
       peak_order_nineteen(1,loop_vnin) = 2;
    end
end

% Q: 1. how many events in V1/A1
%    2. how many independant events in V1/A19
%    3. events originate first in V1/A19

data_onset(1,1) = length(peaks_der_one); % number of events in V1
data_onset(1,2) = length(peaks_der_nin); % number of events in A19

%1 = individual peak,V1. 2 = individual peak,A19. 3 = V1 first. 4 = A19
%first. 5 = simultaneous onset at the same frame
count_ons_one = accumarray(peak_order_one',1); 
count_ons_nin = accumarray(peak_order_nineteen',1);

%number of peaks in the high-pass filter
 [peaks,locations] = findpeaks(a_nin_sum); 
   
    peaks_std = std(peaks);
    peask_mean = mean(peaks);
    threshold = peask_mean + .05 * peaks_std;
    
    peak_thresholded = find(peaks>threshold);
    peak_thresholded = locations(peak_thresholded,1);

   


%plotting
figure; plot(x_vec,v_one_sum,'r')
hold on 
plot(x_vec,a_nin_sum, 'b')
hold on
plot(x_vec,onset_nin, 'b*')
hold on
plot(x_vec, onset_one,'r+')



end