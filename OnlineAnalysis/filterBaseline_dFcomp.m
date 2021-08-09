%remove slow baselin in raw F traces
%(1) cut off large events 
%(2) 99 pt medfilt for low-pass trace
%(3) calc initial F value (median)
%(4) subtract and add back
function [raw_new]=filterBaseline_dFcomp(raw,pt)
if nargin<2
    pt= 99;
end
raw_new = raw;

% % % % raw_new(raw_new > std(raw)+median(raw)) = median(raw);

% pt-order (window size),one dimensional median filter
raw_new = medfilt1(raw_new,pt); 

% raw_new = cat(1,median(raw).*ones(100,1),raw_new);
% raw_new = cat(1,median(raw).*ones(100,1),raw_new);
% raw_new = prctfilt1(raw_new,90);
% raw_new = raw_new(101:end-100);

raw_new = (raw - raw_new)./raw_new;



