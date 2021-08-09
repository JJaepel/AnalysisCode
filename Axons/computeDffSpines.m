function data = computeDffSpines(data, pt)

%remove slow baselin in raw F traces
if nargin <2
    pt = 99;
end

for i = 1:length(data.roi)
    raw_new = data.roi(i).rawF;
    
    % (1) cut off large events
    % raw_new(raw_new > std(raw)+median(raw)) = median(raw);
    
    %(2) 99 pt medfilt for low-pass trace
    % pt-order (window size),one dimensional median filter
    raw_new = medfilt1(raw_new,pt); 
    
    %(3) calc initial F value (median)
    % raw_new = cat(1,median(raw).*ones(100,1),raw_new);
    % raw_new = cat(1,median(raw).*ones(100,1),raw_new);
    % raw_new = prctfilt1(raw_new,90);
    % raw_new = raw_new(101:end-100);
    
    %(4) subtract and add back
    data.roi(i).dff = (data.roi(i).rawF - raw_new)./raw_new;

end