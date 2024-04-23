function [dff]=filterBaseline_dFcomp(raw,pt)

% Calculates the dff by getting a sliding raw Traces and usinig that as F0
%
%Input:
% - raw: raw trace
% - pt: what percentage to use as filter, if no set it is 99
%
%Ouput:
% - dff: dff trace

if nargin<2
    pt= 99;
end

%% Step 1: Do a medfilter
slidingF0 = raw;

%(1) Optional: cut off large events 
% raw_new(raw_new > std(raw)+median(raw)) = median(raw);

%(2) 99 pt medfilt for low-pass trace
slidingF0 = medfilt1(slidingF0,pt); 

%(3) Optional: Calc initial F value (median)
% raw_new = cat(1,median(raw).*ones(100,1),raw_new);
% raw_new = cat(1,median(raw).*ones(100,1),raw_new);
% raw_new = prctfilt1(raw_new,90);
% raw_new = raw_new(101:end-100);

%% Step 2: Get dff: dff = (F-F0)/F0
dff = (raw - slidingF0)./slidingF0;



