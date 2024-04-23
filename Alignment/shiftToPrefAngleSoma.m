function tuningCurveNew = shiftToPrefAngleSoma(prefInd, tuningCurve)
%Shifts tuning curve so that the prefered orientation or direction is in
%the center of it
%Input:
% - prefInd: Index of the prefered orietation
% - tuningCurve: values of responses to all orientations or directions, e.
% g. a 1xnumber of orientations or directions vector
%
%Output:
% - tuningCurveNew: shifted tuning curve with the prefered
% orientation/direction in the center, e. g. a 1xnumber of orientations or
% directions + 1 vector

%depending on whether the prefered stimulus is in the first or second half
%of the tuning curve
if prefInd < length(tuningCurve)/2+1
    %how many stimuli from the end are added before the old tuning curve
    addBefore = length(tuningCurve)/2-prefInd;
    %until when do you add the tuning curve
    addUntil =  length(tuningCurve)/2 + prefInd;
    tuningCurveNew = [tuningCurve(end-addBefore:end), tuningCurve(1:addUntil)];
else
    add = prefInd - length(tuningCurve)/2;
    tuningCurveNew = ([tuningCurve(add:end), tuningCurve(1:add)]);
end
