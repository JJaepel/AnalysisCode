function [data] = computeDff(data, sourceSignal, sourceBaseline, fieldTarget)
    for i=1:length(data.roi)
        data.roi(i).(fieldTarget)= (data.roi(i).(sourceSignal)-data.roi(i).(sourceBaseline))./abs(data.roi(i).(sourceBaseline));
        %data.roi(i).(fieldTarget)= data.roi(i).(sourceSignal)-data.roi(i).(sourceBaseline);
    end
end