function analysis = normResponse(analysis, data, field, ROI)

%get trace data
Response = data.roi(ROI).dff;

%normalize Response
xmin=min(Response(:));
xmax=max(Response(:));
normResponse = (Response-xmin)/(xmax-xmin);

%save in analysis struct
analysis.(field).roi(ROI).normResponse = normResponse;
