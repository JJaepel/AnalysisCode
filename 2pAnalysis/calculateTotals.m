function All_OFFx=calculateTotals(ROI,NumROIs,varargin)
All_OFFx=[];
if nargin>2
    numIndex = find(cellfun('isclass', varargin(1:end), 'char'));
    if isfield(ROI,(varargin{numIndex}))
        %    data= getfield(ROI,(varargin{numIndex}));
        for nr=1:NumROIs
            if isempty(getfield(ROI(nr),(varargin{numIndex})))
                All_OFFx(nr)=0;
            else
                All_OFFx(nr)= getfield(ROI(nr),(varargin{numIndex}));
            end
        end
    end
end
