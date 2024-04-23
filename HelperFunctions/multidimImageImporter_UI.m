function multidimImageImporter_UI(mainObj, handles)
% Import multidimensional images to workspace using BioFormats
% DESCRIPTION: The GUI uses BioFormats matlab files to import multidimensional 
% arrays (XYZCT) and OME-XML metadata to the MATLAB workspace.  
% See https://docs.openmicroscopy.org/bio-formats/5.7.0/developers/matlab-dev.html
% for further documentation regarding 'bfmatlab' functions.
%
% OUTPUT:   Structure (I) containing three fields...
%               I.name: name of the image file with XYZCT dimensions
%               I.img: image array
%               I.meta: OME-XML metadata
%           New images can be appended to the structure in seqential elements
%           (i.e. I(1), I(2), .... I(n))
% 
% NAVIGATE GUI: 
%   1. Select 'Import' found under the File menu.
%   2. Navigate to your image file of interest and open.  
%   (the list_box will be populated with metadata for images including the 
%   series #, X, Y, Z, channel and time dimensions).  If image file 
%   contains multiple series, each series will be listed seperately.
%   3. Click on an image in the list_box to import image data 
%   directly to the workspace. 
%   4. Select as many images as you wish to be appended to the
%   structure, I.
%   5. Select 'Clear' under the File menu to reset list_box menu and start
%   fresh.
%
% EXAMPLES:
% Import images in sequential order then run script to view image arrays
% 1. mri-stack.tif (3D array) 
% 2. confocal-series.tif (4D array) 
% 3. mitosis.tif (5D array) then run script to view image array
% 
%       % view mri-stack.tif in Z
%       implay(I(1).img)
%       % seperate red and gree channels from confocal-series.tif then 
%       %view in Z
%       red_channel = I(2).img(:,:,:,1);
%       implay(red_channel)
%       green_channel = I(2).img(:,:,:,2);
%       implay(green_channel)
% 
%       % process mitosis.tif to view multicolor timelapse 
%       maxI = max(I(3).img, [], 3); % max project the z-dimension
%       maxI = squeeze(maxI); % squeeze the array to remove single dimension
%       b_plane = zeros(196,171, 1, 51, 'uint16');  % create emtpy blue channel array
%       maxI = cat(3, maxI, b_plane);   % concatenate blue channel
%       maxI = im2double(maxI); % convert image to type double 
%       mov = zeros(196,171,3,51);  % empty array to populate 
% 
%       % adjust scaling of images for better visualization
%       % perform imadjust on R and G channels
%       for c = 1:2
%           for t = 1:size(maxI, 4)
%               img = maxI(:,:,c,t);
%               mov(:,:,c,t) =  imadjust(img, [min(img(:)) max(img(:))]);
%           end
%       end
%       implay(mov)
%
% NOTE: 
%   Create a new figure window to display image as outlined below...
%       figure; imshow(I(i).img) %to display image
%
% Comments and suggestions are greatly appreciated!
    N = 300;
    M = 300;
    
    mainObj = figure('Units', 'Pixels',  ...
        'Position', [100 100 N M], 'Menubar', 'none',...
        'NumberTitle','off', 'Name', 'Image Array Importer');
    handles = guihandles(mainObj);
    
    menu_File = uimenu('Label', 'File');
    uimenu('Parent', menu_File, ...
        'Label', 'Import', ...
        'Callback', {@menu_Import_callback, handles});
    
    handles.reset = uimenu('Parent', menu_File, ...
    'Label', 'Clear list', 'Enable', 'off', ...
    'Callback', {@menu_resetListBox_callback, handles});
    
    handles.reader = [];
    handles.directory = [];
    handles.list_str = {};
    handles.structIdx = 0;
    
    handles.list_box = uicontrol(mainObj, ...
        'style', 'listbox', 'string', handles.list_str, ...
        'HorizontalAlignment', 'left', 'Position', [0 0 N M]);
    
    set(handles.list_box, 'callback', {@list_box_callback, handles});
    % access guidata by calling 'handle = guidata(figure(1))'
    guidata(mainObj, handles)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function list_box_callback(mainObj, event_data, handles)
    handles = guidata(gcbo);
    
    list_index = get(handles.list_box, 'value');
    img_dir = handles.directory{list_index};
    reader = bfGetReader(img_dir);
    omeMeta = reader.getMetadataStore();
    series = omeMeta.getImageCount(); 
    handles.structIdx = handles.structIdx + 1;
    
    % s is the index for series
    for s = series
            % define series of interest and acquire needed dimensions
            reader.setSeries(s-1);
            X = omeMeta.getPixelsSizeX(s-1); X = str2num(X);
            Y = omeMeta.getPixelsSizeY(s-1); Y = str2num(Y);
            Z = omeMeta.getPixelsSizeZ(s-1); Z = str2num(Z);
            T = omeMeta.getPixelsSizeT(s-1); T = str2num(T);
            C = omeMeta.getChannelCount(s-1);
            title = omeMeta.getImageName(s-1);
            if isempty(title) == 1
                title = 'No name';
            else
                title = title.toCharArray';
            end
            B = omeMeta.getPixelsType(s-1);
            B = char(B);
            str = [title, ' ', num2str(Y), 'x', num2str(X), 'x' num2str(Z),...
                '; (', num2str(C),'C x ', num2str(T),'T)'];
            M = zeros(Y, X, Z, C, T, B);
        % c is the index for color
        for c = 1:C
            % t is the index for time 
            for t = 1:T
                % z is the index for Z plane
                for z = 1:Z
                    % create multidimensional matrix one image plane at a time
                    iPlane = reader.getIndex(z - 1, c - 1, t - 1) + 1;
                    M(:, :, z, c, t) = bfGetPlane(reader, iPlane);
                end        
            end
        end
    end
    
    handles.I(handles.structIdx).name = str;
    handles.I(handles.structIdx).img = M;
    
    omeXML = char(omeMeta.dumpXML());
    % regular expression to isolate metadata from long string
    expression = '<[^>]*>';
    matches = regexp(omeXML,expression,'match');
    handles.I(handles.structIdx).meta = cell2table(matches', 'VariableNames', {'MetaData'});
    assignin('base', 'I', handles.I);
    
guidata(mainObj, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function menu_Import_callback(mainOjb, event_data, handles)
    handles = guidata(gcbo);
    % navigate to image file containing data to import
    [filename path] = uigetfile('*.*');
    % use BioFormats functions to import image and metadata 
    reader = bfGetReader(fullfile(path, filename));
    omeMeta = reader.getMetadataStore();
    
    % empty cell to populate with image series metadata that will appear in
    % list_box
    str = {};
    series = omeMeta.getImageCount(); 
    for i = 1:series
        X = omeMeta.getPixelsSizeX(i-1); X = str2num(X);
        Y = omeMeta.getPixelsSizeY(i-1); Y = str2num(Y);
        Z = omeMeta.getPixelsSizeZ(i-1); Z = str2num(Z);
        T = omeMeta.getPixelsSizeT(i-1); T = str2num(T);
        C = omeMeta.getChannelCount(i-1);
        title = omeMeta.getImageName(i-1);
        if isempty(title) == 1
            title = 'No name';
        else
            title = title.toCharArray';
        end
        str{i, 1} = [title, ' ', num2str(X), 'x', num2str(Y),'x' num2str(Z),...
            '; (', num2str(C),'C x ', num2str(T),'T)'];
    end
 
    handles.directory = [handles.directory; repmat(cellstr(fullfile(path,filename)), series, 1)];
    handles.list_str = [handles.list_str; str];
    set(handles.list_box, 'string', handles.list_str);
    set(handles.reset, 'enable', 'on');
    guidata(mainOjb, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function menu_resetListBox_callback(mainObj, event_data, handles) 
    
    handles = guidata(gcbo);
    
    handles.directory = [];
    handles.list_str = {};
    handles.structIdx = 0;
    handles.I = [];
    set(handles.list_box, 'value', 1);
    set(handles.list_box, 'string', handles.list_str);
    
    set(handles.reset, 'enable', 'off');
    
    guidata(mainObj, handles)