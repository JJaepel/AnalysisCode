function view_tiff_stack(data, metadata, data2, map)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tiff viewer
% 
% data - our stack
% metadata - information about imaging rate
% data2 - optional another stack with same dimensions
%         displayed next to first one
% map - optional change of colormap, default gray
%
% key control:
% 'space' - start movie
% 's' - scale to square
% 'm' - save as movie
% arrows up and down - speed up and slow down movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 2
    disp('no speed info');
    metadata = struct;
    metadata.Imaging.rate = 1;
end

if nargin < 3
    data2 = 0;
end

if nargin < 4
    map = 'gray';
end

hf=figure;
clf;
params.button_down=0;
params.movie_status=0;
params.frame_pause=1/(metadata.Imaging.rate); % pause in seconds between frames in playback
params.playback_spacing=1; % only every n-th frame is shown during playback default=1 (every frame)
params.as_ind=0; % index over area-selected for plotting activity

set(hf,'UserData',params);
dimensions = size(data);
dimensions2 = size(data2);

if length(dimensions2) > 2
    if dimensions(3) ~= dimensions2(3) || dimensions(2) ~= dimensions2(2)
        data2 = 0; 
        disp('Warning! Inputs do not match in size, only displaying first!');
    end
end

dimensions2 = size(data2);
offset = 100;

%set size of video
if length(dimensions2)>2
    set(hf,'Position',[offset offset dimensions(2)*2+offset 1.03*dimensions(1)+offset]);
else
    set(hf,'Position',[offset offset dimensions(2)+offset 1.03*dimensions(1)+offset]);
end
frame_ind=1;

%set contrast
template=data(:,:,round(size(data,3)/2));
low_contrast_lim=prctile(reshape(template,numel(template),1),10);
high_contrast_lim=prctile(reshape(template,numel(template),1),99.99);
params.clim=[low_contrast_lim high_contrast_lim];

%start with first frame
if length(dimensions2)>2
    params.h_f_ax=axes('position',[0 0.03 1 0.97],'color','k','clim',params.clim); 
    params.h_im_data=imagesc([data(:,:,frame_ind) data2(:,:,frame_ind)]);
else
    params.h_f_ax=axes('position',[0 0.03 1 0.97],'color','k','clim',params.clim);
    params.h_im_data=imagesc(data(:,:,frame_ind));
end

%add in text showing the time stamp
yl=ylim;xl=xlim;
params.h_txt=text(xl(2)/20,yl(2)/20,'0 s','fontsize',12,'color','w','fontweight','bold');

colormap(map);
h_sl_ax=axes('position',[0 0 1 0.03],'color','k');
xlim([0 1]);
ylim([0 1]);
hold on;

%add in positioned bar
sl_x_pos=0;
h_sl_bar=plot([1 1]*sl_x_pos,[0 1],'b','linewidth',15);

set(hf,'windowbuttondownfcn',@view_tiff_stack_buttdofcn);
set(hf,'windowbuttonupfcn',@view_tiff_stack_buttupfcn);
set(hf,'windowbuttonmotionfcn',{@view_tiff_stack_winmotfcn,h_sl_ax,h_sl_bar,data,data2, metadata});
set(hf,'KeyPressFcn' ,{@view_tiff_stack_keypress,h_sl_bar,data,data2,length(dimensions2)>2,metadata});
set(hf,'CloseRequestFcn',@close_fcn);
set(hf,'menubar','none','color','k');
set(params.h_f_ax,'clim',params.clim);
set(hf,'UserData',params);


function view_tiff_stack_winmotfcn(hf,e,h_sl_ax,h_sl_bar,data,data2, metadata)
params=get(hf,'UserData');
if params.button_down
    cp=get(h_sl_ax,'currentpoint');
    cp=cp(1);
    cp=min(cp,1);
    cp=max(cp,0);
    set(h_sl_bar,'Xdata',[1 1]*cp);
    if length(size(data2))>2
        set(params.h_im_data,'CData',[data(:,:,round(cp*(size(data,3)-1))+1) data2(:,:,round(cp*(size(data,3)-1))+1)]);
    else
        set(params.h_im_data,'CData',data(:,:,round(cp*(size(data,3)-1))+1));
    end
    set(params.h_txt,'string',[num2str(round(cp*(size(data,3)-1)/metadata.Imaging.rate)+1) ' s']);
end

function view_tiff_stack_buttdofcn(hf,e)
params=get(hf,'UserData');
params.button_down=1;
set(hf,'UserData',params);

function view_tiff_stack_buttupfcn(hf,e)
params=get(hf,'UserData');
params.button_down=0;
set(hf,'UserData',params);

function view_tiff_stack_keypress(hf,event,h_sl_bar,data,data2,twosided,metadata)
params=get(hf,'UserData');
nFrames = size(data,3);
spacing = 1/nFrames;
switch event.Character
    case 's'  % scale
        temp = get(hf);
        figLength = temp.Position(3);
        if twosided
            set(hf,'Position',[temp.Position(1) temp.Position(2) figLength 1.06*figLength/2]);
        else
            set(hf,'Position',[temp.Position(1) temp.Position(2) figLength 1.03*figLength*size(data,1)/size(data,2)]);
        end
    case 'm' % make a movie
        acqfps =1;
        [avi_fname,avi_path]=uiputfile('tiff_stack.avi','save stack as');
        writerObj = VideoWriter([avi_path avi_fname],'MPEG-4');
        writerObj.FrameRate = acqfps;
        open(writerObj);
        allFrames=isequal(input(['Do you want to save all individual frames (Y) or a time speed up?: '],'s'),'Y');
        if allFrames
            speed_factor = 1;
        else
            writerObj.FrameRate =23;
            speed_factor=input('Select speed to save movie at: ');
        end
        frame_bounds=input('Select [start_frame stop_fram], 0 for all frames: ');
        if frame_bounds==0
            frame_bounds=[0 nFrames];
        end
        for cp=frame_bounds(1)/nFrames:spacing*speed_factor:(frame_bounds(2)+1)/nFrames-spacing*speed_factor
            set(h_sl_bar,'Xdata',[1 1]*cp);
            if length(size(data2))>2
                set(params.h_im_data,'CData',[data(:,:,round(cp*(nFrames-1))+1) data2(:,:,round(cp*(nFrames-1))+1)]);
            else
                set(params.h_im_data,'CData',mean(data(:,:,round(cp*(nFrames-1)+1):min(size(data,3),round(cp*(nFrames-1))+speed_factor)),3));
            end
            set(params.h_txt,'string',[num2str(round(cp*(nFrames-1)/metadata.Imaging.rate)+1) ' s']);
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
        end
        close(writerObj);
    case 'a' % show activity of selected area
        params.as_ind=params.as_ind+1;
        color_ind='mgycrbw';
        as=ginput(2);
        hold on
        plot([as(2) as(2) as(1) as(1) as(2)],[as(4) as(3) as(3) as(4) as(4)],color_ind(params.as_ind));
        as_x=round(sort([as(3) as(4)]));
        as_y=round(sort([as(1) as(2)]));
        fig_pos=get(gcf,'position');
        if isfield(params,'sup_fig_h') && ishandle(params.sup_fig_h)
            set(params.sup_fig_h,'position',[fig_pos(1) fig_pos(2)-100-33 fig_pos(3) 100]);
        else
            params.sup_fig_h=figure('position',[fig_pos(1) fig_pos(2)-100-33 fig_pos(3) 100],'menubar','none','color','k');
            axes('position',[0 0 1 1],'color','k');
            hold on
        end
        prev_fig_handle=gcf;
        figure(params.sup_fig_h);
        try
            raw_act_trace=squeeze(mean(mean(data(as_x(1):as_x(2),as_y(1):as_y(2),:),2),1));
            plot(raw_act_trace,color_ind(params.as_ind));
            assignin('base','F',raw_act_trace);
        catch
            disp('try again - selection not valid');
        end
        axis tight;
        figure(prev_fig_handle);
        set(hf,'UserData',params);
    case 't' % show template
        fig_pos=get(gcf,'position');
        prev_fig_handle=gcf;
        if isfield(params,'sup_fig2_h') && ishandle(params.sup_fig2_h)
            set(params.sup_fig2_h,'position',[fig_pos(1)+fig_pos(3) fig_pos(2) fig_pos(3) fig_pos(4)]);
            figure(params.sup_fig2_h);
        else
            params.sup_fig2_h=figure('position',[fig_pos(1)+fig_pos(3) fig_pos(2) fig_pos(3) fig_pos(4)],'menubar','none');
        end
        axes('position',[0 0 1 1]);
        imagesc(mean(data,3));colormap gray;
% %         set(gca,'clim',params.clim/2)
        set(gca,'clim',params.clim/2)
        axis off
        figure(prev_fig_handle);
        set(hf,'UserData',params);
    case 'c' % clear
        hold off
        cp=get(h_sl_bar,'Xdata');
        cp=cp(1);
        axes(params.h_f_ax);
        params.h_im_data=imagesc(data(:,:,round(cp*(size(data,3)-1))+1));
        set(params.h_f_ax,'clim',params.clim,'xtick',[]);
        yl=ylim;xl=xlim;
        params.h_txt=text(xl(2)/20,yl(2)/20,num2str(round(cp*(size(data,3)-1)/metadata.Imaging.rate)+1),'fontsize',12,'color','w','fontweight','bold');
        if isfield(params,'sup_fig_h') && ishandle(params.sup_fig_h)
            close(params.sup_fig_h);
        end
        if isfield(params,'sup_fig2_h') && ishandle(params.sup_fig2_h)
            close(params.sup_fig2_h);
        end
        params.as_ind=0;
        set(hf,'UserData',params);
        
end

switch event.Key
    case 'rightarrow'
        if strcmp(event.Modifier,'shift')
            params.clim(1)=params.clim(1)*1.1;
            set(params.h_f_ax,'clim',params.clim)
            params
            set(hf,'UserData',params);
        else
            cp=get(h_sl_bar,'Xdata');
            cp=cp(1);
            cp = cp + spacing;
            if cp > 1, cp = 1; end
            set(h_sl_bar,'Xdata',[1 1]*cp);
            if length(size(data2))>2
                set(params.h_im_data,'CData',[data(:,:,round(cp*(size(data,3)-1))+1) data2(:,:,round(cp*(size(data,3)-1))+1)]);
            else
                set(params.h_im_data,'CData',data(:,:,round(cp*(size(data,3)-1))+1));
            end
            set(params.h_txt,'string',[num2str(round(cp*(size(data,3)-1)/metadata.Imaging.rate)+1) ' s']);
        end
    case 'leftarrow'
        if strcmp(event.Modifier,'shift')
            params.clim(1)=params.clim(1)*0.9;
            set(params.h_f_ax,'clim',params.clim)
            set(hf,'UserData',params);
        else
            cp=get(h_sl_bar,'Xdata');
            cp=cp(1);
            cp = cp - spacing;
            if cp < 0, cp = 0; end;
            set(h_sl_bar,'Xdata',[1 1]*cp);
            if length(size(data2))>2
                set(params.h_im_data,'CData',[data(:,:,round(cp*(size(data,3)-1))+1) data2(:,:,round(cp*(size(data,3)-1))+1)]);
            else
                set(params.h_im_data,'CData',data(:,:,round(cp*(size(data,3)-1))+1));
            end
            set(params.h_txt,'string',[num2str(round(cp*(size(data,3)-1)/metadata.Imaging.rate)+1) ' s']);
        end
    case 'uparrow'
        if strcmp(event.Modifier,'shift')
            params.clim(2)=params.clim(2)*1.1;
            set(params.h_f_ax,'clim',params.clim)
        else
            params=get(hf,'UserData');
            if params.frame_pause<=0.01
                params.playback_spacing=params.playback_spacing*2;
            else
                params.frame_pause=params.frame_pause*0.5;
            end
        end
            set(hf,'UserData',params);

    case 'downarrow'
        if strcmp(event.Modifier,'shift')
            params.clim(2)=params.clim(2)*0.9;
            set(params.h_f_ax,'clim',params.clim)
        else
            params=get(hf,'UserData');
            if params.playback_spacing>1
                params.playback_spacing=params.playback_spacing/2;
            else
                params.frame_pause=params.frame_pause*2;
            end
        end
            set(hf,'UserData',params);
    case 'space'   % play as movie
        params.movie_status=~params.movie_status;
        set(hf,'UserData',params);
        cp=get(h_sl_bar,'Xdata');
        cp=cp(1);
        while params.movie_status
            params=get(hf,'UserData');
            cp = cp + spacing*params.playback_spacing;
            if cp > 1, cp = 0; end
            set(h_sl_bar,'Xdata',[1 1]*cp);
            if length(size(data2))>2
                set(params.h_im_data,'CData',[data(:,:,round(cp*(size(data,3)-1))+1) data2(:,:,round(cp*(size(data,3)-1))+1)]);
            else
                set(params.h_im_data,'CData',data(:,:,round(cp*(size(data,3)-1))+1));
            end
            set(params.h_txt,'string',[num2str(round(cp*(size(data,3)-1)/metadata.Imaging.rate)+1) ' s']);
            pause(params.frame_pause);
        end
end

function close_fcn(hf,e)
params=get(hf,'UserData');
if isfield(params,'sup_fig_h') && ishandle(params.sup_fig_h)
    close(params.sup_fig_h);
end
if isfield(params,'sup_fig2_h') && ishandle(params.sup_fig2_h)
    close(params.sup_fig2_h);
end
delete(hf)






