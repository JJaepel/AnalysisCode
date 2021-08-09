function semi_manual_axon_ROIs(datapath, filename)

avg = mijread([datapath filename]);
MIJ.run('Duplicate...', 'title=actthresh');
MIJ.run("8-bit");
MIJ.run("Bandpass Filter...", "filter_large=40 filter_small=3 suppress=None tolerance=5 autoscale saturate");
MIJ.run("Threshold...", "dark background");
MIJ.run("Convert to Mask");
%f = figure('Position', [40 400 210 50],'menuBar', 'none', 'name', 'execution paused');
%h = uicontrol('Position',[10 10 190 30],'String','Threshold OK?',...
%                'Callback','uiresume(gcbf)');
%uiwait(gcf);
%close(f);

MIJ.run('Watershed');
pause(0.5);
MIJ.run('Save', ['path=[' datapath filesep 'act_thresh.tif]'])
pause(0.5);
MIJ.run("Analyze Particles...", "size=50-500 circularity=0.25-1.00 show=Masks display exclude clear add");
pause(0.5);
MIJ.run('Save', ['path=[' datapath filesep 'final_thresh.tif]'])
pause(0.5);

%f = figure('Position', [40 400 210 50],'menuBar', 'none', 'name', 'execution paused');
%h = uicontrol('Position',[10 10 190 30],'String','Save and Next Experiment?','Callback','uiresume(gcbf)');
%uiwait(gcf);
%MIJ.run('Save ROIs [k]')
MIJ.run('AxonROIs')
MIJ.run('Close All');
%close(f);