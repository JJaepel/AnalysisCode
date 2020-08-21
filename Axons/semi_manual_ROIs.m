function semi_manual_ROIs(datapath, filename)

avg = mijread([datapath filename]);
MIJ.run('Duplicate...', 'title=actthresh');
MIJ.run("8-bit");
MIJ.run("Bandpass Filter...", "filter_large=40 filter_small=3 suppress=None tolerance=5 autoscale saturate");
MIJ.run("Threshold...");
f = figure('Position', [40 400 210 50],'menuBar', 'none', 'name', 'execution paused');
h = uicontrol('Position',[10 10 190 30],'String','Threshold OK?',...
                'Callback','uiresume(gcbf)');
%WinOnTop(gcf);
uiwait(gcf);
close(f);

MIJ.run('Watershed');
MIJ.run('Save', ['path = [' datapath filesep 'avg_thresh.tif]']);
MIJ.run("Analyze Particles...", "size=50-500 circularity=0.50-1.00 show=Masks display exclude clear add");
MIJ.run('Save', ['path=[' datpath filesep 'final_thresh.tif]'])

f = figure('Position', [40 400 210 50],'menuBar', 'none', 'name', 'execution paused');
h = uicontrol('Position',[10 10 190 30],'String','Save and Next Experiment?','Callback','uiresume(gcbf)');
WinOnTop(gcf);
uiwait(gcf);
MIJ.roiManager('Save', ['path = [' datapath filesep 'ROI.zip]']);
MIJ.run('Close All');
close(f);