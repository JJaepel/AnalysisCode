function [data]= LoadRoisS2p(analysisParams)
    basedirectory = analysisParams.baseDirectory;
    expt_id = analysisParams.expID;
    if analysisParams.level == 1
        Suite2pFile = [basedirectory filesep expt_id filesep 'suite2p' filesep 'combined' filesep 'Fall.mat'];
    else
        Suite2pFile = [basedirectory filesep expt_id filesep 'suite2p' filesep 'plane0' filesep 'Fall.mat'];
    end
    Suite2p = load(Suite2pFile);
    load(Suite2pFile);
    disp('loaded Suite2pFile');
    cell_selector = logical(Suite2p.iscell(:,1));
    data.roi = [];
    counter = 1;
    if analysisParams.level == 1
        for i = 1:length(cell_selector)
            if cell_selector(i) == 1
                    data.roi(counter).plane = Suite2p.stat{i}.iplane+1;
                if Suite2p.stat{i}.iplane == 0 || Suite2p.stat{i}.iplane == 1
                    data.roi(counter).xPos = Suite2p.stat{i}.med(2);
                    data.roi(counter).yPos = Suite2p.stat{i}.med(1);
                    data.roi(counter).mask = [Suite2p.stat{i}.xpix; Suite2p.stat{i}.ypix]';
                elseif Suite2p.stat{i}.iplane == 2
                    data.roi(counter).xPos = Suite2p.stat{i}.med(2)-1024;
                    data.roi(counter).yPos = Suite2p.stat{i}.med(1)+512;
                    data.roi(counter).mask = [Suite2p.stat{i}.xpix-1024; Suite2p.stat{i}.ypix+512]';
                elseif Suite2p.stat{i}.iplane == 3
                    data.roi(counter).xPos = Suite2p.stat{i}.med(2)+512;
                    data.roi(counter).yPos = Suite2p.stat{i}.med(1);
                    data.roi(counter).mask = [Suite2p.stat{i}.xpix+512; Suite2p.stat{i}.ypix]';
                end
                data.roi(counter).name = i;
                data.roi(counter).rawF = double(Suite2p.F(i,:));
                data.roi(counter).spks = double(Suite2p.spks(i,:));
                counter = counter +1;
            end
        end
        if size(Suite2p.ops.meanImg,2) == 1024
            plane0 = Suite2p.ops.meanImg(1:512,1:512); plane1 = Suite2p.ops.meanImg(1:512, 513:1024); plane2 = Suite2p.ops.meanImg(513:1024,1:512); plane3 = Suite2p.ops.meanImg(513:1024,513:1024);
        elseif size(Suite2p.ops.meanImg,2) == 1536
            plane0 = Suite2p.ops.meanImg(1:512,1:512); plane1 = Suite2p.ops.meanImg(1:512, 513:1024); plane2 = Suite2p.ops.meanImg(1:512, 1025:1536); plane3 = Suite2p.ops.meanImg(513:1024,1:512);
        end
        template = [plane0 plane1; plane2 plane3];
    else
        for i = 1:length(cell_selector)
            if cell_selector(i) == 1
                data.roi(counter).xPos = Suite2p.stat{i}.med(2);
                data.roi(counter).yPos = Suite2p.stat{i}.med(1);
                data.roi(counter).mask = [Suite2p.stat{i}.xpix; Suite2p.stat{i}.ypix]';
                data.roi(counter).name = i;
                data.roi(counter).rawF = double(Suite2p.F(i,:));
                data.roi(counter).spks = double(Suite2p.spks(i,:));
                counter = counter +1;
            end
        end
        template = Suite2p.ops.meanImg(1:512,1:512);
    end
    data.template = template./prctile(template(:),99.9);
    clear Suite2p
    disp('loaded data from Suite2p')
end