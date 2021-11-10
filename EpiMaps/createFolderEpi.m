function analysisParams = createFolderEpi(analysisParams)

if analysisParams.server == 0
    drive = 'F:\';
else 
    drive = 'Z:\Juliane\';
end

EpiDir = [drive 'Data\Epi\'];
Sp2Dir = [drive 'Data\Spike2Data\'];
SaveDir = [drive 'Data\ImageAnalysis\'];

analysisParams.EpiDirectory = [EpiDir filesep analysisParams.animal filesep analysisParams.expID filesep];
analysisParams.Sp2dDirectory = [Sp2Dir analysisParams.animal filesep analysisParams.sp2ID filesep];
analysisParams.saveDirectory = [SaveDir analysisParams.animal filesep analysisParams.expID filesep];

if ~exist(analysisParams.saveDirectory, 'dir')
    mkdir(analysisParams.saveDirectory);  
end