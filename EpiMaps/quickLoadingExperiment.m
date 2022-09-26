close all
clear all

% select data
%animal = 'F2532_2021-06-23'; expt_id =1; 
%animal = 'F2532_2021-06-23'; expt_id =12; % lots of movement!!!
%animal = 'F2537_2021-06-30'; expt_id =3;
%animal = 'F2538_2021-07-02'; expt_id =1;
%animal = 'F2563_2021-08-23'; expt_id =4;
%animal = 'F2564_2021-08-26'; expt_id =2;
%animal = 'F2569_2021-09-03'; expt_id =3;
%animal = 'F2570_2021-09-08'; expt_id =3;
%animal = 'F2573_2021-09-14'; expt_id =3;
%animal = 'F2574_2021-09-13'; expt_id =3;
%animal = 'F2636_2022-02-14'; expt_id =2;
%%animal = 'F2636_2022-02-14'; expt_id =14;
%animal = 'F2632_2022-02-16'; expt_id =3;
%%animal = 'F2635_2022-02-17'; expt_id =1;
%animal = 'F2635_2022-02-17'; expt_id =13;
%animal = 'F2656_2022-03-25'; expt_id =2;
%%animal = 'F2656_2022-03-25'; expt_id =13;
animal = 'F2657_2022-03-28'; expt_id =1;
%%animal = 'F2657_2022-03-28'; expt_id =14;

downsample = 2;

EpiDir = 'Z:\Juliane\Data\Epi\';
Sp2Dir = 'Z:\Juliane\Data\Spike2Data\';
saveDir = 'F:\Data\Movies';

EpiDirectory = [EpiDir filesep animal filesep 'tseries_' num2str(expt_id) filesep];
if expt_id > 9
    Sp2dDirectory = [Sp2Dir animal filesep 't000' num2str(expt_id) filesep];
else
    Sp2dDirectory = [Sp2Dir animal filesep 't0000' num2str(expt_id) filesep];
end
try
    metadata.Imaging=LoadFrameTimes(Sp2dDirectory);
catch
    metadata.Imaging.rate = 14.9993;
end

% load data
t0=tic;
disp('Loading data from tif files')
rawF = readingImagingDataMaskRestricted(EpiDirectory, downsample, 1);
toc(t0)
disp('Normalizing stack')
dff = NormalizeStack(metadata,rawF);

%view Data
view_tiff_stack(dff)