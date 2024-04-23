%% Master script for going through the data
%For each cell (and slice), there is a folder in InputAnalysis, with the
%following substructure
% A - 2p Imaging: includes all the analyzed 2p data
% B - Confocal: includes all confocal information
%               - Dendrites: Dendrite Tracing with SNT
%               - Spines: Spine & Soma marking with RoiManage
%               - Dendrite reconstruction: Tracing of all dendrites with
%               its spines
%               - ....tif: Confocal stack of cell
%               - ....jpg: Projection of confocal stack of cell
% C - STED: includes all STED imaging files
% D - Alignments: includes files for alignment
%               - Alignement 2p: Files for 2p confocal grouping
%               - Alignment STED: Illustrator files for STED alignments
%               - STED coverage: Files for STED coverage
% E - Analysis: includes all analysis files for each animal
%               - 01_Morphological analysis
%               - 02_Functional cell reconstruction
%               - 03_Anatomical input analysis
%               - 04 Functional input analysis


%% Step 0: Run the parser to get the most up to date info
cellInfo = animalParser;
i = 1; %select your cell of interest or in later instances set that to all cells

%% Step 1: Analyse the individual parts

%A) 2p imaging --> generate completeCell.mat
%1) run through all experiments, make sure to mark the required experiments in the
%file SpinePerAnimal.xslx
SpineImagingAnalysis([], 'F2712')
%2) combine both to generate the completeCell.mat
spineToCellComputations('F2712')

%B) Confocal tracing -> generate cellReconstruction.mat
%1) Trace the dendrites with ImageJ -> Neuroanatomy -> SNT, make sure that
%each branch starts at the soma -> Save as
%{Animal}_{Cell}_{Slice}_DendriteXX-yyy.swc
%2) Mark all the spines for each branch with ImageJ -> ROIManager and select
%the CellMagicWand tool for it -> Save as
%RoiSet_{Animal}_{Cell}_{Slice}_DendriteXX.zip
%3) Markt the soma with ImageJ -> ROIManager and select the CellMagicWand tool
%for it -> Save as {Animal}_{Cell}_{Slice}_Soma.roi
%4) Put everything together to generate the cell -> cellReconstruction.mat
confocalTracingtoPlot(cellInfo, i)
%5) Analyze the cell
AnatomicalAnalysisCell('F2688') %also add cellName if multiple cells per animal

%C) STED imaging
% 1) Do the imaging & 2) Look through which ones are inputs

%% Step 2: Align 2p & STED imaging to confocal imaging

%A) 2p and confocal imaging -> generate FunctionalCellReconstruction.mat
%1) Rough mathing with powerpoint and adding them on top of each other in
%illustrator
%2) Fine alignment and grouping
ConfTo2pMatching
%3) Run the functional analysis of the cell
FunctionalConfocalRepresentation('F2688') %also add cellName if multiple cells per animal

%B) STED and confocal -> generate InputCellReconstruction.mat
%1) Align all STED dendrites to the confocal on the illustrator files and
%mark the inputs there
%2) Using the ROI manager and excel, mark which confocal ROIs are inputs in
%the excel file STED Inputs.xsxl for each branch individually
%3) Select for each dendrite whether it was covered by STED and to what
%percentage
STEDcoverage_V3
%4) Run the anatomical input analysis to combine both modalities and get
%some first anatomical measurements -> Input
AnatomicalInputAnalysis('F2688') %also add cellName if multiple cells per animal

%% Step 3: Combine information of both alignments to get functional properties of the input
MultiModalAnalysisNew('F2688') %also add cellName if multiple cells per animal

%% Step 4: Combine all cells to look at functional properties
MultiCellAnalysis

%% Step 5: Compare input sources
MultiCellInputAnalysis

%% Step 6: Rerun all steps automatically, assuming that the manual ones are done
UpdateInputAnalysis
