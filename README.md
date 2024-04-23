# AnalysisCode

#Spine Analysis & Combining modalities

Step 1: Analyse functional data 
    1.1 Of cell -> CellAnalysis.m
    1.2 Of spines -> SpineAnalysis.m
    1.3 Compute relationship between cell and spines -> spineToCellCompuations.m

Step 2: Divide cell into its major dendrites

Step 3: Match in vivo 2p and confocal data
    3.1: Roughly align 2p and confocal in power point -> which experiments belong to which dendrite?
    3.2: For each dendrite, match confocal data with all relevant experiments -> ConfTo2pMatching.mlapp
    3.3: Build the funcitonal representation of the cell -> FunctionalConfocalRepresentation.m

Step 4: Match Confocal and STED data
    3.1: Prepare confocal data
            3.1.1 Trace cell in ImageJ with Neurite Tracer: mark by dendrite Nr, save in Z:/Juliane/Data/Confocal/ and make a presentation
            3.1.2 For each dendrite, mark the spines with ROI manager, save in  Z:/Juliane/Data/Confocal/
            3.1.3 Combine all data -> confocalTracingtoPlot.m, saved as cellReconstruction.mat in Z:/Juliane/Data/Confocal/
    3.2: STED Data
            3.2.1 Match all STED dendrites to the confocal dendrites in Illustrator
            3.2.2 Mark the inputs in all STED dendrites in a file
    3.3: Matching confocal & STED
            3.3.1 Open stack and for each dendrite the ROI files to see the positions
            3.3.2 Mark the ROINr for each input for each dendrite in an excel file
    
Step 5: Analyze anatomical data -> AnatomicalAnalysisCell.m

Step 6: Combine anatomical and functional data -> 

