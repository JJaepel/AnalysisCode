function [decodingData] = decoderOri(analysis, data, field)
% Computs performance of a decoder based on significantly responding cells,
% standardly using 50 cells per trial
% Based on script from David Whitney
%
% Input:
% - data: structure containing for each ROI (data.roi) the dff trace
% - analysis: structure containing for each ROI (analysis.roi) the
% stim response trace and whether they are significant responders
%
% Ouput:
% - decodingData: decoder performance over all trials (nBootstraps)
%


    interval = (analysis.(field).windowStart:analysis.(field).windowStop)';
    nCells = length(data.roi);
    [nConds, nTrials, ~] = size(analysis.(field).roi(1).stimResponseTrace);
    isSig     = false(nCells,1);
    cellResps = zeros(nCells,(nConds-1)/2,2*nTrials);
    for n = 1:nCells
        dff = data.roi(n).(field);
        trialResps = analysis.(field).roi(n).stimResponseTrace;
        trialResps = nanmean(trialResps(1:(end-1),:,interval),3); %this is for decoding the direction, size numberDir x numTrials
        OriTrialResps = cat(2,trialResps(1:(nConds-1)/2,:),trialResps(((nConds-1)/2+1):end,:)); %this is for decoding the orientation, size numberOri x numTrials
        zscoredOriResps = (OriTrialResps-nanmean(dff))./nanstd(dff);
        cellResps(n,:,:) = zscoredOriResps; 
        isSig(n) = analysis.(field).roi(n).isResponseSignificant;
    end

    % Setup MATLAB anomynous functions to evalute decoding on real/shuffled data
    decodingStyle = 'similarity'; % 1-using similarity metric from Vince et al.,2-using correlation similarity
    switch decodingStyle
        case 'similarity'
            decoder = @(x,y) nansum(x(:).*y(:))./(sqrt(nansum(x(:).^2))*sqrt(nansum(y(:).^2)));
        case 'correlation'
            decoder = @(x,y) nansum(zscore(x(:)).*zscore(y(:)));
    end
    collapse3dMatrix = @(matrix,matrixSize) reshape(matrix,[matrixSize(1) matrixSize(2)*matrixSize(3)]);
    shuffle2dMatrix  = @(matrix) matrix(:,randperm(size(matrix,2)));
    restore3dMatrix  = @(matrix,matrixSize) reshape(matrix,[matrixSize(1) matrixSize(2) matrixSize(3)]);
    shuffle3dMatrix  = @(matrix) restore3dMatrix(shuffle2dMatrix(collapse3dMatrix(matrix,size(matrix))),size(matrix));

    % Process data with loops
    nBootstraps = 1000;
    nCellsUsed  = 50; 
    rng(1);
    medianTemplateResps = nanmedian(cellResps,3);

    decodingData  = {};
    for isShuffle = 0:1
        decoderPerformance  = zeros(nBootstraps,1);
        parfor(loop = 1:nBootstraps) % Can make a par for loopp
            % Pick a random subset of cells
            cellSubset = ceil(nCells*rand(nCells,1));
            cellSubset = cellSubset(1:nCellsUsed); 

            % Apply decoder on the selected population responses to each
            % stimulus trial
            if isShuffle
                x = shuffle3dMatrix(cellResps);
                y = nanmedian(shuffle3dMatrix(cellResps),3);
                %disp('Shuffled data')
            else
                x = cellResps;
                y = medianTemplateResps;
                %disp('Actual Data')
            end
            decodingMatrix = zeros((nConds-1)/2,nTrials,(nConds-1)/2);
            for cond = 1:(nConds-1)/2
                for trial = 1:nTrials
                    for refCond = 1:(nConds-1)/2
                        decodingMatrix(cond,trial,refCond) = decoder(x(cellSubset,cond,trial),y(cellSubset,refCond));
                    end
                end
            end

            % Evaluate decoding performance across all stimulus trials
            [~,predictedTrialLabel] = max(decodingMatrix,[],3);
            actualTrialLabel = repmat(1:(nConds-1)/2,[nTrials 1])';
            correctTrials = predictedTrialLabel==actualTrialLabel;  
            decoderPerformance(loop) = mean(correctTrials(:));
        end
        decodingData{isShuffle+1}=decoderPerformance; % Save the decoder for actual data and shuffle
    end
end