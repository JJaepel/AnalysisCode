function constructInputCell(allSpines, allBranches, inputType, property)

% Constructs a summed cell of all dendrites that have an input from that
% type and marks them as the property that is desired
%
% Input:
% - Spines: Data of the spines, including their functional data
% - Dendrites: Tracing of the dendrites
% - inputType: is it A19 or V1?
% - property: which property should be plotted? E. g. deltaOri, OSI, ...
% 
% Steps:
%   1) Find all input spines - choose with kind of requirements based on
%   property and change the coordinates relative to the cell body so that
%   the soma is in the center at 0/0
%   2) Find the corresponding correct dendrites plus all their way to the
%   cell body, change the coordinates relative to the cell body so that
%   the soma is in the center at 0/0 and plot them
%   3) Plot the soma in white
%   4) Plot the spines and their property, changing coordinates relative to
%   the cell body so that the soma is in the center at 0/0
%   5) Add colorbar as legend and a title
%
% Output: 
% - figure: constructed cell with all inputs from a specific cell
%
% Written by Juliane Jaepel
% Max Planck Florida Institute for Neuroscience
% Version 1.0: March, 2024

coc = cbrewer('qual', 'Set1', 8);

%% Step 1: Find all relevant input spines

%depending on input type, get the inputROIs  
switch inputType
    case 'A19'
        inputNr = 1;
    case 'V1'
        inputNr = 2;
end

inputROIs = find([allSpines.input] == inputNr); %which ones are the right input

%get all relevant spines
TwoPROIsNR = find([allSpines.TwoPMatch] ==1); %which ones are 2p spines
goodROIs = find([allSpines.good] == 1); %which of those are responsive
oriTwoPROIs = find([allSpines.OSI] > 0.1); %all oriselectROIs, size = allfuncSpines
oriSelect = TwoPROIsNR(oriTwoPROIs);
oriGood = intersect(oriSelect, goodROIs);
dirTwoPROIs = find([allSpines.DSIvect] > 0.1); %all oriselectROIs, size = allfuncSpines
dirSelect = TwoPROIsNR(dirTwoPROIs);
dirGood = intersect(dirSelect, goodROIs);

%intersect with input
goodInputs = intersect(goodROIs, inputROIs);
oriInputs = intersect(oriGood, inputROIs);
dirInputs = intersect(dirGood, inputROIs);

switch property
    case 'deltaOri'
        inputs = oriInputs;
   case 'deltaDir'
        inputs = dirInputs;
    case 'OSI'
        inputs = goodInputs;
    case 'DSI'
        inputs = goodInputs;
    case 'struct'
        inputs = inputROIs;
    case 'CellNr'
        inputs = inputROIs;
end

%% Step 2: Find the corresponding correct dendrites and plot them
hold on
for i = 1:length(inputs)
    % which cell is it and what is the soma pos
    cellNr = allSpines(inputs(i)).CellNr;
    
    %find the dendriteNr of the input
    dendriteNr = allSpines(inputs(i)).denOnBranch;
    branchNr = allSpines(inputs(i)).Branch;
    %get which position it is on on the allBranches
    dendrite = intersect(find([allBranches.CellNr] == cellNr), intersect(find([allBranches.Branch] == branchNr), find([allBranches.denOnBranch] == dendriteNr)));
    %get the dendriteOrder
    dendriteOrder = allBranches(dendrite).BranchOrder;
    
    %plot the dendrite
    if contains(property, 'CellNr')
        plot(allBranches(dendrite).normCoord(:,1),allBranches(dendrite).normCoord(:,2), 'LineWidth',1.5, 'color', coc(cellNr,:))
    else 
        plot(allBranches(dendrite).normCoord(:,1),allBranches(dendrite).normCoord(:,2), 'LineWidth',1.5, 'color', 'black');
    end
    
    %while the dendrite Order is not one, meaning it is not the closest to
    %the soma, find the previous dendrites and plot them
    while dendriteOrder > 1
        
       %get the new dendrite nr and hwere it is on the allBranches 
       newDendriteNr = allBranches(dendrite).startPoint;
       if newDendriteNr == 0
           dendriteOrder = 1;
       else
           newDendrite = intersect(find([allBranches.CellNr] == cellNr), intersect(find([allBranches.Branch] == branchNr), find([allBranches.denOnBranch] == newDendriteNr)));

           %set this as the "old" dendrite and get the branchOrder
           dendrite = newDendrite;
           dendriteOrder = allBranches(dendrite).BranchOrder;

           %plot the dendrite
           if contains(property, 'CellNr')
               plot(allBranches(dendrite).normCoord(:,1),allBranches(dendrite).normCoord(:,2), 'LineWidth',1.5, 'color', coc(cellNr,:))
           else 
               plot(allBranches(dendrite).normCoord(:,1),allBranches(dendrite).normCoord(:,2), 'LineWidth',1.5, 'color', 'black');
           end
       end
       
    end
end

%% Step 3: Plot the soma in the center
plot(0,0,'ok','MarkerSize',15,'MarkerFaceColor', [1,1,1]);
hold on

%% Step 4: Plot the spines and their property

for i = 1:length(inputs)
    % get the new position relative to the centered soma
    newPos = [allSpines(inputs(i)).xPosNorm, allSpines(inputs(i)).yPosNorm];
    switch property
        case 'deltaOri'
            fieldContent = allSpines(inputs(i)).funcData.deltaOri;
            LUT = inferno(9000);
            plot(newPos(1), newPos(2),'ok','MarkerSize',5,'MarkerFaceColor', LUT((1+floor(fieldContent*100)),:))
        case 'deltaDir'
            fieldContent = allSpines(inputs(i)).funcData.deltaDir;
            LUT = inferno(36000);
            plot(newPos(1), newPos(2),'ok','MarkerSize',5,'MarkerFaceColor', LUT((1+floor(fieldContent*100)),:))
        case 'OSI'
            fieldContent = allSpines(inputs(i)).funcData.OSI;
            LUT = inferno(100);
            plot(newPos(1), newPos(2),'ok','MarkerSize',5,'MarkerFaceColor', LUT((1+floor(fieldContent*100)),:))
        case 'DSI'
            fieldContent = allSpines(inputs(i)).funcData.DSI;
            LUT = viridis(100);
            plot(newPos(1), newPos(2),'ok','MarkerSize',5,'MarkerFaceColor', LUT((1+floor(fieldContent*100)),:))
        case 'struct'
            switch allSpines(inputs(i)).type
                case 'basal'
                    plot(newPos(1), newPos(2), 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'green');
                case 'apical'
                    plot(newPos(1), newPos(2), 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'red');
            end
        case 'CellNr'
           cellNr = allSpines(inputs(i)).CellNr;
           plot(newPos(1), newPos(2),'ok','MarkerSize',5,'MarkerFaceColor', coc(cellNr,:)) 
    end
            
end

%% Step 5: Add colorbar as legend and title
axis off %remove axis
switch property
    case 'struct'
        title('Anatomical distribution', 'FontSize', 14)
    case 'CellNr'
        title('Cell Nr', 'Fontsize', 14)
    otherwise      
        colormap(LUT) 
        colorbar %get the legend
        caxis([0 size(LUT,1)/100]) %define the axis in colorrange
        title(property, 'FontSize', 14) %add the field as the title
end
set(gcf, 'color', 'w'); %set background to white
