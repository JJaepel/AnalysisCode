function [dsi]= computeDSI(numStims, responses, type)
    %numStims = number of directions
    %responses is a 1xnumStims array containing the average response
    if nargin < 3
        type = 'preferred';
    end
    
    
    [~, prefIdx]= max(responses);
    nullIdx = mod(numStims/2 + prefIdx -1, numStims)+1; %this wraps our indices from 1-N to find the null directio
    
    switch type
        case 'preferred'
            Rnull = responses(nullIdx);
            Rpref = responses(prefIdx);

        case 'hemifield'
            
            hemifield = -numStims/4 : numStims/4;
            
            nullHemi = mod(nullIdx + hemifield -1, numStims)+1;
            prefHemi = mod(prefIdx + hemifield -1, numStims)+1;
            
                       
            Rpref = sum(responses(prefHemi));
            Rnull = sum(responses(nullHemi));
    end
    
    dsi= (Rpref-Rnull)/(Rpref+Rnull);
    
    %restrict dsi [-1,1]
    dsi(dsi>1)= 1;
    dsi(dsi<-1) = -1;
    
end
