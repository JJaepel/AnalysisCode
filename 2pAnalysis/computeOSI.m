function [osi] = computeOSI(numStims, responses)
    %numStims = number of stimulations
    %responses is a 1xnumStims array containing the average response

    [Rpref, preferredDirInd]= max(responses);
    orthIdx=[ (mod(numStims/4 + preferredDirInd -1, numStims) +1), (mod( preferredDirInd -1- numStims/4 , numStims) +1)];

    Rorth= mean(responses(orthIdx));
    osi= (Rpref-Rorth)/(Rpref+Rorth);    
    
    %restrict osi [-1,1]
    osi(osi>1)= 1;
    osi(osi<-1) = -1;
    
end
