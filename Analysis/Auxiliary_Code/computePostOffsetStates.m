function postRunningStates = computePostOffsetStates(runningBlocks,windows,mod)
% For each mouse, this function computes the average of resting-state correlation
% maps across all legit running blocks. It returns 2 matrices, which are
% initial- and sustained-rest states, respectively.
% 
% Author: Somayeh "Bahar" Shahsavarani
% email: bahar@huskers.unl.edu

% convert struct to table
runningBlocksTable = struct2table(runningBlocks);


for w  = 1:length(windows)
    
    for n = 1:height(runningBlocksTable)
        
        if strcmp(mod,'jRGECO')
            C = runningBlocksTable.correlations_neural(n);
            C = C{1};
      
        elseif strcmp(mod,'HbT')
            C = runningBlocksTable.correlations_HbT(n);
            C = C{1};
            
        end
        
        mytemp(:,:,n) = C(windows(w),:,:);
        
    end
    FCmaps{w} = mytemp;
    
    
    clear C mytemp myts
end

for w = 1:length(windows)
    mytemp = FCmaps{w};
    myvar = mean(mytemp,3);
    meanvar(:,:,w) = myvar;
    %stdvar(:,:,w) = std(mytemp,[],3);
    
    clear mytemp myvar myts myvar_mean myvar_se
end
postRunningStates = meanvar;
%stdCorrelations = stdvar;

end