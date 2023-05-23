function onsetState = computeOnsetState(runningBlocks,mod)
% For each mouse, this function computes the average of onset correlation
% maps across all legit running blocks. This returns a 92 x 92 matrix, 
% which is considered as the onset state.
% 
% Author: Somayeh "Bahar" Shahsavarani
% email: bahar@huskers.unl.edu

% convert struct to table
runningBlocksTable = struct2table(runningBlocks);


for n = 1:height(runningBlocksTable)
    
    if strcmp(mod,'jRGECO')
        C = runningBlocksTable.correlations_neural(n);
        FCmaps(:,:,n) = C{1};
        
    elseif strcmp(mod,'HbT')
        C = runningBlocksTable.correlations_HbT(n);
        FCmaps(:,:,n) = C{1};
    end
    
end

% mean of FCmaps
onsetState = mean(FCmaps,3);

end