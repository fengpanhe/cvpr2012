function [result, dispim] = mergeStageBp(X, bndinfo, dtBnd, dtCont, thresh)
% [result, dispim] = mergeStage(X, bndinfo, dtBnd, thresh)

global DO_DISPLAY;

for f = 1:numel(X)

    if DO_DISPLAY
        disp(num2str(f))                         
        im = im2double(imread(['./iccvGroundTruth/images/' bndinfo(f).imname]));        
        gtim = drawBoundaries(im, bndinfo(f), bndinfo(f).edges.boundaryType);
        figure(2), imagesc(gtim), hold off, axis image        
    end   
    
    [pB, bndx] = useBoundaryClassifier(bndinfo(f), X(f), dtBnd, 'c');
    pB = [pB(:, [1 2]) ; pB(:, [1 3])];

    pB(:, 1) = pB(:, 1) * exp(-thresh);
    
    [pC, adjlist] = getContinuityLikelihood2(X(f), bndinfo(f), ...
                 (1:bndinfo(f).ne*2), pB, dtCont, bndx);                         
             
    initnsp = bndinfo(f).nseg;
      
    
    [factors, f2var] = getContourPotentials(pB, pC, adjlist, bndinfo(f));
    nnodes = 3*ones(bndinfo(f).ne, 1);
    [mllab, bel] = maxBeliefPropBethe(factors, f2var, nnodes, 1E-2, 0.25, Inf);
    bel = cell2num(bel')';          
    
%     [lab, result(f)] = ...
%            mergeEnergyBased(pB, pC, adjlist, bndinfo(f), thresh);    
    
%     [tmp1, tmp2, result(f), tmp3, dispim] = ...
%         occlusionInferenceMerging2(pB, bndinfo(f), thresh, 0);               
    
%    finalnsp = numel(result(f).regions);
%    disp(['Regions: ' num2str(initnsp) ' --> ' num2str(finalnsp)]);
    
    
    if DO_DISPLAY
        pcim = zeros(size(bndinfo(f).wseg));   
        for k = 1:size(pB,1)/2
            pcim(bndinfo(f).edges.indices{k}) = 1-bel(k,1);% 1-pB(k, 1);
        end
        figure(3), hold off, imagesc(ordfilt2(pcim,9,ones(3))), axis image, colormap gray    
    end   
    
    keyboard;
    
end


