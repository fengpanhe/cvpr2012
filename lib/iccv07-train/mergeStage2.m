function [result, dispim] = mergeStage2(X, bndinfo, dtBnd, dtCont, thresh)
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

    [pC, adjlist] = getContinuityLikelihood2(X(f), bndinfo(f), ...
                 (1:bndinfo(f).ne*2), pB, dtCont, bndx);                         
             
    initnsp = bndinfo(f).nseg;
    
    [lab, result(f)] = ...
            mergeEnergyBased(pB, pC, adjlist, bndinfo(f), thresh);    
    
%     [tmp1, tmp2, result(f), tmp3, dispim] = ...
%         occlusionInferenceMerging2(pB, bndinfo(f), thresh, 0);               
    
    finalnsp = numel(result(f).regions);
    disp(['Regions: ' num2str(initnsp) ' --> ' num2str(finalnsp)]);
    
    
    if DO_DISPLAY
        pcim = zeros(size(bndinfo(f).wseg));   
        for k = 1:size(pB,1)/2
            pcim(bndinfo(f).edges.indices{k}) = 1-pB(k, 1);
        end
        figure(3), imagesc(ordfilt2(pcim,9,ones(3))), axis image, colormap gray    
    end   
    
    keyboard;
end


