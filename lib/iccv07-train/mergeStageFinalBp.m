function [bndinfo2, lab, plab] = mergeStageFinalBp(bndinfo, dtBnd, dtCont, bias, ...
    imdir, pbdir, gcdir)
% [result, dispim] = mergeStage(X, bndinfo, dtBnd, thresh)
% bias sets the bias for turning an edge off (e.g., bias = 60% means that
% we want to be 60% confident that the edge is off)

global DO_DISPLAY;

thresh = 3.0;

for f = 1:numel(bndinfo)

    if 0 || DO_DISPLAY
        disp(num2str(f))                         
        im = im2double(imread(['./iccvGroundTruth/images/' bndinfo(f).imname]));        
        %gtim = drawBoundaries(im, bndinfo(f), bndinfo(f).edges.boundaryType);
        %figure(2), imagesc(gtim), hold off, axis image        
        figure(2), imagesc(im), axis image
    end                               
    
    while 1    
                
        X = getAllFeatures(bndinfo(f), imdir, pbdir, gcdir);
        
        initnsp = bndinfo(f).nseg;  
        [pB, bndx] = useBoundaryClassifier(bndinfo(f), X(f), dtBnd);
        pB = [pB(:, [1 2]) ; pB(:, [1 3])];
                
        
        [pC, adjlist] = getContinuityLikelihood2(X(f), bndinfo(f), ...
                 (1:bndinfo(f).ne*2), pB, dtCont, bndx);                                                            
             
        [factors, f2var] = getContourPotentials(pB, pC, adjlist, bndinfo(f));
        nnodes = 3*ones(bndinfo(f).ne, 1);
        [mllab, bel] = maxBeliefPropBethe(factors, f2var, nnodes, 0.01, 0.15, 50);
        bel = cell2num(bel')';          
        tmpbel = [bel(:, 1:2) ; bel(:, [1 3])];  
    
        [tmp1, result] = mergeMin(tmpbel, bndinfo(f), -thresh, 2, 0);            
        
        finalnsp = numel(result.regions);
        disp(['Regions: ' num2str(initnsp) ' --> ' num2str(finalnsp)]);            
    
        if 0 || DO_DISPLAY
            pcim = zeros(size(bndinfo(f).wseg));   
            for k = 1:size(pB,1)/2
                pcim(bndinfo(f).edges.indices{k}) = 1-bel(k,1);% 1-pB(k, 1);
            end
            figure(3), hold off, imagesc(ordfilt2(pcim,9,ones(3))), axis image, colormap gray    
        end   
               
        drawnow;
        
        if finalnsp==initnsp
            lab{f} = [(bel(:, 2)>bel(:, 3)) ; (bel(:, 3) >= bel(:, 2))];
            plab{f} = bel;                        
            break;
        else
            bndinfo2(f) = updateBoundaryInfo2(bndinfo(f), result(f));
            bndinfo(f) = bndinfo2(f);
        end
    end
        
    
end