function [bndinfo2, lab, plab_e, plab_g] = mergeStageFinalBpGeometry(...
    bndinfo, dtBnd, dtBnd_fast, dtCont, thresh, imdir, pbdir, gcdir, gdatadir, cvnum, ...
    dorestimate)
% [result, dispim] = mergeStage(X, bndinfo, dtBnd, thresh)
% bias sets the bias for turning an edge off (e.g., bias = 60% means that
% we want to be 60% confident that the edge is off)

global DO_DISPLAY;

segdir = './iccvGroundTruth/segs/';

for f = 1:numel(bndinfo)

    if 0 || DO_DISPLAY
        disp(num2str(f))                         
        im = im2double(imread(['./iccvGroundTruth/images/' bndinfo(f).imname]));        
        %gtim = drawBoundaries(im, bndinfo(f), bndinfo(f).edges.boundaryType);
        %figure(2), imagesc(gtim), hold off, axis image        
        figure(2), imagesc(im), axis image
    end                               
    
    while 1    
        if dorestimate
            segmaps = bndinfo2segmaps(bndinfo(f), segdir);
            X = getAllFeatures(bndinfo(f), imdir, pbdir, gcdir, gdatadir, segmaps, cvnum);                
        else
            X = getAllFeatures(bndinfo(f), imdir, pbdir, gcdir);
        end
        
        initnsp = bndinfo(f).nseg;  
        [pB, bndx] = useBoundaryClassifier(bndinfo(f), X(f), dtBnd);
        pB = [pB(:, [1 2]) ; pB(:, [1 3])];
                        
        [pC, adjlist] = getContinuityLikelihood2(X(f), bndinfo(f), ...
                 (1:bndinfo(f).ne*2), pB, dtCont, bndx);                                                            
             
        pg = X(f).region.geomContext;
        [factors, f2var] = getContourGeometryPotentials3(pB, pC, adjlist, pg, bndinfo(f));
        nnodes = [3*ones(bndinfo(f).ne, 1) ; 5*ones(bndinfo(f).nseg, 1)];
        [mllab, bel] = maxBeliefPropBethe(factors, f2var, nnodes, 0.025, 0.5, Inf);
        ebel = cell2mat(bel(1:bndinfo(f).ne)')';          
        gbel = cell2mat(bel(bndinfo(f).ne+1:end)')';
        tmpbel = [ebel(:, 1:2) ; ebel(:, [1 3])];  

        result = mergeMinSafe(tmpbel, bndinfo(f), X(f), dtBnd_fast, thresh, 2, 0);               
        
        finalnsp = numel(result.regions);
        disp(['Regions: ' num2str(initnsp) ' --> ' num2str(finalnsp)]);            
    
        if 0 || DO_DISPLAY
            pcim = zeros(size(bndinfo(f).wseg));   
            for k = 1:size(pB,1)/2
                pcim(bndinfo(f).edges.indices{k}) = 1-bel(k,1);% 1-pB(k, 1);
            end
            figure(3), hold off, imagesc(ordfilt2(pcim,9,ones(3))), axis image, colormap gray    
            drawnow;            
        end                  
        
        if finalnsp==initnsp
            lab{f} = [(ebel(:, 2)>ebel(:, 3)) ; (ebel(:, 3) >= ebel(:, 2))];
            plab_e{f} = ebel;                        
            plab_g{f} = gbel;
            break;
        else
            bndinfo2(f) = updateBoundaryInfo2(bndinfo(f), result(f));
            bndinfo(f) = bndinfo2(f);
        end
    end
        
    
end