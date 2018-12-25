function [thresh, valdata] = validateMergeBpGeometry(bndinfo, X, dtBnd, dtBnd_fast, dtCont, errThresh)
% [thresh, valdata] = validateClassifier(X, bndinfo, dtBnd, errThresh)

global DO_DISPLAY;

for f = 1:numel(bndinfo)

    disp(num2str(f)) 
    if DO_DISPLAY                                
        im = im2double(imread(['./iccvGroundTruth/images/' bndinfo(f).imname]));        
        figure(2), hold off, imagesc(im), axis image
        %gtim = drawBoundaries(im, bndinfo(f), bndinfo(f).edges.boundaryType);
        %figure(2), imagesc(gtim), hold off, axis image        
    end   
    
    [pB, bndx] = useBoundaryClassifier(bndinfo(f), X(f), dtBnd);                
    pB = [pB(:, [1 2]) ; pB(:, [1 3])];
    [pC, adjlist] = getContinuityLikelihood2(X(f), bndinfo(f), ...
             (1:bndinfo(f).ne*2), pB, dtCont, bndx);        
    
    pg = X(f).region.geomContext;     
    [factors, f2var] = getContourGeometryPotentials3(pB, pC, adjlist, pg, bndinfo(f));
    nnodes = [3*ones(bndinfo(f).ne, 1) ; 5*ones(bndinfo(f).nseg, 1)];
    [mllab, bel] = maxBeliefPropBethe(factors, f2var, nnodes, 0.025, 0.5, Inf);
    bel = cell2mat(bel(1:bndinfo(f).ne)')';
    tmpbel = [bel(:, 1:2) ; bel(:, [1 3])];
    
    %[tmp1, tmp2, valdata(f)] = mergeMin(tmpbel, bndinfo(f), 0, 2, 1);           
    [tmp1, valdata(f)] = mergeMinSafe(tmpbel, bndinfo(f), X(f), dtBnd_fast, 0.5, 2, 1);  
    
    if DO_DISPLAY    
        pcim = zeros(size(bndinfo(f).wseg));   
        for k = 1:size(tmpbel,1)/2
            pcim(bndinfo(f).edges.indices{k}) = 1-tmpbel(k, 1);
        end
        figure(3), imagesc(ordfilt2(pcim, 9, ones(3))), axis image, colormap gray
    end
   
    drawnow;
    
end

cthresh = (0.001:0.001:0.5);
err = zeros(size(cthresh));
nregions = zeros(size(cthresh));
for f = 1:numel(valdata)
    for k = 1:numel(cthresh)
        [tmp, minind] = min(abs(valdata(f).mergeCost-cthresh(k)));
        err(k) = err(k) + valdata(f).segError(minind)/numel(valdata);
        nregions(k) = nregions(k) + valdata(f).nregions(minind)/numel(valdata);
    end
end
figure(1), subplot(2,2,1), hold off, plot(cthresh, err), title('bias vs. err')
figure(1), subplot(2,2,2), plot(err, nregions), title('err vs. nregions')
cerr = err;

rthresh = (5:5:1000);
err = zeros(size(rthresh));
nregions = zeros(size(rthresh));
for f = 1:numel(valdata)
    for k = 1:numel(rthresh)
        [tmp, minind] = min(abs(valdata(f).nregions-rthresh(k)));
        err(k) = err(k) + valdata(f).segError(minind)/numel(valdata);
        nregions(k) = nregions(k) + valdata(f).nregions(minind)/numel(valdata);
    end
end
figure(1), subplot(2,2,3), plot(rthresh, err),  title('minregion vs. err')
figure(1), subplot(2,2,4), plot(err, nregions), title('err vs. nregions')



[tmp, minind] = min(abs(err-errThresh));

thresh = -cthresh(minind);

tmp = valdata;
clear valdata;
valdata.all = tmp;
valdata.thresh = cthresh;
valdata.err = cerr;

