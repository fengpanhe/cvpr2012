function [thresh, valdata] = validateClassifier2(X, bndinfo, dtBnd, dtCont, errThresh)
% [thresh, valdata] = validateClassifier(X, bndinfo, dtBnd, errThresh)

global DO_DISPLAY;

if 0 % for testing
    ind = (4:6);
    disp(['Processing only ' num2str(ind)]);
    X = X(ind);  
    bndinfo = bndinfo(ind);
end

bias = (1.5:-0.1:0);
err = zeros(numel(X), numel(bias));
nregions = zeros(numel(X), numel(bias));

for f = 1:numel(X)

    disp(num2str(f))
    
    if DO_DISPLAY                                 
        im = im2double(imread(['./iccvGroundTruth/images/' bndinfo(f).imname]));        
        gtim = drawBoundaries(im, bndinfo(f), bndinfo(f).edges.boundaryType);
        figure(1), imagesc(gtim), hold off, axis image 
        drawnow;
    end   
    
    [pB, bndx] = useBoundaryClassifier(bndinfo(f), X(f), dtBnd, 'c');
    pB = [pB(:, [1 2]) ; pB(:, [1 3])];

    [pC, adjlist] = getContinuityLikelihood2(X(f), bndinfo(f), ...
                 (1:bndinfo(f).ne*2), pB, dtCont, bndx);    
    
%     [factors, f2var] = getContourPotentials(pB, pC, adjlist, bndinfo(f));
%     nnodes = 3*ones(bndinfo(f).ne, 1);
%     [mllab, bel] = maxBeliefPropBethe(factors, f2var, nnodes, 1E-2, 0.25, Inf);
%     bel2 = cell2num(bel')';  

    stats = regionprops(bndinfo(f).wseg, 'Area');
    areas = cat(1, stats(:).Area);
    nsp = bndinfo.nseg;    
    
    for k = 1:numel(bias)

        [lab, result] = ...
            mergeEnergyBased(pB, pC, adjlist, bndinfo(f), bias(k));             
        
        err(f,k) = getSegmentationErrorFast(bndinfo(f).labels, nsp, result.regions, areas);
        nregions(f,k) = numel(result.regions);           
    
    end
    
   
    
end

errave = mean(err, 1);
nregionsAve = mean(nregions, 1);
figure(3), hold off, plot(bias, errave)
figure(4), hold off, plot(bias, nregionsAve)

[tmp, minind] = min(abs(errave-errThresh) + 1000*(errave>errThresh));

thresh = bias(minind);
valdata.thresh = bias(minind);
valdata.bias = bias;
valdata.err = err;
valdata.nregions = nregions;
valdata.errAve = errave;
valdata.nregionsAve = nregionsAve;

