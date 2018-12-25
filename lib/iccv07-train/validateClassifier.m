function [thresh, valdata] = validateClassifier(X, bndinfo, dtBnd, errThresh)
% [thresh, valdata] = validateClassifier(X, bndinfo, dtBnd, errThresh)

global DO_DISPLAY;

for f = 1:numel(X)

    disp(num2str(f)) 
    if DO_DISPLAY                                
        im = im2double(imread(['./iccvGroundTruth/images/' bndinfo(f).imname]));        
        gtim = drawBoundaries(im, bndinfo(f), bndinfo(f).edges.boundaryType);
        figure(2), imagesc(gtim), hold off, axis image        
    end   
    
    pB = useBoundaryClassifier(bndinfo(f), X(f), dtBnd, 'c');
    pB = [pB(:, [1 2]) ; pB(:, [1 3])];

    [tmp1, tmp2, valdata(f)] = mergeMin(pB, bndinfo(f), 0, 2, 1);           
    
    pcim = zeros(size(bndinfo(f).wseg));   
    for k = 1:size(pB,1)/2
        pcim(bndinfo(f).edges.indices{k}) = 1-pB(k, 1);
    end
    figure(3), imagesc(ordfilt2(pcim,9,ones(3))), axis image, colormap gray
   
    drawnow;
    
end

cthresh = (-10:0.1:0);
err = zeros(size(cthresh));
nregions = zeros(size(cthresh));
for f = 1:numel(valdata)
    for k = 1:numel(cthresh)
        [tmp, minind] = min(abs(valdata(f).mergeCost-cthresh(k)));
        err(k) = err(k) + valdata(f).segError(minind)/numel(valdata);
        nregions(k) = nregions(k) + valdata(f).nregions(minind)/numel(valdata);
    end
end
figure(1), hold off, plot(cthresh, err)
figure(2), hold off, plot(err, nregions)

rthresh = (25:25:1000);
err = zeros(size(rthresh));
nregions = zeros(size(rthresh));
for f = 1:numel(valdata)
    for k = 1:numel(rthresh)
        [tmp, minind] = min(abs(valdata(f).nregions-rthresh(k)));
        err(k) = err(k) + valdata(f).segError(minind)/numel(valdata);
        nregions(k) = nregions(k) + valdata(f).nregions(minind)/numel(valdata);
    end
end
figure(3), hold off, plot(rthresh, err)
figure(4), hold off, plot(err, nregions)



[tmp, minind] = min(abs(err-errThresh));

thresh = -cthresh(minind);

tmp = valdata;
clear valdata;
valdata.all = tmp;
valdata.thresh = cthresh;
valdata.err = err;

