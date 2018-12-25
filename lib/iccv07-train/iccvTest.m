% iccvTest

DO_SAVE = 0;

outdir = './results/continuity1/';
%testind = [clusterind(21:23)];
%testind = trainind(1:5);
testind = ftrain(41:46);  %cv_images(21:25);

for f = testind    
    
    disp(num2str(find(f==testind)))
    tmp = load([gtdir strtok(fn{f}, '.') '_gt.mat']);
    bndinfo(f) = orderfields(tmp.bndinfo);                              
   
    
    % compute features
    if ~exist('X', 'var') || numel(X) < f || isempty(X(f).edge)    
        %[X(f), Y(f)] = getAllBoundaryFeatures(bndinfoTest(f), pbdir, gcdir);
        X(f) = getAllFeatures(bndinfo(f), imdir, pbdir, gcdir);
    end
    
    % compute contour likelihoods
    if ~exist('pB', 'var') || numel(pB) < f || isempty(pB{f})
        %pB{f} = useBoundaryClassifier(X{f}(1:end/2, :), dtBnd, 'c');
        pB{f} = useBoundaryClassifier(bndinfo2(f), X2(f), dtBnd4, 'c');
        pB{f} = [pB{f}(:, [1 2]) ; pB{f}(:, [1 3])];

    end        
% 
%     % compute continuity likelihood    
%     if 0 && (~exist('pC', 'var') || numel(pC) < f || isempty(pC{f}))
%         if ~isempty(bndinfo(f).imname)
%             [pC{f}, adjlist{f}] = ...
%                 getContinuityLikelihood(X(f), bndinfo(f), ...
%                 (1:bndinfo(f).ne*2), pB{f}, [], dtCont);
%         else
%             [pC{f}, adjlist{f}] = ...
%                 getContinuityLikelihood(X(f), bndinfoTest(f), ...
%                 (1:bndinfoTest(f).ne*2), pB{f}, [], dtCont);
%         end
%     end
    
    [clab{f}, approxEnergy(f), result(f)] = ...
        occlusionInferenceMerging2(pB{f}, bndinfo2(f), 3.5, 0);    
    
    figure(4), hold off, plot(valdata(f).mergeCost, valdata(f).segError);
    
    %[factors, vars] = getContourPotentials(pB{f}, pC{f}, adjlist{f}, bndinfoTest(f));
    %energy(f) = computeEnergy(vars, factors, clab{f});
    
    imname = bndinfo(f).imname;   
    im = imread([imdir imname]);      
        
    %[pcim, pcim2, pim] = displayContinuityLikelihoods(dtBnd, dtCont, bndinfoTest(f), X{f}, Y{f});
    
    figure(2), imagesc(im), axis image
    
    pcim = zeros(size(bndinfo2(f).wseg));   
    for k = 1:size(pB{f},1)/2
        pcim(bndinfo2(f).edges.indices{k}) = 1-pB{f}(k, 1);
    end
    figure(3), imagesc(ordfilt2(pcim,9,ones(3))), axis image, colormap gray
    
    
%     lab = zeros([size(pContour{f}, 1) 1]);
%     plab = zeros(size(lab));
%     for k = 1:numel(lab)
%         if pContour{f}(k, 1) < 0.75
%             [plab(k), lab(k)] = max(pContour{f}(k, 2:end));
%             %lab(k) = lab(k) + 1;            
%         end
%     end
%     
%     % only show most likely direction
%     ne = bndinfoTest(f).ne;
%     ind = find((lab(1:ne)>0) & (lab(ne+1:end)>0));
%     [minval, minind] = min([plab(ind) plab(ne+ind)], [], 2);
%     lab(ind + (minind-1)*ne) = 0;
%     
%     cim = drawBoundaries(im, bndinfoTest(f), lab);
    
    %[pbim, pcim, ptim] = displayContourResults(im, bndinfoTest(f), ...
    %    Xtest{f}(:, 21), pContour{f}, pType{f});

    if ~isempty(outdir) && 0
        bn = strtok(imname, '.');
        tmp = load([gcdir bn '.c.mat']);    
        gconf = tmp.cim_geom_7;    
        gconf = cat(3, sum(gconf(:, :, 2:6), 3), gconf(:, :, 1), gconf(:, :, 7));

        imwrite(im, [outdir imname]);
        imwrite(pim, [outdir bn '_u.jpg']);
        imwrite(pcim, [outdir bn '_c1.jpg']);
        imwrite(pcim2, [outdir bn '_c2.jpg']);
        
        %imwrite(pbim, [outdir strtok(imname, '.') '_pb.jpg'], 'Quality', 95);
        %imwrite(cim, [outdir strtok(imname, '.') '_c.jpg'], 'Quality', 95);
        %imwrite(pcim, [outdir strtok(imname, '.') '_pc.jpg'], 'Quality', 95);
        %imwrite(ptim, [outdir strtok(imname, '.') '_pt.jpg'], 'Quality', 95);
        imwrite(gconf, [outdir strtok(imname, '.'), '_pg.jpg'], 'Quality', 90);
    end

end 
if DO_SAVE
    save('./data/testBndinfo.mat', 'bndinfoTest');
    save('./data/testData.mat', 'Xtest', 'pContour');
end





