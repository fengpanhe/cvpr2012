function testMultipleSegs(fn, segdir, outdir, gtdir, imdir, pbdir, gcdir, ...
    dtBnd, nsegs, nregions)

for f = 1:numel(fn)
    
    disp(num2str(f));
    
    bn = strtok(fn{f}, '.');
    
    load([segdir bn '_seg.mat']);
    gt = load([gtdir bn '_gt.mat']);
    
    X = getAllFeatures(bndinfo, imdir, pbdir, gcdir);

    pB = useBoundaryClassifier(bndinfo, X, dtBnd);
    pB = [pB(:, [1 2]) ; pB(:, [1 3])];        
    
    for k = 1:nsegs

        tmppb = pB;
        if k > 1
            tmppb = tmppb + rand(size(pB))*0.05;
            tmppb = tmppb ./ repmat(sum(tmppb, 2), [1 size(pB, 2)]);
        end
                
        [lab, result] = mergeMin(tmppb, bndinfo, 0, nregions, 0);                                   
        bndinfo2(k) = updateBoundaryInfo2(bndinfo, result);     
              
        energy(k) = sum(-log(pB(lab==0, 1))) + sum(-log(1-pB(lab>0, 1)));
        
        tmp = evaluateSegmentationError(gt.bndinfo, bndinfo2(k).wseg); 
        pixacc(k) = tmp.conservation;        
        
        disp(['  ' num2str(k) ': pixacc = ' num2str(pixacc(k)) '  energy = ' ...
            num2str(energy(k))]);
        
        if exist([outdir 'display'], 'file')
            imwrite(label2rgb(bndinfo2(k).wseg), ...
                [outdir 'display/' bn '_seg' num2str(k) '.jpg'], 'Quality', 99);
        end
        
    end
    
    bndinfo = bndinfo2;
    clear bndinfo2;
    save([outdir bn '_seg.mat'], 'bndinfo', 'energy', 'pixacc');    

end    

