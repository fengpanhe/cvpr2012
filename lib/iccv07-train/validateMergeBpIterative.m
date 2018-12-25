function valdata = validateMergeBpIterative(bndinfo, dtBnd, dtBnd_fast, ...
    dtCont, imdir, pbdir, gcdir, gtdir)
% [thresh, valdata] = validateClassifier(X, bndinfo, dtBnd, errThresh)

for f = 1:numel(bndinfo)

    disp(num2str(f))          
    
    cthresh = (0.10:0.05:0.30);
    for c = 1:numel(cthresh)
        
        bndinfo2 = bndinfo(f);
    
        for iter = 1:3
    
            initRegions = bndinfo2.nseg;
            X = getAllFeatures(bndinfo2, imdir, pbdir, gcdir);
    
            [pB, bndx] = useBoundaryClassifier(bndinfo2, X, dtBnd);                
            pB = [pB(:, [1 2]) ; pB(:, [1 3])];
            [pC, adjlist] = getContinuityLikelihood2(X, bndinfo2, ...
                 (1:bndinfo2.ne*2), pB, dtCont, bndx);        
    
            [factors, f2var] = getContourPotentials(pB, pC, adjlist, bndinfo2);
            nnodes = 3*ones(bndinfo(f).ne, 1);
            [mllab, bel] = maxBeliefPropBethe(factors, f2var, nnodes, 0.025, 0.5, Inf);
            bel = cell2mat(bel')';  
            tmpbel = [bel(:, 1:2) ; bel(:, [1 3])];                    
            
            %[tmp1, tmp2, valdata(f)] = mergeMin(tmpbel, bndinfo(f), 0, 2, 1);           
            result = mergeMinSafe(tmpbel, bndinfo2, X, dtBnd_fast, cthresh(c), 2, 0);  
            
            bndinfo2 = updateBoundaryInfo2(bndinfo2, result);
            
            finalRegions = bndinfo2.nseg;
            disp(['Iter ' num2str(iter) ' init ' num2str(initRegions) ' --> ' ' final ' num2str(finalRegions)]);
            
        end
        
        gt = load([gtdir strtok(bndinfo(f).imname, '.') '_gt.mat']);
        tmp = evaluateSegmentationError(gt.bndinfo, bndinfo2.wseg);
        tmp.cthresh = cthresh(c);
        tmp.niter = iter;
        valdata(f, c) = tmp;                        
        disp(num2str([tmp.conservation tmp.efficiency]));
    end
    
    plotIterativeValdata(valdata, 3);
            
end



