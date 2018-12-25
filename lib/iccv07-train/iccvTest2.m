


for f = testind
    
    disp(num2str(find(f==testind)))        
    
    
    
    [X, Y(f)] = getAllFeatures(bndinfo(f), imdir, pbdir, gcdir);

    [pB{f}, bndx] = useBoundaryClassifier(bndinfo(f), X, dtBnd(end));
    pB{f} = [pB{f}(:, [1 2]) ; pB{f}(:, [1 3])];

    [pC{f}, adjlist{f}] = getContinuityLikelihood2(X, bndinfo(f), ...
             (1:bndinfo(f).ne*2), pB{f}, dtCont(end), bndx);        
             
             
             
    for f = ftrain
        disp(num2str(f))
        if 0 && ismember(f, ftrain)
            tmp = load([gtdir strtok(fn{f}, '.') '_gt.mat']);
            bndinfo = tmp.bndinfo;
        else
            tmp = load([smallsegdir '/2/' strtok(fn{f}, '.') '_seg.mat']);
            bndinfo = tmp.bndinfo;            
        end

        X = getAllFeatures(bndinfo, imdir, pbdir, gcdir);
        thresh = 2.6;
        result = mergeStageBp2(X, bndinfo, dtBnd(end), dtCont(end), thresh);                  
        bndinfo = updateBoundaryInfo2(bndinfo, result);    
        
        %save([smallsegdir '2/' strtok(fn{f}, '.') '_seg.mat'], 'bndinfo');    
    end   
    
    
    
    
    for f = testind(32:end)
        disp(num2str(find(f==testind)))
        if 0 && ismember(f, ftrain)
            tmp = load([gtdir strtok(fn{f}, '.') '_gt.mat']);
            bndinfo = tmp.bndinfo;
        else
            tmp = load([smallsegdir '/2/' strtok(fn{f}, '.') '_seg.mat']);
            bndinfo = tmp.bndinfo;            
        end

        [bndinfo, lab] = mergeStageFinalBp(bndinfo, dtBnd(end), dtCont(end), 1, imdir, pbdir, gcdir);  
        lab = lab{1};      
        tmplab = lab*2; 
        
        %bndinfo2.edges.boundaryType = tmplab;
        %bndinfo2 = transferSuperpixelLabels(bndinfo2, bndinfo.wseg);
        
        
        save([smallsegdir '/3/' strtok(fn{f}, '.') '_seg.mat'], 'bndinfo');    
        
        im = im2double(imread([imdir bndinfo.imname]));        
        
        gtim = drawBoundaries(im, bndinfo, tmplab);
        %figure(2), imagesc(gtim), hold off, axis image                                
        imwrite(im, [smallsegdir '/3/display/' fn{f}]);
        imwrite(gtim, [smallsegdir '/3/display/' strtok(fn{f}, '.') '_res.jpg'], 'Quality', 99);
        imwrite(label2rgb(bndinfo.wseg), [smallsegdir '/3/display/' strtok(fn{f}, '.') '_seg.jpg'], 'Quality', 99);
    end       