% iccvTrainMRF

%% Parameters

DO_LOAD         = 0;
DO_MRF_DATA     = 0;
DO_MRF_MERGE    = 0;
DO_MRF_MERGE_FINAL    = 1;

%% Processes

if DO_LOAD
    load([datadir 'bndClassifiersAll.mat']);
    load([datadir 'contClassifiersAll.mat']);
end

if DO_MRF_DATA
    
    for f = ftrain
        disp(num2str(find(f==ftrain)))        
        [X, Y(f)] = getAllFeatures(bndinfo(f), imdir, pbdir, gcdir);
        
        [pB{f}, bndx] = useBoundaryClassifier(bndinfo(f), X, dtBnd(end));
        pB{f} = [pB{f}(:, [1 2]) ; pB{f}(:, [1 3])];

        [pC{f}, adjlist{f}] = getContinuityLikelihood2(X, bndinfo(f), ...
                 (1:bndinfo(f).ne*2), pB{f}, dtCont(end), bndx);         
    end
    save('./data/initMrfData2.mat', 'Y', 'pB', 'pC', 'adjlist');
end

if DO_MRF_MERGE
         
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
end


if DO_MRF_MERGE_FINAL
         
    for f = testind(15:end)
        disp(num2str(find(f==testind)))
        if 0 && ismember(f, ftrain)
            tmp = load([gtdir strtok(fn{f}, '.') '_gt.mat']);
            bndinfo = tmp.bndinfo;
        else
            tmp = load([smallsegdir '/2/' strtok(fn{f}, '.') '_seg.mat']);
            bndinfo = tmp.bndinfo;            
        end

        [bndinfo, lab] = mergeStageFinalBp(bndinfo, dtBnd(end), dtCont(end), 0.70, imdir, pbdir, gcdir);  
        lab = lab{1};               
        
        bndinfo.edges.boundaryType = lab;
        %bndinfo2 = transferSuperpixelLabels(bndinfo2, bndinfo.wseg);
        
        
        save([smallsegdir '/3_2/' strtok(fn{f}, '.') '_seg.mat'], 'bndinfo');    
        
        im = im2double(imread([imdir bndinfo.imname]));        
        
        gtim = drawBoundaries(im, bndinfo, tmplab);
        %figure(2), imagesc(gtim), hold off, axis image                                
        outdir = [smallsegdir '/3_2/display/'];
        imwrite(im, [outdir fn{f}]);
        printOcclusionResult(im, bndinfo, lab, [outdir strtok(fn{f}, '.') '_res.jpg'], 1);
        %imwrite(gtim, [outdir strtok(fn{f}, '.') '_res.jpg'], 'Quality', 99);
        %imwrite(label2rgb(bndinfo.wseg), [smallsegdir '/3/display/' strtok(fn{f}, '.') '_seg.jpg'], 'Quality', 99);
    end            
end
