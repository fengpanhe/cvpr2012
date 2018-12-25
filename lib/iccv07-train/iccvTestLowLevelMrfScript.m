% iccvTestLowLevelScript

DO_READ         = 0;
DO_FEATURES     = 0;
DO_TRAIN        = 0;
DO_TEST         = 0;
DO_TEST_GEOM    = 1;

testind = cv_images(1:50);
trainind = setdiff(clusterind, 230);


global DO_DISPLAY;
DO_DISPLAY = 0;

%% Read train data and compute features

if DO_READ
    disp('Reading...');
    for f = trainind 
        tmp = load([gtdir strtok(fn{f}, '.') '_gt.mat']);
        f2 = find(trainind==f);
        bndinfo(f2) = tmp.bndinfo;
    end
end

if DO_FEATURES
    disp('Computing features...')
    [X, Y] = getAllFeatures(bndinfo, imdir, pbdir, gcdir);
end


%% Train mrf classifier
if DO_TRAIN

    disp('Training continuity classifier...');   
    for f = 1:numel(X)       
        [pB{f}, trainx{f}] = useBoundaryClassifier(bndinfo(f), X(f), dtBnd_all);
        pB{f} = [pB{f}(:, [1 2]) ; pB{f}(:, [1 3])];
    end
    dtCont = trainBoundaryContinuityClassifier(bndinfo, X, Y, pB, trainx);    

end


%% Test mrf classifier

if DO_TEST        
    
    disp('Getting test likelihoods ...')
    for f = testind         
        
        tmp = load([gtdir strtok(fn{f}, '.') '_gt.mat']);
        bndinfo = tmp.bndinfo;
        f2 = find(f==testind);
        disp(num2str(f2))

        if 0 % no geometry update
            [X, Y] = getAllFeatures(bndinfo, imdir, pbdir, gcdir);
        else % update geometry
            segmaps = bndinfo2segmaps(bndinfo, segdir);
            [X, Y] = getAllFeatures(bndinfo, imdir, pbdir, gcdir, gdatadir, segmaps, 1);
        end          
        Y = Y{1};        
        w = X.edge.length; w = w / sum(w);
        
        figure(1), hold off, plot(result_all(f2).pr.r, result_all(f2).pr.p, '--b');

        [pB, bndx] = useBoundaryClassifier(bndinfo, X, dtBnd_all);                
        pB = [pB(:, [1 2]) ; pB(:, [1 3])];
        [pC, adjlist] = getContinuityLikelihood2(X, bndinfo, ...
                 (1:bndinfo.ne*2), pB, dtCont(end), bndx);                                      
      
    
        [factors, f2var] = getContourPotentials(pB, pC, adjlist, bndinfo);
        nnodes = 3*ones(bndinfo.ne, 1);
        [mllab, bel] = maxBeliefPropBethe(factors, f2var, nnodes, 0.025, 0.5, Inf);
        pBmrf = cell2mat(bel')';                         
        
        tmp = evaluateLowLevelPerformance(Y, pBmrf, w);
        tmp.imname = fn{f};
        tmp.imnum = f;        
        tmp.pB = pBmrf;
        tmp.labels = Y;        
        result_mrf(f2) = tmp;
        figure(1), hold on, plot(result_mrf(f2).pr.r, result_mrf(f2).pr.p, 'g');                
        drawnow;
        
        disp(['Figure/Ground accuracy: ' num2str([result_all(f2).fgaccuracy ...
            result_mrf(f2).fgaccuracy])]);
        
    end
    
    [pr_mrf, roc_mrf, fgacc_mrf] = summarizeLowLevelPerformance(result_mrf);
    
    save('./results/classifierComparisonResult_stage3mrf.mat', ...
        'dtBnd_all', 'dtCont', 'result_mrf', 'pr_mrf', 'roc_mrf', 'fgacc_mrf');        
    
end
   
if DO_TEST_GEOM        
    
    disp('Getting test likelihoods ...')
    for f = testind         
        
        tmp = load([gtdir strtok(fn{f}, '.') '_gt.mat']);
        bndinfo = tmp.bndinfo;
        f2 = find(f==testind);
        disp(num2str(f2))

        if 0 % no geometry update
            [X, Y] = getAllFeatures(bndinfo, imdir, pbdir, gcdir);
        else % update geometry
            segmaps = bndinfo2segmaps(bndinfo, segdir);
            [X, Y] = getAllFeatures(bndinfo, imdir, pbdir, gcdir, gdatadir, segmaps, 1);
        end          
        Y = Y{1};        
        w = X.edge.length; w = w / sum(w);
        
        figure(1), hold off, plot(result_mrfg3(f2).pr.r, result_mrfg3(f2).pr.p, '--b');
        axis([0 1 0 1])
        
        [pB, bndx] = useBoundaryClassifier(bndinfo, X, dtBnd_all);                
        pB = [pB(:, [1 2]) ; pB(:, [1 3])];
        [pC, adjlist] = getContinuityLikelihood2(X, bndinfo, ...
                 (1:bndinfo.ne*2), pB, dtCont(end), bndx);                                      
      
        %pg = [X.region.geomContext(:, 1) sum(X.region.geomContext(:, 2:4), 2)  X.region.geomContext(:, 5)];
        pg = X.region.geomContext;
        [factors, f2var] = getContourGeometryPotentials3(pB, pC, adjlist, pg, bndinfo);
        nnodes = [3*ones(bndinfo.ne, 1) ; 5*ones(bndinfo.nseg, 1)];
        [mllab, bel] = maxBeliefPropBethe(factors, f2var, nnodes, 0.025, 0.5, Inf);                
        pBmrf = cell2mat(bel(1:bndinfo.ne)')';                         
        
        tmp = evaluateLowLevelPerformance(Y, pBmrf(1:bndinfo.ne, :), w);
        tmp.imname = fn{f};
        tmp.imnum = f;        
        tmp.pB = bel;
        tmp.labels = Y;        
        result_mrfg4(f2) = tmp;
        figure(1), hold on, plot(result_mrfg4(f2).pr.r, result_mrfg4(f2).pr.p, 'g');                
        drawnow;
        
        disp(['Figure/Ground accuracy: ' num2str([result_mrfg3(f2).fgaccuracy ...
            result_mrfg4(f2).fgaccuracy])]);
        
    end
    
    [pr_mrfg4, roc_mrfg4, fgacc_mrfg4] = summarizeLowLevelPerformance(result_mrfg4);
    
    save('./results/classifierComparisonResult_stage3mrfg4.mat', ...
        'dtBnd_all', 'dtCont', 'result_mrfg4', 'pr_mrfg4', 'roc_mrfg4', 'fgacc_mrfg4');        
    
end

%% Clean-up
clear DO_READ DO_FEATURES DO_TRAIN DO_TEST