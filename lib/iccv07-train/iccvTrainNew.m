% iccvTrain


%% Settings

STAGES = (3);

DO_LOAD                  = 0;
DO_READ                  = 0;
DO_RESTART               = 0;
DO_SMALL_SEGS            = 0;
DO_FEATURES              = 0;
DO_FEATURES_GEOMETRY     = 0;
DO_GEOMETRY_PARAMS       = 0;
DO_BOUNDARY_CLASSIFIER   = 0;
DO_CONTINUITY_CLASSIFIER = 0;
DO_VALIDATION            = 0;
DO_MERGE_MIN             = 0;
DO_MERGE_MRF             = 0;
DO_MERGE_MRF_FINAL       = 0;
DO_UNARY_VALIDATION      = 0;
DO_BP_VALIDATION         = 0;
DO_GEOMETRY_ONLY_SEG     = 0;
DO_NCUT_SEG              = 1;

ERR_THRESH = [0.01 ; 0.005 ; 0.005 ; 0.005 ; 0.005 ; 0.005];

ftrain = setdiff(clusterind, 230);
testind = cv_images;
geomcvtrain = cv_images;
geomncv = 5;
NSMALLSEGS = 300;  

thresh = [0.105 0.25];

datadir = '/usr1/projects/dhoiem/iccv07/data4/';

global DO_DISPLAY;
DO_DISPLAY = 0;


%% Training

s = STAGES(1);

if DO_LOAD
    load([datadir '/bndClassifiersAll3.mat']);
    load([datadir '/contClassifiersAll3.mat']);
end

if DO_READ
    disp(['Reading from ' gtdir]);
    for f = ftrain
        tmp = load([gtdir strtok(fn{f}, '.') '_gt.mat']);
        bndinfo(f) = tmp.bndinfo;
    end
    clear tmp;
end

if DO_RESTART
    s = STAGES(1);
    disp(['Restarting from after stage ' num2str(s-1)]);
    %load(['./data/tmp/segmaps' num2str(s-1) '.mat'], 'segmaps'); 
    load([datadir 'trainBndinfo' num2str(s-1) '.mat']);
    bndinfo1 = bndinfo2; 
    clear bndinfo2;    
end

if DO_FEATURES
    disp('Computing features...');    
    [X(ftrain), Y(ftrain)] = getAllFeatures(bndinfo(ftrain), imdir, pbdir, gcdir);
end    

if DO_FEATURES_GEOMETRY
    disp('Computing features with geometry...');    
    segmaps = bndinfo2segmaps(bndinfo(ftrain), segdir);
    [X(ftrain), Y(ftrain)] = getAllFeatures(bndinfo(ftrain), imdir, pbdir, ...
        gcdir, gdatadir, segmaps, 1);
end

if DO_GEOMETRY_PARAMS
    disp('Estimating parameters for geometry regression');
    geomparams = estimateGeometryRegression(bndinfo(ftrain), X(ftrain), gdatadir);
end

if DO_BOUNDARY_CLASSIFIER
    disp('Training boundary...');
    [dtBnd(s), xtrain] = trainBoundaryClassifier4(X(ftrain), Y(ftrain), bndinfo(ftrain));
    dtBnd_fast(s) = trainFastBoundaryClassifier(X(ftrain), Y(ftrain), bndinfo(ftrain), 1);
    save([datadir 'bndClassifiersAllTmp2.mat'], 'dtBnd', 'dtBnd_fast');
end

if DO_CONTINUITY_CLASSIFIER

    disp('Training continuity classifier...');   
    for f = ftrain      
        [pB{f}, trainx{f}] = useBoundaryClassifier(bndinfo(f), X(f), dtBnd(s));
        pB{f} = [pB{f}(:, [1 2]) ; pB{f}(:, [1 3])];
    end
    dtCont(s) = trainBoundaryContinuityClassifier(bndinfo(ftrain), X(ftrain), Y(ftrain), pB(ftrain), trainx(ftrain));    
    save([datadir 'contClassifiersAllTmp2.mat'], 'dtCont');
end

if DO_UNARY_VALIDATION       
    [thresh(s), valdata(s)] = validateMergeMinSafe(X(ftrain), ...
        bndinfo(ftrain), dtBnd(s), dtBnd_fast(s), ERR_THRESH(s));
    save([datadir 'validationAll.mat'], 'thresh', 'valdata');    
end

if DO_MERGE_MIN
    smallsegdir = [datadir '/smallsegs/' num2str(s) '/'];
    for f = testind %numel(fn)

        disp(num2str(f))
        tmp = load([segdir strtok(fn{f}, '.') '_seg.mat']);
        bndinfo = tmp.bndinfo;            

        X = getAllFeatures(bndinfo, imdir, pbdir, gcdir);
        result = mergeStageMin(X, bndinfo, dtBnd(end), dtBnd_fast(end), 2, thresh(s));
        bndinfo = updateBoundaryInfo2(bndinfo, result);     
        
        save([smallsegdir strtok(fn{f}, '.') '_seg.mat'], 'bndinfo');    
        
        im = im2double(imread([imdir fn{f}]));
        lim = displayOcclusionResult(im, bndinfo, [], []);
        imwrite(lim, [smallsegdir 'display/' strtok(fn{f}, '.') '_res.jpg']);
    end
end
     
if DO_BP_VALIDATION       
    [thresh(s), valdata(s)] = validateMergeBp(X(ftrain), ...
        bndinfo1(ftrain), dtBnd(s), dtBnd_fast(s), dtCont(s), ERR_THRESH(s));
    save([datadir 'validationAll.mat'], 'thresh', 'valdata');    
end


if DO_MERGE_MRF    
    lastdir = [datadir '/smallsegs/' num2str(s-1) '/'];
    smallsegdir = [datadir '/smallsegs/' num2str(s) '/'];
    for f = ftrain
        
        disp(num2str(find(f==ftrain)))
        tmp = load([lastdir strtok(fn{f}, '.') '_seg.mat']);
        bndinfo = tmp.bndinfo;                   

        X = getAllFeatures(bndinfo, imdir, pbdir, gcdir);
        result = mergeStageBp2(X, bndinfo, dtBnd(end), dtBnd_fast(end), dtCont(end), thresh(s));           
        bndinfo = updateBoundaryInfo2(bndinfo, result);    
        
        save([smallsegdir strtok(fn{f}, '.') '_seg.mat'], 'bndinfo');    
        im = im2double(imread([imdir fn{f}]));
        lim = displayOcclusionResult(im, bndinfo, [], []);
        imwrite(lim, [smallsegdir 'display/' strtok(fn{f}, '.') '_res.jpg']);        
    end            
end    


if DO_MERGE_MRF_FINAL
    lastdir = [datadir '/smallsegs/' num2str(s-1) '/'];   
    smallsegdir = [datadir '/smallsegs/' num2str(s) '_2' '/'];
    for f = testind(151:250)
        disp(num2str(find(f==testind)))
        tmp = load([lastdir strtok(fn{f}, '.') '_seg.mat']);
            bndinfo = tmp.bndinfo;                   

        thresh = 0.6;       
        
        cvnum = 1;
        if any(f==testind)
            cvnum = ceil(find(f==testind)/50);
        end
            
        dorestimate = 1;
        [bndinfo, lab, plab_e, plab_g] = ...
            mergeStageFinalBpGeometry(bndinfo, dtBnd(s), dtBnd_fast(s), dtCont(s), ...
            thresh, imdir, pbdir, gcdir, gdatadir, cvnum, dorestimate);  
        lab = lab{1};               
        
        bndinfo.edges.boundaryType = lab;
        bndinfo.result.edgeProb = plab_e{1};
        bndinfo.result.geomProb = plab_g{1};
        bndinfo.result.boundaries = lab;
        bndinfo.result.thresh = thresh;
        %bndinfo2 = transferSuperpixelLabels(bndinfo2, bndinfo.wseg);
                
        save([smallsegdir strtok(fn{f}, '.') '_seg.mat'], 'bndinfo');    
        
        
        im = im2double(imread([imdir bndinfo.imname]));        
                              
        outdir = [smallsegdir 'display/'];
        imwrite(im, [outdir fn{f}]);
        printOcclusionResult(im, bndinfo, lab, [outdir strtok(fn{f}, '.') '_res.jpg'], 1);
        %imwrite(gtim, [outdir strtok(fn{f}, '.') '_res.jpg'], 'Quality', 99);
        %imwrite(label2rgb(bndinfo.wseg), [smallsegdir '/3/display/' strtok(fn{f}, '.') '_seg.jpg'], 'Quality', 99);
    end            
end


if DO_GEOMETRY_ONLY_SEG
    lastdir = './iccvGroundTruth/segs/';   
    outdir = [datadir '/geomsegs/'];
    for f = testind(1:50)
        disp(num2str(find(f==testind)))
        tmp = load([lastdir strtok(fn{f}, '.') '_seg.mat']);
        bndinfo = tmp.bndinfo;                       
                    
        bndinfo = geometry2seg(bndinfo, gcdir);                
        save([outdir strtok(fn{f}, '.') '_seg.mat'], 'bndinfo');    
                
        im = im2double(imread([imdir bndinfo.imname]));        
                              
        dispdir = [outdir 'display/'];
        printOcclusionResult(im, bndinfo, [], [dispdir strtok(fn{f}, '.') '_res.jpg'], 1);

    end            
end

if DO_NCUT_SEG
    ncutdir = [datadir '/ncutsegs100/'];
    outdir = [datadir '/ncutsegs100/'];
    for f = testind(1:50)
        disp(num2str(find(f==testind)))
     
        load([ncutdir strtok(fn{f}, '.') '_ncut.mat']);
        [tmp1, tmp2, segimage(:)] = unique(double(segimage(:)));
        
%         segimage = double(segimage);
%         segimage = segimage - min(segimage(:)) + 1;
        [edges, juncts, neighbors, wseg] = seg2fragments(segimage, [], 0);
        bndinfo = processBoundaryInfo3(wseg, edges, neighbors);
        bndinfo.imname = fn{f};        
                     
        save([outdir strtok(fn{f}, '.') '_seg.mat'], 'bndinfo');    
                
        im = im2double(imread([imdir bndinfo.imname]));        
                              
        dispdir = [outdir 'display/'];
        %imwrite(im, [outdir fn{f}]);
        printOcclusionResult(im, bndinfo, [], [dispdir strtok(fn{f}, '.') '_res.jpg'], 1);
        %imwrite(gtim, [outdir strtok(fn{f}, '.') '_res.jpg'], 'Quality', 99);
        %imwrite(label2rgb(bndinfo.wseg), [smallsegdir '/3/display/' strtok(fn{f}, '.') '_seg.jpg'], 'Quality', 99);
    end            
end

clear DO_LOAD DO_READ  DO_RESTART DO_SMALL_SEGS DO_FEATURES DO_FEATURES_GEOMETRY 
clear DO_GEOMETRY_PARAMS  DO_BOUNDARY_CLASSIFIER DO_CONTINUITY_CLASSIFIER
clear DO_VALIDATION   DO_MERGE_MIN   DO_MERGE_MRF DO_MERGE_MRF_FINAL 
clear DO_UNARY_VALIDATION DO_BP_VALIDATION DO_GEOMETRY_ONLY_SEG 

