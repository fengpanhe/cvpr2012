% iccvTestLowLevelScript

DO_READ                  = 0;
DO_FEATURES              = 0;
DO_TRAIN                 = 0;
DO_TEST                  = 1;

%iccvLoad;

testind = cv_images(1:50);
trainind = setdiff(clusterind, 230);



global DO_DISPLAY;
DO_DISPLAY = 0;

%% Read train data and compute features

for STAGE = 2:3
if STAGE==1
    indir = gtdir;  % for stage 1
else
    indir = [smallsegdir '/' num2str(STAGE-1) '/'];  % for stage 3 
end
if DO_READ
    disp(['Reading from ' indir ' ...']);
    for f = trainind 
        tmp = load([indir strtok(fn{f}, '.') '_gt.mat']);
        f2 = find(trainind==f);
        bndinfo(f2) = tmp.bndinfo;
    end
end

if DO_FEATURES && STAGE==1
    disp('Computing features...')
    [X, Y] = getAllFeatures(bndinfo, imdir, pbdir, gcdir);
end


if DO_FEATURES && STAGE>1
    disp('Computing features with geometry...');    
    segmaps = bndinfo2segmaps(bndinfo, segdir);
    [X, Y, gparams] = getAllFeatures(bndinfo, imdir, pbdir, ...
        gcdir, gdatadir, segmaps, 1);
end


%% Train each classifier
% tx:  (from getBoundaryClassifierFeatures5 -- now outdated)
%   Edge features (1-6)
%        1:  Pb
%        2:  Length / Perimeter
%        3:  Smoothness
%        4:  Angle
%      5-6:  Continuity
%      7-8:  Convexity (area and ratio)
%        9:  Chain length
%   Region features (7-17)+3
%    10-11:  Area
%       12:  Color Mean Difference
%       13:  Color Entropy Difference
%       14:  Gradient Entropy Difference
%    15-16:  Position (x,y)
%    17-18:  Extent Overlap (x, y)
%   Geometry features (16-39)+3
%    19-28:  Geometric Context Mean
%    29-33:  Geometric Context Difference
%    34   :  Geometric Context Sum Abs Difference
%    35-36:  Geometric Context Most Likely Label (G V or S)
%    37-40:  Depth under- and over-estimates for each side
%    31-43:  Depth, min1-min2, max1-max2, min(max12) - max(min12)
%    44-47:  Depthcol, each sp, diff, abs diff

if DO_TRAIN

    disp('Training 2d + gc classifier...');
    featureset{3} = (1:45);
    dtBnd_gc = trainBoundaryClassifier4(X, Y, bndinfo, featureset{3});    
    
    if 0
    disp('Training pb-based classifier...');
    featureset{1} = 1;
    %dtBnd_pb = trainBoundaryClassifier4(X, Y, bndinfo, featureset{1});

    disp('Training image-based classifier...');
    featureset{2} = (1:28);
    dtBnd_image = trainBoundaryClassifier4(X, Y, bndinfo, featureset{2});            
    
    disp('Training with all features....');
    featureset{4} = (1:60);
    dtBnd_all = trainBoundaryClassifier4(X, Y, bndinfo, featureset{4});
    end
    
    if 0
    save(['./results_new/classifierComparisonTrain_stage' num2str(STAGE) '.mat'], ...
        'dtBnd_image', 'dtBnd_gc', 'dtBnd_all', 'featureset');    
    else
        save(['./results_new/classifierComparisonTrain_stage' num2str(STAGE) '.mat'], ...
        'dtBnd_gc', 'featureset');           
    end
    clear X Y bndinfo
end


%% Test each classifier

if DO_TEST        
    
    disp('Getting test likelihoods ...')
    for f = testind         
        
        tmp = load([indir strtok(fn{f}, '.') '_gt.mat']);
        bndinfo = tmp.bndinfo;
        f2 = find(f==testind);
        disp(num2str(f2))
        
        if 1 || STAGE==1 % no geometry update
            [X, Y] = getAllFeatures(bndinfo, imdir, pbdir, gcdir);
        else % update geometry
            segmaps = bndinfo2segmaps(bndinfo, segdir);
            [X, Y] = getAllFeatures(bndinfo, imdir, pbdir, gcdir, gdatadir, segmaps, 1, gparams);
        end
        Y = Y{1};        
        w = X.edge.length; w = w / sum(w);
        
        pB_nrt = useBoundaryClassifier(bndinfo, X, ...
            dtBnd_all, featureset{4});
        tmp = evaluateLowLevelPerformance(Y, pB_nrt, w);
        tmp.imname = fn{f};
        tmp.imnum = f;        
        tmp.featureset = featureset{4};
        tmp.pB = pB_nrt;
        tmp.labels = Y;
        result_nrt(f2) = tmp;
        figure(1), hold on, plot(result_nrt(f2).pr.r, result_nrt(f2).pr.p, 'c');                        
        
        if 0
        pB_pb = [1-X.edge.pb X.edge.pb/2 X.edge.pb/2];
        tmp = evaluateLowLevelPerformance(Y, pB_pb, w);
        tmp.imname = fn{f};
        tmp.imnum = f;
        tmp.featureset = featureset{1};        
        tmp.pB = pB_pb;
        tmp.labels = Y;        
        tmp.w = w;
        result_pb(f2) = tmp;
        figure(1), hold off, plot(result_pb(f2).pr.r, result_pb(f2).pr.p, 'r');

        pB_image = useBoundaryClassifier(bndinfo, X, ...
            dtBnd_image, featureset{2});
        tmp = evaluateLowLevelPerformance(Y, pB_image, w);
        tmp.imname = fn{f};
        tmp.imnum = f;        
        tmp.featureset = featureset{2};
        tmp.pB = pB_image;
        tmp.labels = Y;
        result_im(f2) = tmp;
        figure(1), hold on, plot(result_im(f2).pr.r, result_im(f2).pr.p, 'g');        
        

        
        pB_all = useBoundaryClassifier(bndinfo, X, ...
            dtBnd_all, featureset{4});
        tmp = evaluateLowLevelPerformance(Y, pB_all, w); %#ok<AGROW>
        tmp.imname = fn{f};
        tmp.imnum = f;        
        tmp.featureset = featureset{4};
        tmp.pB = pB_all;
        tmp.labels = Y;        
        result_all(f2) = tmp;
        figure(1), hold on, plot(result_all(f2).pr.r, result_all(f2).pr.p, 'b');                
        drawnow;
        
        disp(['Figure/Ground accuracy: ' num2str([result_pb(f2).fgaccuracy ...
            result_im(f2).fgaccuracy result_gc(f2).fgaccuracy  result_all(f2).fgaccuracy])]);
        end
    end
    
    [pr_nrt, roc_nrt, fgacc_nrt] = summarizeLowLevelPerformance(result_nrt);
    
    if 0
    [pr_pb, roc_pb, fgacc_pb] = summarizeLowLevelPerformance(result_pb);
    [pr_im, roc_im, fgacc_im] = summarizeLowLevelPerformance(result_im);    
    [pr_all2, roc_all2, fgacc_all2] = summarizeLowLevelPerformance(result_all2);
    end
    
    disp(['Figure/Ground accuracy: ' num2str(mean([result_nrt.fgaccuracy]))]);
    for k = 1:50, ap(k) = averagePrecision(result_nrt(k).pr, (0:0.01:1)); end; 
    ap2 = averagePrecision(pr_nrt, (0:0.01:1));
    disp(['Average precision: ' num2str([mean(ap) ap2])]);
    
    save(['./results_new/classifierComparisonResult_stage' num2str(STAGE) '_nrt.mat'], ...
        'dtBnd_all', 'featureset', 'pr_nrt', 'roc_nrt', 'fgacc_nrt');      
    
    if 0 
    save(['./results_new/classifierComparisonResult_stage' num2str(STAGE) '.mat'], ...
        'dtBnd_image', 'dtBnd_gc', 'dtBnd_all', 'featureset', 'result_pb', 'result_im', ...
        'result_gc', 'result_all', 'pr_pb', 'roc_pb', 'pr_im', 'roc_im', 'pr_all', 'roc_all', ...
        'fgacc_pb', 'fgacc_im', 'fgacc_all', 'pr_gc', 'roc_gc', 'fgacc_gc');        
    end
end
   
end

%% Clean-up
clear DO_READ DO_FEATURES DO_TRAIN DO_TEST