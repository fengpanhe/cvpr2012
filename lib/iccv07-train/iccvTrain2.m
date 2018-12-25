% iccvTrain


%% Settings

STAGES = (3);

DO_READ                  = 0;
DO_RESTART               = 0;
DO_SMALL_SEGS            = 0;
DO_FEATURES              = 0;
DO_BOUNDARY_CLASSIFIER   = 0;
DO_CONTINUITY_CLASSIFIER = 1;
DO_VALIDATION            = 0;
DO_NEXT_STAGE            = 0;

ERR_THRESH = [0.01 ; 0.005 ; 0.005 ; 0.005 ; 0.005 ; 0.005];

ftrain = setdiff(clusterind, 230);
geomcvtrain = cv_images;
geomncv = 5;
NSMALLSEGS = 300;  

datadir = '/usr1/projects/dhoiem/iccv07/data3/';

global DO_DISPLAY;
DO_DISPLAY = 1;

%% Training

if DO_READ
    disp('Reading...');
    for f = ftrain
        tmp = load([gtdir strtok(fn{f}, '.') '_gt.mat']);
        bndinfo(f) = tmp.bndinfo;
    end
    bndinfo1 = bndinfo;
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

if DO_SMALL_SEGS
    for f = 190:300 %numel(fn)
        disp(num2str(f))
        if ismember(f, ftrain)
            tmp = load([gtdir strtok(fn{f}, '.') '_gt.mat']);
            bndinfo = tmp.bndinfo;
        else
            tmp = load([segdir strtok(fn{f}, '.') '_seg.mat']);
            bndinfo = tmp.bndinfo;            
        end

        X = getAllFeatures(bndinfo, imdir, pbdir, gcdir);
        result = mergeStageMin(X, bndinfo, dtBnd, NSMALLSEGS);                  
        bndinfo_sm = updateBoundaryInfo2(bndinfo, result);     
        
        save([smallsegdir strtok(fn{f}, '.') '_smallseg.mat'], 'bndinfo_sm');    
    end
end
           
    

for s = STAGES

    disp(['STAGE: ' num2str(s)]);
        
    if DO_FEATURES
        disp('Computing features...');    
        [X(ftrain), Y(ftrain)] = getAllFeatures(bndinfo1(ftrain), ...
            imdir, pbdir, gcdir);
    end    
    
    if DO_BOUNDARY_CLASSIFIER
        cind = (1:numel(ftrain));
        disp('Training boundary...');
        [dtBnd(s), xtrain] = ...
            trainBoundaryClassifier4(X(ftrain(cind)), Y(ftrain(cind)), bndinfo1(ftrain(cind)));
        save([datadir 'bndClassifiersAll.mat'], 'dtBnd');
    end

    if DO_CONTINUITY_CLASSIFIER
        disp('Training continuity...');
        for f = ftrain
            if ~exist('pB', 'var') || numel(pB)<f || isempty(pB{f})
                pB{f} = useBoundaryClassifier(bndinfo1(f), X(f), dtBnd(s));
                pB{f} = [pB{f}(:, [1 2]) ; pB{f}(:, [1 3])];
            end                
        end        
        dtCont(s) = trainBoundaryContinuityClassifier(...
            bndinfo1(ftrain), X(ftrain), Y(ftrain), pB(ftrain), xtrain);
        save([datadir 'contClassifiersAll.mat'], 'dtCont');
    end    
    
    if DO_VALIDATION
        vind = (1:numel(ftrain));
        [thresh(s), valdata(s)] = validateClassifier2(X(ftrain(vind)), ...
            bndinfo1(ftrain(vind)), dtBnd(s), dtCont(s), ERR_THRESH(s));
        save([datadir 'validationAll.mat'], 'thresh', 'valdata');
    end        
    
    if DO_NEXT_STAGE
        disp('Merging and updating bndinfo...')
        result(ftrain,s) = mergeStageBp2(X(ftrain), bndinfo1(ftrain), dtBnd(s), dtCont(s), thresh(s));
        for f = ftrain                  
            [bndinfo2(f), terr(f)] = updateBoundaryInfo2(bndinfo1(f), result(f,s));      
        end         
        disp(['Mean transfer error: ' num2str(mean(terr(ftrain)))])
        save([datadir 'trainBndinfo' num2str(s) '.mat'], 'bndinfo2', 'terr', 'result');
        bndinfo1 = bndinfo2;
        clear bndinfo2;
           
    end
   
end
    
    
    
    
%% Unused code    
% if DO_CONTINUITY_CLASSIFIER
%     disp('Training continuity...');
%     for f = ftrain
%         if ~exist('pB', 'var') || numel(pB)<f || isempty(pB{f})
%             %pB{f} = useBoundaryClassifier(X{f}(1:end/2, :), dtBnd, 'c');
%             pB{f} = useBoundaryClassifier(X{f}, dtBnd, 'c');
%             pB{f} = [pB{f}(:, [1 2]) ; pB{f}(:, [1 3])];
%         end
%         
%         % make it so that most likely contour type has confidence equal to
%         % probability of contour (important for MAP estimation)
%         % pC{f}(:, 2:5) = pC{f}(:, 2:5) ./ ...
%         %    repmat(max(pC{f}(:, 2:5), [], 2), [1 4]) .* ...
%         %    repmat(1-pC{f}(:, 1), [1 4]);                    
%     end        
%     dtCont = trainBoundaryContinuityClassifier(bndinfo(ftrain), X(ftrain), Y(ftrain), pB(ftrain));
%     save('./data/contClassifier2.mat', 'dtCont');
% end