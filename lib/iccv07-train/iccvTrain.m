% iccvTrain


%% Settings

STAGES = [3 4 5 6];

DO_READ                  = 0;
DO_RESTART               = 0;
DO_GEOMETRY              = 1;
DO_FEATURES              = 1;
DO_BOUNDARY_CLASSIFIER   = 1;
DO_VALIDATION            = 1;
DO_NEXT_STAGE            = 1;

ERR_THRESH = [0.01 ; 0.005 ; 0.005 ; 0.005 ; 0.005 ; 0.005];

ftrain = setdiff(clusterind, 230);
geomcvtrain = cv_images;
geomncv = 5;
      
global DO_DISPLAY;
DO_DISPLAY = 0;

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
    load(['./data/tmp/segmaps' num2str(s-1) '.mat'], 'segmaps'); 
    load(['./data/trainBndinfo' num2str(s-1) '.mat']);
    bndinfo1 = bndinfo2; 
    clear bndinfo2;
    
    for f = geomcvtrain
        disp(num2str(find(f==geomcvtrain)));
        tmp = load([segdir strtok(fn{f}, '.') '_seg.mat']);
        load(['./data/tmp/' strtok(fn{f}, '.'), '_lastbndinfo.mat']);
        savedata.result = savedata.result(1:s-1);
        savedata.bndinfo1 = tmp.bndinfo;
        for k = 1:s-1
            savedata.bndinfo1 = updateBoundaryInfo2(savedata.bndinfo1, savedata.result(k));
        end
        save(['./data/tmp/' strtok(fn{f}, '.'), '_lastbndinfo.mat'], 'savedata');
    end
    clear tmp;
end

for s = STAGES

    disp(['STAGE: ' num2str(s)]);
    
    if DO_GEOMETRY
        disp('Training geometry...');    
        if s==1,    segmaps = cell(size(fn));   end        
        gclassifiers(s) = iccvTrainGeometry(fn(geomcvtrain), ...
            segmaps(geomcvtrain), geomncv, gdatadir);
        save('./data/geomClassifiersAll.mat', 'gclassifiers');                  
    end
    
    if DO_FEATURES
        disp('Computing features...');    
        tmpgclass.vclassifier = gclassifiers(s).vclassifier(1);
        tmpgclass.hclassifier = gclassifiers(s).hclassifier(1);
        [X(ftrain), Y(ftrain)] = getAllFeatures(bndinfo1(ftrain), ...
            imdir, pbdir, gcdir, gdatadir, tmpgclass, segmaps(ftrain));
    end    
    
    if DO_BOUNDARY_CLASSIFIER
        cind = (1:numel(ftrain));
        disp('Training boundary...');
        dtBnd(s) = trainBoundaryClassifier4(X(ftrain(cind)), Y(ftrain(cind)), bndinfo1(ftrain(cind)));
        save('./data/bndClassifiersAll.mat', 'dtBnd');
    end

    if DO_VALIDATION
        vind = (1:numel(ftrain));
        [thresh(s), valdata(s)] = validateClassifier(X(ftrain(vind)), ...
            bndinfo1(ftrain(vind)), dtBnd(s), ERR_THRESH(s));
        save('./data/validationAll.mat', 'thresh', 'valdata');
    end        
    
    if DO_NEXT_STAGE
        disp('Merging and updating bndinfo...')
        result(ftrain,s) = mergeStage(X(ftrain), bndinfo1(ftrain), dtBnd(s), thresh(s));
        for f = ftrain                  
            [bndinfo2(f), terr(f)] = updateBoundaryInfo2(bndinfo1(f), result(f,s));      
            segmaps{f} = getSegmentationMap(bndinfo(f).nseg, result(f, 1:s));
        end         
        disp(['Mean transfer error: ' num2str(mean(terr(ftrain)))])
        save(['./data/trainBndinfo' num2str(s) '.mat'], 'bndinfo2', 'terr', 'result');
        bndinfo1 = bndinfo2;
        clear bndinfo2;
    
        % Load existing info for each geometry training image, compute
        % features, and apply next merging stage.
        for f = geomcvtrain
            disp(['next stage: ' num2str(find(f==geomcvtrain))])
            if s > 1
                load(['./data/tmp/' strtok(fn{f}, '.'), '_lastbndinfo.mat']);
            else
                tmp = load([segdir strtok(fn{f}, '.') '_seg.mat']);
                savedata.bndinfo1 = orderfields(tmp.bndinfo);
                savedata.nsp = tmp.bndinfo.nseg;
            end
            c = ceil(find(f==geomcvtrain)/(numel(geomcvtrain)/geomncv));
            tmpgclass.vclassifier = gclassifiers(s).vclassifier(c);
            tmpgclass.hclassifier = gclassifiers(s).hclassifier(c);
            tmpX = getAllFeatures(savedata.bndinfo1, ...
                imdir, pbdir, gcdir, gdatadir, tmpgclass, segmaps(f));            
            savedata.result(s) = mergeStage(tmpX, savedata.bndinfo1, dtBnd(s), thresh(s));
            savedata.bndinfo1 = updateBoundaryInfo2(savedata.bndinfo1, savedata.result(s));
            save(['./data/tmp/' strtok(fn{f}, '.'), '_lastbndinfo.mat'], 'savedata');            
            segmaps{f} = getSegmentationMap(savedata.nsp, savedata.result(1:s));
        end
        save(['./data/tmp/segmaps' num2str(s) '.mat'], 'segmaps');            
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