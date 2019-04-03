%{
% Filename: trainBoundaryClassifier.m
% Project: boundary
% Created Date: Tuesday March 26th 2019
% Author: Feng Panhe
% -----
% Last Modified:
% Modified By:
% -----
% Copyright (c) 2019 Feng Panhe
%}
function trainBoundaryClassifier()
    %myFun - Description
    %
    % Syntax: trainBoundaryClassifier(input)
    %
    % Long description
    % iccvTrain

    %% Settings

    gcdLoad;

    STAGES = (1:3);

    DO_LOAD = 0;
    DO_READ = 1;
    DO_RESTART = 0;
    DO_FEATURES = 1;
    DO_FEATURES_GEOMETRY = 1;
    DO_GEOMETRY_PARAMS = 1;
    DO_BOUNDARY_CLASSIFIER = 1;
    DO_CONTINUITY_CLASSIFIER = 1;
    DO_UNARY_VALIDATION = 1;
    DO_MERGE_MIN = 1;
    DO_MERGE_MRF = 0;
    DO_MERGE_MRF_FINAL = 0;
    DO_BP_VALIDATION = 0;

    ERR_THRESH = [0.01; 0.005; 0.005; 0.005; 0.005; 0.005];

    ftrain = setdiff(clusterind, 230);
    testind = cv_images(1:50);
    geomcvtrain = cv_images;
    geomncv = 5;
    NSMALLSEGS = 300;

    thresh = [0.105 0.25];

    % data_dir = fullfile('result', 'ClassifierData');

    if ~exist(data_dir, 'dir')
        mkdir(data_dir);
    end

    global DO_DISPLAY;
    DO_DISPLAY = 0;

    %% Training

    if DO_LOAD
        load([data_dir '/bndClassifiersAll3.mat']);
        load([data_dir '/contClassifiersAll3.mat']);
    end

    if DO_READ
        disp(['Reading from ' gt_dir]);

        for f = ftrain
            tmp = load(fullfile(gt_dir, [strtok(fn{f}, '.') '_gt.mat']));
            bndinfo_list(f) = tmp.bndinfo;
        end

        clear tmp;
    end

    if DO_RESTART
        s = STAGES(1);
        disp(['Restarting from after stage ' num2str(s - 1)]);
        %load(['./data/tmp/segmaps' num2str(s-1) '.mat'], 'segmaps');
        load([data_dir 'trainBndinfo' num2str(s - 1) '.mat']);
        bndinfo1 = bndinfo2;
        clear bndinfo2;
    end

    for s = STAGES

        disp(['STAGE ' num2str(s)]);

        if s > 0

            if DO_FEATURES && s < 3
                disp('Computing features...');
                [X(ftrain), Y(ftrain)] = getAllFeatures(bndinfo_list(ftrain), images_dir, pb_dir, gc_dir);
            end

            geomparams = [];

            if DO_GEOMETRY_PARAMS && s == 3
                disp('Estimating parameters for geometry regression');
                geomparams = estimateGeometryRegression(bndinfo_list(ftrain), X(ftrain), gdatadir);
                disp(num2str(geomparams));
            end

            if DO_FEATURES_GEOMETRY && s == 3
                disp('Computing features with geometry...');
                segmaps = bndinfo2segmaps(bndinfo_list(ftrain), segdir);
                [X(ftrain), Y(ftrain)] = getAllFeatures(bndinfo_list(ftrain), images_dir, pb_dir, ...
                    gc_dir, gdatadir, segmaps, 1, geomparams);
            end

        end

        if DO_BOUNDARY_CLASSIFIER
            disp('Training boundary...');
            [dtBnd{s}, xtrain] = trainBoundaryClassifier4(X(ftrain), Y(ftrain), bndinfo_list(ftrain));
            dtBnd_fast{s} = trainFastBoundaryClassifier(X(ftrain), Y(ftrain), bndinfo_list(ftrain), 1);
            save(fullfile(data_dir, ['bndClassifiersAll' num2str(s) '.mat']), 'dtBnd', 'dtBnd_fast', 'geomparams');
        end

        if DO_CONTINUITY_CLASSIFIER

            disp('Training continuity classifier...');

            for f = ftrain
                [pB{f}, trainx{f}] = useBoundaryClassifier(bndinfo_list(f), X(f), dtBnd{s});
                pB{f} = [pB{f}(:, [1 2]); pB{f}(:, [1 3])];
            end

            dtCont{s} = trainBoundaryContinuityClassifier(bndinfo_list(ftrain), X(ftrain), Y(ftrain), pB(ftrain), trainx(ftrain));
            save(fullfile(data_dir, ['contClassifiersAll' num2str(s) '.mat']), 'dtCont');
        end

        if DO_UNARY_VALIDATION && s < 3
            [thresh(s), valdata(s)] = validateMergeMinSafe(X(ftrain), ...
                bndinfo_list(ftrain), dtBnd{s}, dtBnd_fast{s}, ERR_THRESH(s));
            save(fullfile(data_dir, 'validationAll.mat'), 'thresh', 'valdata');
        end

        if s == 3
            thresh(s) = 0.5;
        end

        if DO_MERGE_MIN

            if s == 1
                lastdir = segdir;
            else
                lastdir = fullfile(data_dir, 'smallsegs', num2str(s - 1));
            end

            smallsegdir = fullfile(data_dir, 'smallsegs', num2str(s));

            if ~exist(lastdir, 'dir')
                mkdir(lastdir);
            end

            if ~exist(smallsegdir, 'dir')
                mkdir(smallsegdir);
            end

            imind = [ftrain testind];
            if s == 3, imind = testind; end

            for f = imind%numel(fn)

                disp(num2str(f))
                tmp = load(fullfile(lastdir, [strtok(fn{f}, '.') '_seg.mat']));
                bndinfo = tmp.bndinfo;

                if s < 3
                    X(f) = getAllFeatures(bndinfo, images_dir, pb_dir, gc_dir);
                    result = mergeStageMin(X(f), bndinfo, dtBnd{end}, dtBnd_fast{end}, 2, thresh(s));
                    bndinfo = updateBoundaryInfo2(bndinfo, result);
                else % final stage

                    while 1
                        initnsp = bndinfo.nseg;
                        segmaps = bndinfo2segmaps(bndinfo, segdir);
                        X(f) = getAllFeatures(bndinfo(ftrain), images_dir, pb_dir, ...
                            gc_dir, gdatadir, segmaps, 1, geomparams);
                        result = mergeStageMin(X(f), bndinfo, dtBnd{end}, dtBnd_fast{end}, 2, thresh(s));
                        bndinfo = updateBoundaryInfo2(bndinfo, result);
                        finalnsp = numel(result.regions);

                        if finalnsp == initnsp || finalnsp <= 2
                            break;
                        end

                    end

                end

                save(fullfile(smallsegdir, [strtok(fn{f}, '.') '_seg.mat']), 'bndinfo');

                % im = im2double(imread(fullfile(images_dir, fn{f})));
                % lim = displayOcclusionResult(im, bndinfo, [], [], 1);
                % imwrite(lim, fullfile(smallsegdir, 'display', [strtok(fn{f}, '.') '_res.jpg']));
            end

        end

        if DO_BP_VALIDATION
            [thresh(s), valdata(s)] = validateMergeBp(X(ftrain), ...
                bndinfo1(ftrain), dtBnd{s}, dtBnd_fast{s}, dtCont{s}, ERR_THRESH(s));
            save([data_dir 'validationAll.mat'], 'thresh', 'valdata');
        end

        if DO_MERGE_MRF
            lastdir = [data_dir '/smallsegs/' num2str(s - 1) '/'];
            smallsegdir = [data_dir '/smallsegs/' num2str(s) '/'];

            for f = ftrain

                disp(num2str(find(f == ftrain)))
                tmp = load([lastdir strtok(fn{f}, '.') '_seg.mat']);
                bndinfo = tmp.bndinfo;

                X = getAllFeatures(bndinfo, images_dir, pb_dir, gc_dir);
                result = mergeStageBp2(X, bndinfo, dtBnd{end}, dtBnd_fast{end}, dtCont{end}, thresh(s));
                bndinfo = updateBoundaryInfo2(bndinfo, result);

                save([smallsegdir strtok(fn{f}, '.') '_seg.mat'], 'bndinfo');
                im = im2double(imread([images_dir fn{f}]));
                lim = displayOcclusionResult(im, bndinfo, [], []);
                imwrite(lim, [smallsegdir 'display/' strtok(fn{f}, '.') '_res.jpg']);
            end

        end

        if DO_MERGE_MRF_FINAL
            lastdir = [data_dir '/smallsegs/' num2str(s - 1) '/'];
            smallsegdir = [data_dir '/smallsegs/' num2str(s) '_2'' / '];

            for f = testind(151:250)
                disp(num2str(find(f == testind)))
                tmp = load([lastdir strtok(fn{f}, '.') '_seg.mat']);
                bndinfo = tmp.bndinfo;

                thresh = 0.6;

                cvnum = 1;

                if any(f == testind)
                    cvnum = ceil(find(f == testind) / 50);
                end

                dorestimate = 1;
                [bndinfo, lab, plab_e, plab_g] = ...
                    mergeStageFinalBpGeometry(bndinfo, dtBnd{s}, dtBnd_fast{s}, dtCont{s}, ...
                    thresh, images_dir, pb_dir, gc_dir, gdatadir, cvnum, dorestimate);
                lab = lab{1};

                bndinfo.edges.boundaryType = lab;
                bndinfo.result.edgeProb = plab_e{1};
                bndinfo.result.geomProb = plab_g{1};
                bndinfo.result.boundaries = lab;
                bndinfo.result.thresh = thresh;
                %bndinfo2 = transferSuperpixelLabels(bndinfo2, bndinfo.wseg);

                save([smallsegdir strtok(fn{f}, '.') '_seg.mat'], 'bndinfo');

                im = im2double(imread([images_dir bndinfo.imname]));

                outdir = [smallsegdir 'display/'];
                imwrite(im, [outdir fn{f}]);
                printOcclusionResult(im, bndinfo, lab, [outdir strtok(fn{f}, '.') '_res.jpg'], 1);
                %imwrite(gtim, [outdir strtok(fn{f}, '.') '_res.jpg'], 'Quality', 99);
                %imwrite(label2rgb(bndinfo.wseg), [smallsegdir '/3/display/' strtok(fn{f}, '.') '_seg.jpg'], 'Quality', 99);
            end

        end

    end

    clear DO_LOAD DO_READ DO_RESTART DO_SMALL_SEGS DO_FEATURES DO_FEATURES_GEOMETRY
    clear DO_GEOMETRY_PARAMS DO_BOUNDARY_CLASSIFIER DO_CONTINUITY_CLASSIFIER
    clear DO_VALIDATION DO_MERGE_MIN DO_MERGE_MRF DO_MERGE_MRF_FINAL
    clear DO_UNARY_VALIDATION DO_BP_VALIDATION DO_GEOMETRY_ONLY_SEG

end
