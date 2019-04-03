function [X, Y, gparams] = getAllFeatures(bndinfo, imdir, pbdir, gcdir, ...
        gdatadir, segmaps, cvnum, gparams)
    % [X, Y] = getAllFeatures(bndinfo, imdir, pbdir, gcdir, gdatadir, gclassifiers, segmaps)
    disp(nargin);
    for f = 1:numel(bndinfo)

        im = im2double(imread(fullfile(imdir, [bndinfo(f).imname])));

        bn = strtok(bndinfo(f).imname, '.');

        tmp = load(fullfile(pbdir, [bn '_pb.mat']));
        pball = tmp.pbim;

        tmp = load(fullfile(gcdir, [bn '.c.mat']));
        gconf = tmp.cimages{1}( :, :, 1:7);

        [X(f), Y{f}] = getFeatures(bndinfo(f), im, pball, gconf);

    end

    if nargin > 4
        fn = {bndinfo(:).imname};
        X = updateGeometricFeatures2(fn, X, segmaps, cvnum, gdatadir);

        if ~exist('gparams', 'var')
            gparams = estimateGeometryRegression(bndinfo, X, gdatadir);
        end

        %    disp('skipping geom update!')
        % XXX hardcoded for now
        %gparams = [1.36 0.10]; % stage 2
        %gparams = [1.33 0.16]; % stage 3
        for f = 1:numel(X)

            pg1 = X(f).region.pg1;
            pg1 = [pg1(:, 1) sum(pg1(:, 2:4), 2) pg1(:, 5:7)];
            pg2 = X(f).region.pg2;
            pg2 = [pg2(:, 1) sum(pg2(:, 2:4), 2) pg2(:, 5:7)];
            pg = exp(gparams(1) * log(pg1) + gparams(2) * log(pg2));
            pg = pg ./ repmat(sum(pg, 2), [1 size(pg, 2)]);

            X(f).region.geomContext = pg;

        end

    end

end
