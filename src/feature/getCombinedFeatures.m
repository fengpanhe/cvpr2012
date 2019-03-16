%{
% Filename: getCombinedFeatures.m
% Project: feature
% Created Date: Tuesday January 15th 2019
% Author: Feng Panhe
% -----
% Last Modified:
% Modified By:
% -----
% Copyright (c) 2019 Feng Panhe
%}

function combinedFeatures = getCombinedFeatures(bndinfo, im)
    %%
    imsize = bndinfo.imsize;
    eindices = bndinfo.edges.indices;
    ejunction = bndinfo.edges.junctions;
    jPos = bndinfo.junctions.position;
    edges_spLR = bndinfo.edges.spLR;

    edges_XY = cell(numel(eindices), 1);

    for i = 1:numel(eindices)
        [edge_Y, edge_X] = ind2sub(imsize(1:2), double(eindices{i}));
        edges_XY{i} = [edge_X, edge_Y];
    end

    %% T junction feature info

    j_edge_nums = accumarray([ejunction(:, 1); ejunction(:, 2)], 1);
    tjs_id = find(j_edge_nums == 3);
    tj_num = numel(tjs_id);
    % tj_edge_num = 3 * tj_num;
    % edges_id = zeros([tj_edge_num, 1], 'single');
    % edges_filp = zeros([tj_edge_num, 1], 'single');
    % edges_segid = zeros([tj_edge_num, 1], 'single');
    % edges_atan2d_value = zeros([tj_edge_num, 1], 'single');
    % edges_angle = zeros([tj_edge_num, 2], 'single');
    % edges_convexity = zeros([tj_edge_num, 36], 'single');
    tjs_info = cell(tj_num, 1);

    for k = 1:tj_num
        [ej_rows, ej_cols] = find(ejunction == tjs_id(k));
        tjs_info{k} = getTjinfo(ej_rows, ej_cols, edges_XY(ej_rows), edges_spLR(ej_rows, :));
    end

    tjs_cell_len = cellfun(@length, tjs_info);
    tjs_info(~logical(tjs_cell_len)) = [];
    % edge feature info
    edgeFeatures = getEdgeFeatures(bndinfo, bndinfo.pbim);

    combinedFeatures.TJInfo = tjs_info;
    combinedFeatures.TJnum = size(tjs_info, 1);
    combinedFeatures.edgeInfo = edgeFeatures;
end

function atan2d_value = clacAtan2d_value(pos)
    %% Calculate the direction vector of a set of points
    start_i = 1;

    for i = 2:size(pos, 1)

        if pos(i, :) == pos(start_i, :)
            start_i = i;
        end

    end

    if start_i == size(pos, 1)
        atan2d_value = NaN;
    else

        pos = pos(start_i:end, :);
        thetas = atan2d(-(pos(:, 2) - pos(1, 2)), pos(:, 1) - pos(1, 1));
        thetas = round(thetas);
        unique_thetas = unique(thetas);
        thetas_size = numel(thetas);
        max_thetas_score = 0;
        atan2d_value = 0;

        for i = 1:numel(unique_thetas)
            indexs = find(thetas == unique_thetas(i));
            thetas_score = sum(thetas_size - indexs);

            if thetas_score > max_thetas_score
                atan2d_value = unique_thetas(i);
                max_thetas_score = thetas_score;
            end

        end

    end

end

function ECFeature = getEdgeConvexityFeature(edgeXY)
    %% get the convexity feature of edge

    thetas = atan2d(-(edgeXY(:, 2) - edgeXY(1, 2)), edgeXY(:, 1) - edgeXY(1, 1));
    thetas = thetas - thetas(end);
    thetas = mod(thetas + 180, 360) - 180;

    edges = [-180, -170, -160, -150, -140, -130, -120, -110, -100, -90, -80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180];

    ECFeature = histcounts(thetas, edges);
    ECFeature = ECFeature / sum(ECFeature);
end

function tj_info = getTjinfo(ej_rows, ej_cols, edges_XY, edges_spLR)
    %%  getTjinfo
    edgeId = zeros([3, 1], 'single');
    edgeFlip = zeros([3, 1], 'single');
    edge_segid = zeros([3, 1], 'single');
    atan2d_value = zeros([3, 1], 'single');
    angles = zeros([3, 2], 'single');
    edgeConvexityFeature = zeros([3, 36], 'single');

    for k1 = 1:numel(ej_rows)
        edgeXY = edges_XY{k1};

        if ej_cols(k1) == 2
            edgeXY = flip(edgeXY);
            edgeFlip(k1, 1) = 1;
            edge_segid(k1, 1) = edges_spLR(k1, 2);
        else
            edge_segid(k1, 1) = edges_spLR(k1, 1);
        end

        edgeId(k1, 1) = ej_rows(k1);
        atan2d_value(k1, :) = clacAtan2d_value(edgeXY);
        edgeConvexityFeature(k1, :) = getEdgeConvexityFeature(edgeXY);
    end

    if find(isnan(atan2d_value) == 1)
        tj_info = {};
        return;
    end

    % Align three sides clockwise
    clockwiseAngle1_2 = atan2d_value(1) - atan2d_value(2);

    if clockwiseAngle1_2 < 0
        clockwiseAngle1_2 = clockwiseAngle1_2 + 360;
    end

    clockwiseAngle1_3 = atan2d_value(1) - atan2d_value(3);

    if clockwiseAngle1_3 < 0
        clockwiseAngle1_3 = clockwiseAngle1_3 + 360;
    end

    if clockwiseAngle1_3 < clockwiseAngle1_2
        edgeId([2, 3], :) = edgeId([3, 2], :);
        edgeFlip([2, 3], :) = edgeFlip([3, 2], :);
        edge_segid([2, 3], :) = edge_segid([3, 2], :);
        atan2d_value([2, 3], :) = atan2d_value([3, 2], :);
        % angles([2, 3], :) = angles([3, 2], :);
        edgeConvexityFeature([2, 3], :) = edgeConvexityFeature([3, 2], :);
    end

    angles(1, :) = atan2d_value(1) - atan2d_value([2, 3]);
    angles(2, :) = atan2d_value(2) - atan2d_value([3, 1]);
    angles(3, :) = atan2d_value(3) - atan2d_value([1, 2]);
    angles = mod(angles + 360, 360) / 360;
    tj_info.edgeId = num2cell(edgeId);
    tj_info.edgeFlip = num2cell(edgeFlip);
    tj_info.edge_segid = num2cell(edge_segid);
    tj_info.atan2d_value = num2cell(atan2d_value);
    tj_info.angles = num2cell(angles);
    tj_info.edgeConvexityFeature = num2cell(edgeConvexityFeature);
end

function edgeFeatures = getEdgeFeatures(bndinfo, pbim)
    % [X, Y] = getBoundaryFeatures(bndinfo, pbim, gconf)
    %
    % Computes features for each edgelet based on boundary and geometry
    % confidence images.  This version (3) computes geometry confidences as
    % means of superpixel values on either side of the edglet.
    %
    % Input:
    %   bndinfo - structure of superpixel boundaries
    %   in(imh, imw, 3) - color image
    %   pbim(imh, imw, norient) - probability of boundary image at each orientation
    %   gconf(imh, imw, ngeoms) - confidence in each geometric label
    %
    % Output:
    %   X(nboundaries, :) - features for the boundaries
    %   Y(nboundaries, 2) - geometry on each size or 0 for no edge

    edges = bndinfo.edges;
    nbnd = bndinfo.ne * 2; % number of directed edges
    ne = bndinfo.ne;
    [imh, imw] = size(bndinfo.wseg);

    edgeFeatures.pb = zeros(ne, 1);
    edgeFeatures.theta = zeros(ne, 1);
    edgeFeatures.thetaStart = zeros(ne * 2, 1);
    edgeFeatures.thetaEnd = zeros(ne * 2, 1);
    edgeFeatures.smoothness = zeros(ne, 1);
    edgeFeatures.length = zeros(ne, 1);
    edgeFeatures.convArea = zeros(ne, 1);
    edgeFeatures.convRatio = zeros(ne, 1);
    edgeFeatures.chains = zeros(ne, 1); % chains;
    edgeFeatures.edge2chain = zeros(ne, 1);
    edgeFeatures.chainsize = zeros(ne, 1);

    %% Edge statistics

    % get discretize orientation into 1=1|2, 2=1/2, 3=1_2, 4=2\1
    theta = bndinfo.edges.thetaUndirected;
    rels = (theta <- 3 * pi / 8) + (theta <- pi / 8) + (theta < pi / 8) + (theta < 3 * pi / 8);
    rels = mod(rels, 4) + 1;

    edgeFeatures.theta = bndinfo.edges.thetaDirected; % directed edge angle

    % pbim(imh, imw, [ -, \, |, / ]) (direction of edge)
    pbmap = [3 4 1 2];

    %edgeFeatures.pbOrient = zeros(ne,4);
    % compute features
    for k = 1:ne

        eind = edges.indices{k};

        edgeFeatures.length(k) = numel(eind); % edge length

        pbsubi = pbmap(rels(k));
        ind = double(eind + (pbsubi - 1) * imh * imw);
        edgeFeatures.pb(k) = sum(pbim(ind)) / numel(ind); % mean pb

        % short-range angles
        y = mod(ind - 1, imh) + 1;
        x = floor((ind - 1) / imh) + 1;
        ni = numel(ind);
        de = 10; % length of edge used to measure angle
        x1 = x([1 min(de, ni)]);
        y1 = y([1 min(de, ni)]);
        x2 = x([max(ni - de + 1, 1) ni]);
        y2 = y([max(ni - de + 1, 1) ni]);
        edgeFeatures.thetaStart(k) = atan2(-(y1(2) - y1(1)), x1(2) - x1(1));
        edgeFeatures.thetaEnd(k) = atan2(-(y2(2) - y2(1)), x2(2) - x2(1));

        edgeFeatures.smoothness(k) = (numel(ind) - 1) / (abs(x(end) - x(1)) + abs(y(end) - y(1)) + eps);

        convarea = 0;
        %segcount = [0.5 0.5];
        if 0 && numel(x) > 2

            try
                [~, convarea] = convhull(x, y);
                %             mask = poly2mask(x(ch), y(ch), imh, imw);
                %             segnums = bndinfo.wseg(mask);
                %             segcount = [sum(segnums(:)==bndinfo.edges.spLR(k, 1)) ...
                    %                         sum(segnums(:)==bndinfo.edges.spLR(k, 2))];
            catch
            end

        end

        edgeFeatures.convArea(k) = convarea / imh / imw;
        %edgeFeatures.convRatio(k) = (segcount(1)+eps) / (sum(segcount)+2*eps);
        %edgeFeatures.pbOrient(k, pbsubi) = edgeFeatures.pb(k);

    end

    edgeFeatures.thetaStart(ne + 1:end) = edgeFeatures.thetaEnd(1:ne) + pi;
    edgeFeatures.thetaEnd(ne + 1:end) = edgeFeatures.thetaStart(1:ne) + pi;

    thetaEnd = mod(edgeFeatures.thetaEnd * 180 / pi, 360);
    thetaStart = mod(edgeFeatures.thetaStart * 180 / pi, 360);
    % chain together edgelets
    [chains, e2chain, chainsize] = chainEdgelets([edgeFeatures.pb; edgeFeatures.pb], ...
        edges.adjacency, thetaStart, thetaEnd, 0.02, 45);
    edgeFeatures.chains = chains;
    edgeFeatures.edge2chain = single(e2chain);
    edgeFeatures.chainsize = single(chainsize);
end
