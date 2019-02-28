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

    %% T junction feature info
    Tjinfo = {};

    for k = 1:size(jPos, 1)
        adjacentEdgeIndexs = find(ejunction == k);

        if numel(adjacentEdgeIndexs) ~= 3
            % not T-junction
            continue;
        end

        % XY = cell(3, 1)
        edgeId = zeros([3, 1], 'single');
        edgeFlip = zeros([3, 1], 'single');
        atan2d_value = zeros([3, 1], 'single');
        angles = zeros([3, 2], 'single');
        edgeConvexityFeature = zeros([3, 36], 'single');

        for k1 = 1:numel(adjacentEdgeIndexs)
            [i, j] = ind2sub(size(ejunction), adjacentEdgeIndexs(k1));
            [edgeY, edgeX] = ind2sub(imsize(1:2), double(eindices{i}));
            edgeXY = [edgeX, edgeY];

            if j == 2
                edgeXY = flip(edgeXY);
                edgeFlip(k1, 1) = 1;
            end

            edgeId(k1, 1) = i;
            atan2d_value(k1, :) = clacAtan2d_value(edgeXY);
            edgeConvexityFeature(k1, :) = getEdgeConvexityFeature(edgeXY);
        end

        if find(isnan(atan2d_value) == 1)
            continue;
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
            atan2d_value([2, 3], :) = atan2d_value([3, 2], :);
            % angles([2, 3], :) = angles([3, 2], :);
            edgeConvexityFeature([2, 3], :) = edgeConvexityFeature([3, 2], :);
        end

        angles(1, :) = atan2d_value(1) - atan2d_value([2, 3]);
        angles(2, :) = atan2d_value(2) - atan2d_value([3, 1]);
        angles(3, :) = atan2d_value(3) - atan2d_value([1, 2]);
        angles = mod(angles + 360, 360) / 360;                                                                                     
        adjacentEdgeInfo = [];
        adjacentEdgeInfo.edgeId = num2cell(edgeId);
        adjacentEdgeInfo.edgeFlip = num2cell(edgeFlip);
        % adjacentEdgeInfo.XY = num2cell(XY);
        adjacentEdgeInfo.atan2d_value = num2cell(atan2d_value);
        adjacentEdgeInfo.angles = num2cell(angles);
        adjacentEdgeInfo.edgeConvexityFeature = num2cell(edgeConvexityFeature);
        Tjinfo{end + 1, 1} = adjacentEdgeInfo;
    end

    % edge feature info
    edgeFeatures = getEdgeFeatures(bndinfo, bndinfo.pbim);

    combinedFeatures.TJInfo = Tjinfo;
    combinedFeatures.TJnum = size(Tjinfo, 1);
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
    end

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

function angle = clacVectorsAngle(v1, v2)
    %% calclate angle between two vectors
    angle = acos(dot(v1, v2) / norm(v1) * norm(v2));
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
