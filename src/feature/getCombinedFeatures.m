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
        directionVector = zeros([3, 2], 'single');
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

            dv = clacDirectionVector(edgeXY);
            ecf = getEdgeConvexityFeature(edgeXY);

            edgeId(k1, 1) = i;
            % XY(k1, 1)  = edgeXY;
            directionVector(k1, :) = dv;
            edgeConvexityFeature(k1, :) = ecf;
        end

        if find(isnan(directionVector) == 1)
            continue;
        end

        % Align three sides clockwise
        atan2ds = atan2d(directionVector(:, 2), directionVector(:, 1));
        clockwiseAngle1_2 = atan2ds(1) - atan2ds(2);

        if clockwiseAngle1_2 < 0
            clockwiseAngle1_2 = clockwiseAngle1_2 + 360;
        end

        clockwiseAngle1_3 = atan2ds(1) - atan2ds(3);

        if clockwiseAngle1_3 < 0
            clockwiseAngle1_3 = clockwiseAngle1_3 + 360;
        end

        if clockwiseAngle1_3 > clockwiseAngle1_2
            edgeId([2, 3], :) = edgeId([3, 2], :);
            edgeFlip([2, 3], :) = edgeFlip([3, 2], :);
            directionVector([2, 3], :) = directionVector([3, 2], :);
            % angles([2, 3], :) = angles([3, 2], :);
            edgeConvexityFeature([2, 3], :) = edgeConvexityFeature([3, 2], :);
        end

        angle1 = clacVectorsAngle(directionVector(1, :), directionVector(2, :));
        angle2 = clacVectorsAngle(directionVector(1, :), directionVector(3, :));
        angle3 = clacVectorsAngle(directionVector(2, :), directionVector(3, :));

        angles(1, :) = [angle1, angle2];
        angles(2, :) = [angle3, angle1];
        angles(3, :) = [angle2, angle3];

        adjacentEdgeInfo = [];
        adjacentEdgeInfo.edgeId = num2cell(edgeId);
        adjacentEdgeInfo.edgeFlip = num2cell(edgeFlip);
        % adjacentEdgeInfo.XY = num2cell(XY);
        adjacentEdgeInfo.directionVector = num2cell(directionVector);
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

function directionVector = clacDirectionVector(pos)
    %% Calculate the direction vector of a set of points
    pos_x = pos(:, 1);
    pos_y = pos(:, 2);
    diff_x = pos_x - pos_x(1);
    diff_y = pos_y - pos_y(1);
    diff_x_no_zero_indexs = find(diff_x ~= 0);
    diff_y_no_zero_indexs = find(diff_y ~= 0);

    if ~isempty(diff_x_no_zero_indexs)

        index_tmp = diff_x_no_zero_indexs(1);
        k = diff_y(index_tmp) / diff_x(index_tmp);
        vector = [1, k];

        if pos_x(index_tmp) < pos_x(1)
            vector = -vector;
        end

    elseif ~isempty(diff_y_no_zero_indexs)

        index_tmp = diff_y_no_zero_indexs(1);
        vector = [0, 1];

        if diff_y(index_tmp) < 0
            vector = -vector;
        end

    else
        vector = [NaN, NaN];
    end

    directionVector = vector / norm(vector);
end

function angle = clacVectorsAngle(v1, v2)
    %% calclate angle between two vectors
    angle = acos(dot(v1, v2) / norm(v1) * norm(v2));
end

function ECFeature = getEdgeConvexityFeature(edgeXY)
    %% get the convexity feature of edge

    thetas = atan2d(-(edgeXY(:,2) - edgeXY(1,2)), edgeXY(:,1) - edgeXY(1,1));
    thetas = thetas - thetas(1);
    thetas = mod(thetas + 180, 360) - 180;

    edges = [-180, -170, -160, -150, -140, -130, -120, -110, -100, -90, -80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180];

    ECFeature = histcounts(thetas, edges);
    ECFeature = ECFeature / sum(ECFeature);
end
