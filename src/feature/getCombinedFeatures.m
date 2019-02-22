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

            dv = clacDirectionVector(edgeXY(1:2,:));
            ecf = getEdgeConvexityFeature(edgeXY);

            edgeId(k1, 1) = i;
            % XY(k1, 1)  = edgeXY;
            directionVector(k1, :) = dv;
            edgeConvexityFeature(k1, :) = ecf;
        end

        angle1 = clacVectorsAngle(directionVector(1, :), directionVector(2, :));
        angle2 = clacVectorsAngle(directionVector(1, :), directionVector(3, :));
        angle3 = clacVectorsAngle(directionVector(2, :), directionVector(3, :));

        angles(1, :) = [angle1, angle2];
        angles(2, :) = [angle1, angle3];
        angles(3, :) = [angle2, angle3];

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

        if clockwiseAngle1_3 < clockwiseAngle1_2
            edgeId([2, 3], :) = edgeId([3, 2], :);
            edgeFlip([2, 3], :) = edgeFlip([3, 2], :);
            directionVector([2, 3], :) = directionVector([3, 2], :);
            angles([2, 3], :) = angles([3, 2], :);
            edgeConvexityFeature([2, 3], :) = edgeConvexityFeature([3, 2], :);
        end

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
    X = pos(:, 1);
    Y = pos(:, 2);
    p = polyfit(X, Y, 1);
    vector = [1, p(1)];

    if mean(X) < X(1)
        vector = -vector;
    end

    directionVector = vector / norm(vector);
end

function angle = clacVectorsAngle(v1, v2)
    %% calclate angle between two vectors
    angle = acos(dot(v1, v2) / norm(v1) * norm(v2));
end

function ECFeature = getEdgeConvexityFeature(edgeXY)
    %% get the convexity feature of edge

    lb = edgeXY(end, :) - edgeXY(1, :);
    lb_atan2d = atan2d(lb(2), lb(1));

    thetas = zeros(1, size(edgeXY, 1));

    for i = 2:size(edgeXY, 1)
        li = edgeXY(i, :) - edgeXY(1, :);
        li_atan2d = atan2d(li(2), li(1));
        theta_i2b = li_atan2d - lb_atan2d;

        if theta_i2b > 180
            theta_i2b = theta_i2b - 360;
        end

        if theta_i2b <- 180
            theta_i2b = theta_i2b + 360;
        end

        thetas(i) = theta_i2b;
    end

    edges = [-180, -170, -160, -150, -140, -130, -120, -110, -100, -90, -80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180];

    ECFeature = histcounts(thetas, edges);
    ECFeature = ECFeature / sum(ECFeature);
end
