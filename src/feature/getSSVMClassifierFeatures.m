function itemInfos = getSSVMClassifierFeatures(bndinfo, combinedFeatures, type)
    % X is the raw data
    % ind is the set of indices for which the edgeFeatures should be computed
    %
    % infos:
    %   +1: Lable
    %   +1: Edge id
    %   +1: T junction id (group id)
    %
    %   +2: Angles
    %   +36: Convexity feature
    %
    %   +1: Pb
    %   +1: Length / Perimeter
    %   +1: Smoothness
    %   +1: Angle
    %   +2: Continuity
    %   +2: Convexity (area and ratio) - not used
    %   +1: Chain length
    %   +1: seg type

    boundarylabs = [];

    if type == 'train'
        boundarylabs = bndinfo.edges.boundaryType;
        boundarylabs = (boundarylabs(1:end / 2) > 0) + 2 * (boundarylabs(end / 2 + 1:end) > 0);
    end

    ind = (1:bndinfo.ne);
    edgeFeatures = getEdgeFeatures(bndinfo, combinedFeatures.edgeInfo, ind);

    TJInfos = combinedFeatures.TJInfo;
    TJnum = numel(TJInfos);

    features = zeros([TJnum * 3, 49], 'single');
    lables = zeros([TJnum * 3, 1], 'single');

    for k = 1:TJnum
        col = 1;
        row = k * 3 - 2:k * 3;

        tjinfo = TJInfos{k};

        % Edge id
        features(row, col) = cell2mat(tjinfo.edgeId);
        col = col + 1;

        % T junction id, 1-dimension
        features(row, col) = k;
        col = col + 1;

        % Angle, 2-dimension
        features(row, col:col + 1) = cell2mat(tjinfo.angles);
        col = col + 2;

        % convexity feature, 36-dimension
        features(row, col:col + 35) = cell2mat(tjinfo.edgeConvexityFeature);
        col = col + 36;

        %  Hoiem et al. [7] proposes local features fd(e), 9-dimension
        eIds = cell2mat(tjinfo.edgeId);
        features(row, col:col + 8) = edgeFeatures(eIds, 1:9);
        col = col + 9;

        if isfield(tjinfo, 'seg_type')
            features(row, col:col) = cell2mat(tjinfo.seg_type);
            col = col + 1;
        end

        % disp(col);
        if size(boundarylabs) ~= 0
            elabs = boundarylabs(eIds);
            eFlip = cell2mat(tjinfo.edgeFlip);

            for e_k = 1:numel(elabs)

                if eFlip(e_k) == 1
                    elabs(e_k) = (elabs(e_k) == 1) * 2 + (elabs(e_k) == 2);
                end

            end

            evalues = zeros([3, 1], 'single');
            lab_tmp = [elabs(3), elabs(1)];
            evalues(1) = all(lab_tmp == [2, 1]) * 3 + all(lab_tmp == [1, 1]) * 2 + all(lab_tmp == [2, 2]) * 2 + all(lab_tmp == [1, 2]) * 1;
            lab_tmp = [elabs(1), elabs(2)];
            evalues(2) = all(lab_tmp == [2, 1]) * 3 + all(lab_tmp == [1, 1]) * 2 + all(lab_tmp == [2, 2]) * 2 + all(lab_tmp == [1, 2]) * 1;
            lab_tmp = [elabs(2), elabs(3)];
            evalues(3) = all(lab_tmp == [2, 1]) * 3 + all(lab_tmp == [1, 1]) * 2 + all(lab_tmp == [2, 2]) * 2 + all(lab_tmp == [1, 2]) * 1;

            lables(row, 1) = evalues;
        end

    end

    itemInfos = [lables, features];
end

function edgeFeatures = getEdgeFeatures(bndinfo, edgeInfo, ind)
    %% Edge Features
    % edgeFeatures:
    %   Edge edgeFeatures (1-6)
    %        1:  Pb
    %        2:  Length / Perimeter
    %        3:  Smoothness
    %        4:  Angle
    %      5-6:  Continuity
    %      7-8:  Convexity (area and ratio) - not used
    %        9:  Chain length

    ndata = numel(ind);

    edgeFeatures = zeros([ndata 9], 'single');

    if isempty(ind)
        return;
    end

    spLR = bndinfo.edges.spLR;
    s1 = spLR(ind, 1);
    s2 = spLR(ind, 2);

    f = 0;

    edgeFeatures(:, f + 1) = edgeInfo.pb(ind);

    perim = zeros(bndinfo.nseg, 1);

    for k = 1:numel(edgeInfo.length)
        perim(spLR(k, 1)) = perim(spLR(k, 1)) + edgeInfo.length(k);
        perim(spLR(k, 2)) = perim(spLR(k, 2)) + edgeInfo.length(k);
    end

    minperim = min([perim(s1) perim(s2)], [], 2);
    edgeFeatures(:, f + 2) = edgeInfo.length(ind) ./ minperim; % edge length / perim of smaller region

    edgeFeatures(:, f + 3) = edgeInfo.smoothness(ind); % measure of smoothess

    theta = edgeInfo.theta;
    % discrete angle
    edgeFeatures(:, f + 4) = max(ceil((mod(theta(ind), 2 * pi) / (pi * 2) * 16 - 1E-10)), 1);

    % relative angle (continuity)
    theta1 = mod(edgeInfo.thetaStart * 180 / pi, 360);
    theta2 = mod(edgeInfo.thetaEnd * 180 / pi, 360);
    maxc = zeros(ndata, 2);
    eadj = bndinfo.edges.adjacency;
    ne = bndinfo.ne;

    for k = 1:ndata
        ki = ind(k);
        ra = abs(theta2(ki) - theta1(eadj{ki}));
        ra = ra - 180 * (ra > 180);

        if isempty(ra), maxc(k, 1) = 0;
        else maxc(k, 1) = min(ra);
        end

        ra = mod(abs(theta2(ne + ki) - theta1(eadj{ne + ki})), 180 + 1E-5);

        if isempty(ra), maxc(k, 2) = 0;
        else maxc(k, 2) = min(ra);
        end

    end

    edgeFeatures(:, f + (5:6)) = [min(maxc, [], 2) max(maxc, [], 2)];

    edgeFeatures(:, f + 8) = edgeInfo.convRatio;

    edgeFeatures(:, f + 9) = edgeInfo.edge2chain(ind);

    f = f + 9;
end
