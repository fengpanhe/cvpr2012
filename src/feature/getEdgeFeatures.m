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
                [ch, convarea] = convhull(x, y);
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