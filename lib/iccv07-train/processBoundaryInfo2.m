function bndinfo = processBoundaryInfo2(bndinfo)
%
% bndinfo = processBoundaryInfo2(bndinfo)
%
% Adds several fields to the bndinfo structure based on bndinfo.wseg
% 
% Input:
%   bndinfo.wseg(imh,imw) = oversegmentation with edges set to 0; uint16
%
% Output:
%   bndinfo.
%     wseg2(imh+1,imw+1) = wseg after padding with edges and skel; uint16
%     nseg  = number of segments
%     ne = number of edglets
%     nj = number of junctions
%     edges.
%       indices{ne} = cell array of pixel indices for each edglet
%       adjacency{2*ne} = directed adjacency for each edge
%                         1..ne s.t. j1-->j2;   ne+1..2*ne s.t. j2-->j1
%       junctions(ne, 2) = junction indices for each edgelet ([j1 j2])                  
%       spLR(ne, 2) = superpixel indices to [left right] of edglet
%       thetaDirected(ne, 1) = orientation from j1 --> j2 (-pi to pi)
%       thetaUndireted(ne, 1) = undirected orientation (-pi/2 to pi/2)
%     junctions.
%       indices{nj} = cell array of pixel indices for each junction
%       adjacency{nj}(nadj 2) = [edgelet nextJunction] for each adjacency
%       position(nj, 2) = mean [x y] position for each junction
%
% Notes:
%   1) Edglet indices 1..ne traverse from j1 --> j2.  Edglet indices 
%   ne+1 --> 2*ne traverse from j2 --> j1.  Some statistics for latter
%   (e.g. thetaDirected) are not stored, as they can be easily calculated
%   from forward direction.  
%   2) All indices and positions are wrt wseg2.
%


%% Get edgelets and junctions of edgelets

% get junctions (cells of pixel indices), 
% junction image, edge image, expanded wseg
[junctions, jim, eim, wseg] = getJunctions(bndinfo.wseg);

% get edgelets (cells of edge indices)
[edges, eim] = edgelinkRevised2(eim, jim, junctions);
%Arguments:  eim        - Binary edge image, it is assumed that edges
%                          have been thinned.
%             jim        - Junction image (non-zero where juncts exist)   
%             minlength  - Optional minimum edge length of interest, defaults
%                          to 1 if omitted or specified as [].
%             location   - Optional complex valued image holding subpixel
%                          locations of edge points. For any pixel the
%                          real part holds the subpixel row coordinate of
%                          that edge point and the imaginary part holds
%                          the column coordinate.  See NONMAXSUP.  If
%                          this argument is supplied the edgelists will
%                          be formed from the subpixel coordinates,
%                          otherwise the the integer pixel coordinates of
%                          points in 'im' are used.
%
% Returns:  edgelist - a cell array of edge lists in row,column coords in
%                      the form
%                     { [r1 c1   [r1 c1   etc }
%                        r2 c2    ...
%                        ...
%                        rN cN]   ....]   
%
%           edgeim   - Image with pixels labeled with edge number. Note that
%                      this image also includes edges that do not meet the
%                      minimum length specification.  If you want to see just
%                      the edges that meet the specification you should pass
%                      the edgelist to DRAWEDGELIST.
%
%
% This function links edge points together into lists of coordinate pairs.
% Where an edge junction is encountered the list is terminated and a separate
% list is generated for each of the branche

%% Get superpixels adjecent to edglets

[imh, imw] = size(wseg);
nseg = double(max(wseg(:)));

nbi = imh*[-1 0 1 1 1 0 -1 -1] + [-1 -1 -1 0 1 1 1 0];

% get superpixel pair for each edglet
ne = numel(edges);
sp = zeros(ne, 2);
for k = 1:ne
    ei = edges{k}(ceil(end/2));
    ex = floor((ei-1)/imh)+1;
    ey = mod(ei-1, imh)+1;    
    if     ey==1, sp(k, 1) = 0;    sp(k, 2) = wseg(ey+1, ex); 
    elseif ey==imh, sp(k, 1) = 0;  sp(k, 2) = wseg(ey-1, ex);  
    elseif ex==1, sp(k, 1) = 0;    sp(k, 2) = wseg(ey, ex+1);  
    elseif ex==imw, sp(k, 1) = 0;  sp(k, 2) = wseg(ey, ex-1); 
    else % not on border
        nbvals = wseg(ei+nbi);  
        sp(k, :) = unique(nbvals(nbvals>0));
    end
end


%% Get edgelet/junction adjacency

nj = numel(junctions);
% get mean (x,y) position of junctions
jpos = zeros(nj, 2);
for k = 1:nj
    [jy, jx] = ind2sub([imh imw], junctions{k});
    jpos(k, :) = [mean(jx) mean(jy)];
end

% get 1) junctions adjacent to each edglet
%     2) junction adjacency in form of [edglet nextJunction]
ejunctions = zeros([ne 2]);
jadj = cell(nj, 1);
for k = 1:ne
    ejunctions(k, :) = jim(edges{k}([1 end]));
    jadj{ejunctions(k, 1)}(end+1, :) = [k ejunctions(k, 2)];
    jadj{ejunctions(k, 2)}(end+1, :) = [k ejunctions(k, 1)];
    %存储边的id和变边的结束点id，代表边的方向（j1->j2还是j2->j1）
end

% get directed edglet adjacency
eadj = cell(ne*2, 1);
for k = 1:ne
    % forward edge: assign adjecent edges (+ne if adj edge is reverse)
    j = ejunctions(k, 2); 
    %得到edge的第二个点，通过jadj找到和第二个点相邻的点 
    eadj{k} = setdiff(jadj{j}(:, 1), k); 
    %[eid,j;eid,j]找出jadj中和k相关的（eid，j）相关的的数对，然后过滤边k，得到其他的eid
    reverseind = (ejunctions(eadj{k}, 2)==j); % reverse if meet at 2nd junction
    eadj{k}(reverseind) = ne+eadj{k}(reverseind);
    
    % backward edge: assign adjecent edges (+ne if adj edge is reverse)
    j = ejunctions(k, 1);    
    eadj{k+ne} = setdiff(jadj{j}(:, 1), k); 
    reverseind = (ejunctions(eadj{k+ne}, 2)==j); % reverse if meet at 2nd junction first
    eadj{k+ne}(reverseind) = ne+eadj{k+ne}(reverseind);
end


%% Get edglet orientation from j1 to j2 and order sp so left side is first

imx = repmat((1:imw), [imh 1]);
imy = repmat((1:imh)', [1 imw]);
etheta2 = zeros(ne, 1); % directed orientation
etheta = zeros(ne, 1); % undirected orientation
BIM = false(size(eim));
BIM([1 end], :) = true;  BIM(:, [1 end]) = true;
for k = 1:ne    
    
    % get orientation
    jx = jpos(ejunctions(k, :), 1); 
    jy = jpos(ejunctions(k, :), 2); 
    etheta2(k) = atan2(-(jy(2)-jy(1)), jx(2)-jx(1));
    etheta(k) = mod(etheta2(k) + pi/2, pi)-pi/2;
    if etheta(k)==-pi/2, etheta(k) = pi/2; end
    
    % orient so that left side is first, assuming undirected orientation
    ei = edges{k}(ceil(end/2));    
    ex = imx(ei);    ey = imy(ei);    
    if ey==1, % do nothing as 0 is on left and listed first         
    elseif ey==imh, sp(k, :) = sp(k, [2 1]); % reverse so sp is on left  
    elseif ex==1, % do nothing as 0 is on left and listed first   
    elseif ex==imw, sp(k, :) = sp(k, [2 1]); % reverse so sp is on left
    else % not on border
        stepx = round(cos(etheta(k)+pi/2));
        stepy = -round(sin(etheta(k)+pi/2));
        tmpe = edges{k};
        tmpe = tmpe(~BIM(tmpe));
        leftvals = double(wseg(tmpe + stepy + imh*stepx));  
        leftvals = intersect(leftvals, sp(k, :));
        [mcleft, cnt, allm] = mode(leftvals(leftvals>0));
        if ~any(mcleft==sp(k, :)) || numel(allm)>1 % no unique valid mode: use centers
            ind = (wseg==sp(k, 1));
            x1 = mean(imx(ind));  y1 = mean(imy(ind));
            ind = (wseg==sp(k, 2));
            x2 = mean(imx(ind));  y2 = mean(imy(ind));
            tmptheta1 = atan2(-(y1 -ey), (x1 - ex));
            tmptheta2 = atan2(-(y2 - ey), (x2 - ex));
            [tmpval, lind] = min(abs([tmptheta1 tmptheta2]-(etheta(k)+pi/2)));
            if lind==2, sp(k, :) = sp(k, [2 1]); end            
        elseif sp(k, 1) == mcleft, % do nothing
        elseif sp(k, 2) == mcleft, sp(k, :) = sp(k, [2 1]); % reverse
        else error('invalid case');
        end
    end    

    % if directed orientation is outside of -pi/2 to pi/2, reverse order
    if (etheta2(k) <= -pi/2) || (etheta2(k) > pi/2)
        sp(k, :) = sp(k, [2 1]);
    end
end

                  

%% Store data in bndinfo

bndinfo.wseg2 = uint16(wseg);
bndinfo.nseg = nseg;
bndinfo.ne = ne;
bndinfo.nj = nj;

bndinfo.edges.indices = edges;
bndinfo.edges.adjacency = eadj;
bndinfo.edges.junctions = ejunctions;
bndinfo.edges.spLR = sp; % superpixels on left/right side for j1 --> j2
bndinfo.edges.thetaDirected = etheta2;
bndinfo.edges.thetaUndirected = etheta;

bndinfo.junctions.indices = junctions;
bndinfo.junctions.adjacency = jadj;
bndinfo.junctions.position = jpos;





