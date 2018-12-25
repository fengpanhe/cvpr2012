function [leftS, rightS] = getEdgeLeftRightSides(bndinfo, enum, j1)
% Gets superpixel index of left and right side of edge with index enum and
% starting junctions j1

imh = bndinfo.imsize(1);

j2 = setdiff(bndinfo.ejunctions{enum}, j1);


jx1 = mean(floor((j1-1)/imh)+1);
jy1 = mean(mod(j1-1, imh)+1);

if isempty(j2)
    e = bndinfo.edges{enum};
    ex = floor((e-1)/imh)+1;
    ey = mod(e-1, imh)+1;
    [tmp, ind] = max((ex-jx1).^2 + (ey-jy1).^2);
    j2 = e(ind);
end
jx2 = mean(floor((j2-1)/imh)+1);
jy2 = mean(mod(j1-2, imh)+1);

theta = atan2(jy2-jy1, jx2-jx1);

% with adjlist(enum, 1:2) = [A B]: "A left Of B", "A above B", 
% "A above-aeft of B", or "A above-right of B" for undirected edge
if (theta > -pi/2+pi/8) && (theta < pi/2 + pi/8)
    leftS = bndinfo.adjlist(enum, 1);
    rightS = bndinfo.adjlist(enum, 2);
else
    leftS = bndinfo.adjlist(enum, 2);
    rightS = bndinfo.adjlist(enum, 1);    
end