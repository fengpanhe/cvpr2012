%
%    [cur_locs, dist, f] = getNeighbors3(locs2, yx, idx, rad)
%
%    Find edge points in an image which are within a radius rad of a center point.
%
%  Input:
%    locs2: A set of input (x,y) locations
%    yx: The input reference point, we want to find all the points in locs2 within a radius rad of yx
%    idx:  An array of the same size as the underlying image, with idx(i,j) == 1 if there is an edge
%          point at position (i,j). Idx is computed once so that nearest neighbors can be retrieved 
%          efficiently.
%    rad: Radius of the neighborhood
%
% Output:
%    cur_locs: Locations of the edge points within a redius "rad" of yx
%    dist: Distance of the neighbors to the reference point
%    f: Index of the neighbors in the image (idx)
%
% Comments:
%
%
% Exceptions:
%   None
%


function  [cur_locs, dist, f] = getNeighbors3(locs2, yx, idx, rad)

%disp(' getNeighbors.m ');

yxl = round(yx-rad);

if yxl(1) <= 0
  yxl(1) = 1;
end

if yxl(2) <= 0
  yxl(2) = 1;
end

yxu = round(yx+rad);

if yxu(1) > size(idx,1)
     yxu(1) = size(idx,1);
end

if yxu(2) > size(idx,2)
     yxu(2) = size(idx,2);
end

idxs = idx(yxl(1):yxu(1),yxl(2):yxu(2));

ff = find(idxs);

candidates = idxs(ff);

candidates = candidates(:);

rad = rad^2;

yx = repmat(yx,length(ff),1);

locs2(candidates,:);

dist = sum(((yx - locs2(candidates,:)).^2)');

ff = find(dist < rad);

f = candidates(ff);

cur_locs = locs2(f,:);
dist = sqrt(dist(ff));

