function orient = getEdgeOrientation(juncts, imsize, edgepix)
% Returns rough orientation of an edge, given the junction points (if available) 
% and the edge points

[jy, jx] = ind2sub(imsize(1:2), juncts);
jy = double(jy);  jx = double(jx);

n = numel(jy);
if n==0 % use two furthest edge points
    [ey, ex] = ind2sub(imsize(1:2), edgepix);  
    distmat = squareform(pdist(double([ex(:) ey(:)])));
    [mval, mind] = max(distmat(:));
    [ind(1), ind(2)] = ind2sub(size(distmat), mind);
    jx = double(ex(ind));
    jy = double(ey(ind));            
elseif n==1 % use furthest edge point as other junction
    [ey, ex] = ind2sub(imsize(1:2), edgepix);   
    [maxval, maxind] = max((jy-double(ey)).^2 + (jx-double(ex)).^2);
    jy(2) = double(ey(maxind));  jx(2) = double(ex(maxind));    
elseif n>2 % find furthest two junction points
    distmat = squareform(pdist([jx(:) jy(:)]));
    [mval, mind] = max(distmat(:));
    [ind(1), ind(2)] = ind2sub(size(distmat), mind);
    jx = jx(ind);
    jy = jy(ind);
end

if jx(2)==jx(1)
    orient = pi/2;
else
    orient = atan(-(jy(2)-jy(1)) ./ (jx(2)-jx(1))); % remember y2>y1 means y2 is lower
end