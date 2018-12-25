function findOccludedLines(im)

[imh, imw] = size(im);

theta = [0:179];
[R, xp] = radon(im, theta);

% get the best line from the hough transform image 
[maxval, ind] = max(R(:));    

% get line parameters (ax + by = c)
[rind,thind] = ind2sub(size(R),ind);
t = -theta(thind)*pi/180;
r = xp(rind); 
lines = [cos(t) sin(t) -r];
cx = imw/2-1;
cy = imh/2-1;
lines(:,3) = lines(:,3) - lines(:,1)*cx - lines(:,2)*cy;  
a = lines(1, 1);
b2 = lines(1, 2);
c = lines(1, 3);


% get line parameters in y = mx+b form
m = -lines(:, 1)/lines(:, 2);
b = -lines(:, 3)/lines(:, 2);         

% get the distance of each point from the line
pdst = abs(a*tx+b2*ty+c)/sqrt(a^2+b2^2);
inds = find(pdst < distt);