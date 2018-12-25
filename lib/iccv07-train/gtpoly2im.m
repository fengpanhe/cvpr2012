function gtpoly2im(im, imsegs, pball, bgt)

segimage = imsegs.segimage;

pb = max(pball, [], 3);
% compute superpixel image
wseg = watershed(pb.*(pb>0.025));
stats = regionprops(wseg, 'PixelIdxList');

% break image into ground, sky, and objects
oim = zeros(imsegs.imsize(1:2));

oim(imsegs.vert_labels(segimage)==1) = 1;
oim(imsegs.vert_labels(segimage)==3) = 2;

c = 2;
for k = 1:numel(imsegs.horz_names)
    if any(imsegs.horz_labels==k)
        oim(imsegs.horz_labels(segimage)==k) = c + 1;
        c = c+1;
    end
end

oim2 = zeros(size(oim));
pixind = {stats(:).PixelIdxList};
for k = 1:numel(pixind)
   oim2(pixind{k}) = mode(oim(pixind{k}));
end

figure(1), hold off, imagesc(im), axis image
figure(2), hold off, imagesc(oim), axis image, colormap jet
figure(3), hold off, imagesc(oim2), axis image, colormap jet

cim = zeros(size(oim));
c = 0;
for k = 1:numel(bgt.vpoly)
    for j = 2:size(bgt.vpoly{k}, 2)
        c = c+1;
        pts = bgt.vpoly{k}(1:2, j-1:j);
        [xs, ys] = line2pixels(pts(1), pts(2), pts(3), pts(4));
        cim(ys + (xs-1)*size(oim2, 1)) = c;
    end
end
figure(4), hold off, imagesc(cim), axis image, colormap jet
mean(pb(cim>0))
mean(pb(:))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xs, ys] = line2pixels(x1,y1,x2,y2)
% given (x1,y1), (x2,y2), compute each pixel along the way

dy = y2-y1;
dx = x2-x1;
if abs(dy) > abs(dx)
    step = (y2>y1)*2-1;
    ys = (y1:step:y2);
    xs = x1 + (ys-y1)/dy*dx;
else
    step = (x2>x1)*2-1;
    xs = (x1:step:x2);
    ys = y1 + (xs-x1)/dx*dy;
end
xs = round(xs);
ys = round(ys);