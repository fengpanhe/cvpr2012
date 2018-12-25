function [pbim, pcim, ptim] = displayContourResults(im, bndinfo, pb, pc, pt)
% display pB (prob boundary) pC (prob contour)  pT ( prob type)

grayim = im2double(rgb2gray(im));

pbim = zeros(size(im));
pcim = zeros(size(im));
ptim = zeros(size(im));

[imh, imw, tmp] = size(im);

gind = floor(((1:36)-1)/6)+1 == 1 | mod((1:36)-1, 6)+1==1;
sind = floor(((1:36)-1)/6)+1 == 2 | mod((1:36)-1, 6)+1==2;
vind = ~gind & ~sind;

for k = 1:numel(pb)

    ind = bndinfo.edges{k};              
    
    pbim(ind) = pb(k);
    pbim(imh*imw+ind) = 0;
    pbim(imh*imw*2+ind) = 0;
    
    pcim(ind) = pc(k);
    pcim(imh*imw+ind) = 0;
    pcim(imh*imw*2+ind) = 0;    
        
    ptk = reshape(pt(k, :), [6 6]);
    
    ptim(0*imh*imw + ind) = pc(k) * sum(sum(ptk(3:6, 3:6))); % r = vert
    ptim(1*imh*imw + ind) = pc(k) * (sum(ptk(:, 1)) + sum(ptk(1, :))); % g = gnd
    ptim(2*imh*imw + ind) = pc(k) * (sum(ptk(:, 2)) + sum(ptk(2, :))); % b = sky
    
end

pbim = widen(pbim / max(pbim(:)), 1);
pcim = widen(pcim / max(pcim(:)), 1);
ptim = widen(ptim / max(ptim(:)), (1:3));

pbim = hsv2rgb(ones([imh imw]), pbim(:, :, 1), grayim*0.5+pbim(:, :, 1)/2);
pcim = hsv2rgb(ones([imh imw]), pcim(:, :, 1), grayim*0.5+pcim(:, :, 1)/2);
[tmph, tmps, tmpv] = rgb2hsv(ptim);
ptim = hsv2rgb(tmph, max(ptim, [], 3), grayim*0.5+max(ptim, [], 3)/2);

%pbim = fillim(widen(pbim, 1), grayim, 0.025);
%pcim = fillim(widen(pcim, 1), grayim, 0.025);
%ptim = fillim(widen(ptim, (1:3)), grayim, 0.01);



disp(num2str(mean(pt.*repmat(pc, [1 36]), 1)))
disp(num2str(max(pt.*repmat(pc, [1 36]), [], 1)))



figure(1), imagesc(im), axis image
figure(2), imagesc(pbim), axis image, colormap gray
figure(3), imagesc(pcim), axis image, colormap gray
figure(4), imagesc(ptim), axis image


%% widen lines
function pim = widen(pim, b)
pim(2:end-1, :, b) = max(pim(2:end-1, :, b), pim(3:end, :, b));
pim(2:end-1, :, b) = max(pim(2:end-1, :, b), pim(1:end-2, :, b));
pim(:, 2:end-1, b) = max(pim(:, 2:end-1, b), pim(:, 3:end, b));
pim(:, 2:end-1, b) = max(pim(:, 2:end-1, b), pim(:, 1:end-2, b));

%% put gray image where values are less than minval
function pim = fillim(pim, grayim, minval)

tmpim = repmat(all(pim<minval, 3), [1 1 3]);
grayim = repmat(grayim, [1 1 3]);
pim(tmpim) = grayim(tmpim);

