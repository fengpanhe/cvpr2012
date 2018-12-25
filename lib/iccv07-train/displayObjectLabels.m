function lim = displayObjectLabels(im, splab, wseg)

if size(im, 3)==3
    grayim = repmat(rgb2gray(im), [1 1 3]);
else
    grayim = repmat(im, [1 1 3]);
end

tmplim = splab(wseg); 

nlab = max(splab);
colors = hsv2rgb([0.25 0.55 (1:10)/10]', ones(12,1), ones(12,1))*0.5;
lim = im;
for k = 1:nlab
    ind = (tmplim==k);
    nind = sum(ind(:)); 
    cols = colors(mod(k-1, size(colors, 1))+1, :);
    lim(repmat(ind, [1 1 3])) = reshape(repmat(cols, [nind 1]), [nind*3 1]) + ...
        grayim(repmat(ind, [1 1 3]))*0.5;    
end