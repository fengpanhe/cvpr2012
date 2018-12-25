function [seg, score] = getPopout(resdir, imdir, fn, lab, nseg)


co = 0;

imnum = zeros(10000,1);
objnum = zeros(10000,1);
score = zeros(10000,1);

for f = 1:numel(fn)
    
    load([resdir strtok(fn{f}, '.') '_seg.mat']);        
    
    for k = 1:bndinfo.nseg
        co = co + 1;
        imnum(co) = f;
        objnum(co) = k;
        score(co) = bndinfo.result.geomProb(k, lab);
        eprob = bndinfo.result.edgeProb;
        ind = find(bndinfo.edges.spLR(:, 1)==k); 
        score(co) = score(co)*prod(eprob(ind, 2));
        ind = find(bndinfo.edges.spLR(:, 2)==k);
        score(co) = score(co)*prod(eprob(ind, 3));
    end
end

[score, sind] = sort(score(1:co), 'descend');
imnum = imnum(sind);
objnum = objnum(sind);
score = score(sind);


for k = 1:nseg
    load([resdir strtok(fn{imnum(k)}, '.') '_seg.mat']);
    seg{k} = zeros([bndinfo.imsize(1:2) 3]);
    seg{k}(:, :, 3) = 1;
    im = im2double(imread([imdir fn{imnum(k)}]));
    
    pix = find(bndinfo.wseg==objnum(k));
    [py, px] = ind2sub(bndinfo.imsize(1:2), pix);
    
    npix = bndinfo.imsize(1)*bndinfo.imsize(2);
    
    for b = 1:3
        seg{k}(pix+(b-1)*npix) = im(pix+(b-1)*npix);
    end
    
    seg{k} = seg{k}(min(py):max(py), min(px):max(px), :);
end
    
score = score(1:nseg);