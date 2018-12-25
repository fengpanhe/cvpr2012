function bndinfo2 = geometry2seg(bndinfo, gcdir)

for f = 1:numel(bndinfo)   
    
    bn = strtok(bndinfo(f).imname, '.');
        
    stats = regionprops(bndinfo(f).wseg, 'Area', 'PixelIdxList');
    idx = {stats.PixelIdxList};
    area = vertcat(stats.Area);
    nsp = bndinfo(f).nseg;
    imh = bndinfo(f).imsize(1);  imw = bndinfo(f).imsize(2);
    
    tmp = load([gcdir bn '.c.mat']);    
    gconf = tmp.cim_geom_7;               

    splab = zeros(nsp, 1);
    
    % get most likely label for each segment
    pg = zeros(bndinfo(f).nseg, 5);    
    
    gconf = cat(3, gconf(:, :, 1), sum(gconf(:, :, 2:4), 3), gconf(:, :, 5:7));    
    gconf = reshape(gconf, [imh*imw 5]);
    for g = 1:5
        gim = gconf(:, g);
        for k = 1:nsp    
            pg(k,g) = sum(gim(idx{k})) / area(k);
        end
    end        
    gvs = [pg(:, 1) sum(pg(:, 2:4), 2) pg(:, 5)];    
    pps = pg(:, 2:4);
    
    [tmp, glab] = max(gvs, [], 2);
    glab(glab==3) = 5;
    ind = glab==2;
    [tmpval, tmpmax] = max(pps, [], 2);
    glab(ind) = 1+tmpmax(ind);

    gim = glab(bndinfo(f).wseg);
    c = 0;
    for g = 1:5
        bwim = bwlabel(gim==g);
        for k = 1:max(bwim(:))
            c = c + 1;
            splab(unique(bndinfo(f).wseg(bwim==k))) = c; 
        end
    end
    
    bndinfo2(f) = updateBoundaryInfo2(bndinfo, splab);
    
end  