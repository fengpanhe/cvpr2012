function enums = getObjectBoundary(bndinfo, objnum)

elabs = bndinfo.labels(bndinfo.adjlist);

% exactly one side of edge is in object
ind = (sum(elabs==objnum, 2)==1) & ~any(elabs==0, 2); 

enums = find(ind);

keep = true(size(enums));
for k = 1:numel(enums)
    if numel(bndinfo.edges{enums(k)}) < 2
        keep(k) = false;
    end
end
enums = enums(keep);

% get ordering of edglets, starting from lowest and continuing
% counter-clockwise
ordering = zeros(numel(enums), 1);
maxy = zeros(numel(enums), 1); % max y is lowest

for k = 1:numel(enums)
    maxy(k) = max(mod(bndinfo.ejunctions{enums(k)}-1, bndinfo.imsize(1))+1);
end

[tmpval, ordering(1)] = max(maxy);

for k = 2:numel(ordering)
    for j = 1:numel(enums)
        e1 = min(enums(k-1), enums(j));
        e2 = max(enums(k-1), enums(j));
        if bndinfo.edgeAdjacency(e1, e2)
            j1 = intersect(bndinfo.ejunctions{e1}, bndinfo.ejunctions{e2});
            lside = bndinfo.labels(getEdgeLeftRightSides(bndinfo, enums(j), j1));
            if lside==objnum
                ordering(k) = j;
                break;
            end
        end
    end
end



enums = enums(ordering(ordering>0));

tmpim = zeros(bndinfo.imsize(1:2));
for k = 1:numel(enums)
    tmpim(bndinfo.edges{enums(k)}) = k;
end
figure(1), hold off, imagesc(tmpim), axis image, colormap gray
keyboard