function [X, Y] = getAllBoundaryFeatures(bndinfo, pbdir, gcdir)

for f = 1:numel(bndinfo)

    bn = strtok(bndinfo(f).imname, '.');
    tmp = load([pbdir bn '_pb.mat']);
    pball = tmp.pball;
    tmp = load([gcdir bn '.c.mat']);
    
    gconf = tmp.cim_geom_7;
%    gconf = gconf ./ repmat(sum(gconf, 3) , [1 1 size(gconf, 3)]);    
    
    [X{f}, Y{f}] = getBoundaryFeatures3(bndinfo(f), pball, gconf);
end

