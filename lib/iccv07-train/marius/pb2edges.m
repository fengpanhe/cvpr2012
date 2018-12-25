function [eim, ecount, nedge] = pb2edges(bndinfo, X, grayim)

eim = zeros(bndinfo.imsize(1:2));
ecount = [];
nedge = 0;
for j = 1:size(X.edge.pbOrient, 2),     
    pbim = zeros(size(eim));
    for k = 1:numel(bndinfo.edges.indices)
        ind = bndinfo.edges.indices{k};
        pbim(ind) = X.edge.pbOrient(k, j);
    end
    tmpim = getCannyEdgesFast_Dwh(grayim, pbim);  
    ind = tmpim>0;
    ntmp = max(tmpim(ind));
    ecount(nedge+1:nedge+ntmp) = hist(tmpim(ind), (1:ntmp));
    eim(ind) = tmpim(ind) + nedge;    
    nedge = nedge + ntmp;
end