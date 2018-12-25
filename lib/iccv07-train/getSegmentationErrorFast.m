function err = getSegmentationErrorFast(labels, nsp, rsp, area)
% err = getSegmentationError(seg1, seg2, direction)
%
% Returns the oversegmentation error of seg2 with respect to seg1.  

nlab = max(labels);

nr = numel(rsp);

npix = sum((labels>0).*area);

total = 0;

count = zeros(nr,1);
for k = 1:nr
    count(k) = numel(rsp{k});
end
for k = find(count==1)'
    total = total + (labels(rsp{k})>0)*area(rsp{k});
end

for k = find(count>1)'
    rcount = zeros(nlab, 1);
    labs = labels(rsp{k});
    valid = find(labs>0);        
    for k2 = valid(:)'
        rcount(labs(k2)) = rcount(labs(k2)) + area(rsp{k}(k2)); 
    end
    total = total + max(rcount);
end

err = (npix-total) / npix;



    