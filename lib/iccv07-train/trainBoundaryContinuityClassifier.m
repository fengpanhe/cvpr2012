function classifier = trainBoundaryContinuityClassifier(...
    bndinfo, X, Y, pC, bndfeatures)

[tx, ty, categoryFeatures, cn, in] = reformatData(bndinfo, X, Y, pC, bndfeatures);    

classifier = train_boosted_dt_2c(tx, categoryFeatures, ty*2-1, 20, 8);

py = test_boosted_dt_mc(classifier, tx);
py = 1./(1+exp(-py));


pall = cell(100000, 1);
yall = cell(100000, 1);
c = 0;
lastc = 0;  lasti = 0;
for k = 1:numel(ty)
    if cn(k)~=lastc || in(k)~=lasti
        lastc = cn(k);  lasti = in(k);
        c = c + 1;
    end
    pall{c}(end+1) = py(k);
    yall{c}(end+1) = ty(k);
end

iscorrect = zeros(c, 1);
deltap = zeros(c,1);
cempty = zeros(c, 1);
hasmultiple = zeros(c,1);
for k = 1:c
       
    pvals = pall{k};
    cind = find(yall{k});    
    
    if isempty(cind)
        cempty(k) = 1;
    else
    
        if numel(cind)>1
            hasmultiple(k) = 1;
        end
        
        pvals = pvals / sum(pvals);
        cval = pvals(cind);
        
        [svals, sind] = sort(pvals, 'descend');
        
        [cval, tind] = max(cval);
        cind = cind(tind);
        
        iscorrect(k) = (cind(1)==sind(1));

        if iscorrect(k) && numel(svals)>1
            deltap(k) = cval - svals(2);
        else
            deltap(k) = cval - svals(1);
        end
        
    end
end
    
disp(num2str(mean(iscorrect)));
disp(num2str(mean(deltap)));




%% Reformat the data for classifying between boundary and no boundary
function [tx, ty, cf, cnum, inum] = reformatData(bndinfo, X, Y, pC, bndx)

nc = 0; % number of contours
ndata = 0; % number of possible directed transitions from a contour
for f = 1:numel(Y)
    ind = find(Y{f} > 0);
    nc = nc + numel(ind);
    for k = 1:numel(ind)        
        ndata = ndata + numel(bndinfo(f).edges.adjacency{ind(k)});        
    end    
end

disp(num2str([nc ndata]));

ty = zeros(ndata, 1); % whether it is a continuation of the same contour
tx = zeros(ndata, 126); % features: continuity + pContour of same type

cnum = zeros(ndata, 1); % record edge index of contour
inum = zeros(ndata, 1); % record image index of contour

c = 0;
for f = 1:numel(Y)
    
    ind = find(Y{f} > 0);
    spLR = bndinfo(f).edges.spLR;
    spLR = [spLR ; spLR(:, [2 1])];
    oleft = bndinfo(f).labels(spLR(ind, 1));
       
    for k = 1:numel(ind)
        
        adj = bndinfo(f).edges.adjacency{ind(k)};
        nadj = numel(adj);
        
        type = Y{f}(ind(k));
        [tx(c+(1:nadj), :), cf] = ...
            getContinuityFeatures4(X(f), bndinfo(f), ind(k), adj, pC{f}, bndx{f});
        
        cnum(c+(1:nadj)) = ind(k);
        inum(c+(1:nadj)) = f;
        
        tind = find(Y{f}(adj)>0);
        oleft2 = bndinfo(f).labels(spLR(adj(tind), 1));
        ty(c+tind) = (oleft(k) == oleft2);                

        c = c + nadj;
        
    end

end


