function testImage2(bndinfo, dtBnd, thresh, ...
    imdir, pbdir, gcdir, gdatadir, gclassifiers, outdir)

bndinfo1 = bndinfo;
segmap = (1:bndinfo.nseg);

%thresh(end) = 0.75;
for s = 1:numel(dtBnd)
    
    X = getAllFeatures(bndinfo1, imdir, pbdir, gcdir, gdatadir, ...
        gclassifiers(s), {segmap});
        
    [result(s), dispim] = mergeStage(X, bndinfo1, dtBnd(s), thresh(s));                  
    
    [bndinfo1, terr(s)] = updateBoundaryInfo2(bndinfo1, result(s));      
    [segmap, final.regions] = getSegmentationMap(bndinfo.nseg, result);    
     
    disp(['Pixel error: ' num2str(terr(s))])          
    
end

if 1
changed = 1;
count = 0;
while changed
    
    count = count + 1;
    oldnum = bndinfo1.nseg;
    X = getAllFeatures(bndinfo1, imdir, pbdir, gcdir, gdatadir, ...
        gclassifiers(s), {segmap});
        
    [result(s+count), dispim] = mergeStage(X, bndinfo1, dtBnd(s), thresh(s));                  
    
    [bndinfo1, terr(s+count)] = updateBoundaryInfo2(bndinfo1, result(s+count));      
    [segmap, final.regions] = getSegmentationMap(bndinfo.nseg, result);      
    
    newnum = bndinfo1.nseg;
    changed = oldnum~=newnum;
end
end

%% display final results

[tmp, terr] = updateBoundaryInfo2(bndinfo, final);   
disp(['Final pixel error: ' num2str(terr)])  

X = getAllFeatures(bndinfo1, imdir, pbdir, gcdir, gdatadir, ...
    gclassifiers(end), {segmap});
pB = useBoundaryClassifier(bndinfo1, X, dtBnd(end), 'c');

disp(num2str(sum(pB(:, 1)>0.5)))
[pbval, pbind] = max(pB(:, 2:3), [], 2);

pg2 = X.region.geomContext2;
pg2 = [pg2(:, 1) sum(pg2(:, 2:4), 2) pg2(:, 5:7)];

pgvs = [pg2(:, 1) sum(pg2(:, 2:4), 2) pg2(:, 5)];
[pval, gvslab] = max(pgvs, [], 2);
ind = (gvslab==2);

[pval2, sublab] = max(pg2(ind, 2:4), [], 2);

types(gvslab==1) = 1;
types(gvslab==3) = 2;
types(gvslab==2) = sublab+2;

newtypes = zeros(numel(segmap),1);
newtypes(segmap>0) = types(segmap(segmap>0));

spLR = bndinfo.edges.spLR;
spLR = [spLR ; spLR(:, [2 1])];
spLR2 = bndinfo1.edges.spLR;
%spLR2 = [spLR2 ; spLR2(:, [2 1])];
btype = zeros(bndinfo.ne*2, 1);
for k = 1:bndinfo.ne*2
    r1 = segmap(spLR(k, 1));
    r2 = segmap(spLR(k, 2));
    ind = find(spLR2(:,1)==r1 & spLR2(:,2)==r2);
    if ~isempty(ind)
        ind = ind(1);
        if pbind(ind)==1
            lab = newtypes(spLR(k,1));
            btype(k) = (lab==1)*1 + (lab==3)*2 + (lab==4)*3 + (lab==5)*4;
            %btype(k) = max(btype(k), 1);
        end
    else
        ind = find(spLR2(:,2)==r1 & spLR2(:,1)==r2);
        if ~isempty(ind)
            ind = ind(1);
            if pbind(ind)==2
                lab = newtypes(spLR(k,1));
                btype(k) = (lab==1)*1 + (lab==3)*2 + (lab==4)*3 + (lab==5)*4;
                %btype(k) = max(btype(k), 1);
            end
        end
    end
end
        

bn = strtok(bndinfo.imname, '.');

im = imread([imdir bndinfo.imname]);
im = im2double(im);

tmplab = segmap;
tmplab(newtypes==1) = 1;
tmplab(newtypes==2) = 2;
lim = displayObjectLabels(im, tmplab, bndinfo.wseg);
lim = drawBoundaries(lim, bndinfo, btype);
%figure(2), imagesc(lim), axis image
%drawnow;

imwrite(dispim, [outdir bn '_segs.jpg'], 'Quality', 90);

imwrite(lim, [outdir bn '_result.jpg'], 'Quality', 90);
imwrite(im, [outdir bndinfo.imname]);

lim2 = displayObjectLabels(im, bndinfo.labels, bndinfo.wseg);
lim2 = drawBoundaries(lim2, bndinfo, bndinfo.edges.boundaryType);
imwrite(lim2, [outdir bn '_gt.jpg'], 'Quality', 90);

tmp = load([gcdir bn '.c.mat']);
gconf = tmp.cim_geom_7;
pv = cat(3, gconf(:, :, 1), sum(gconf(:, :, 2:6), 3), gconf(:, :, 7));
ph = gconf(:, :, 2:6) ./ repmat(sum(gconf(:, :, 2:6), 3), [1 1 5]);
limg = APPgetLabeledImage2(im, [], pv, ph);
imwrite(limg, [outdir bn '_geom.jpg'], 'Quality', 90);

pg2 = X.region.geomContext2;
pg2 = pg2(segmap, :);
[pv, ph] = splitpg(pg2);

tmp = load([gdatadir bn '_gdata.mat']);

limg2 = APPgetLabeledImage2(im, tmp.imsegs, pv, ph);
imwrite(limg2, [outdir bn '_geom2.jpg'], 'Quality', 90);

