function testImage(bndinfo, dtBnd, thresh, imdir, pbdir, gcdir, cnum)

bndinfo1 = bndinfo;

thresh(end) = 0;
for s = 1:numel(dtBnd)

    X = getAllFeatures(bndinfo1, imdir, pbdir, gcdir);
        
    [result(s), dispim] = mergeStage(X, bndinfo1, dtBnd(s), thresh(s));                  
    [bndinfo1, terr(s)] = updateBoundaryInfo2(bndinfo1, result(s));   
    disp(['Pixel error: ' num2str(terr(s))])   
   
end

%% display final results

regions = result(end).regions;
for s = numel(result)-1:-1:1
    prevregions = result(s).regions;    
    regions2 = cell(size(regions));
    for k = 1:numel(regions)
        for k2 = regions{k}
            regions2{k} = [regions2{k} prevregions{k2}];
        end
    end
    regions = regions2;
end

final.regions = regions;
[tmp, terr] = updateBoundaryInfo2(bndinfo, final);   
disp(['Pixel error: ' num2str(terr)])  

X = getAllFeatures(bndinfo1, imdir, pbdir, gcdir);
splab = zeros(bndinfo.nseg, 1);
for k = 1:numel(regions)
    splab(regions{k}) = k;
end
pB = useBoundaryClassifier(bndinfo1, X, dtBnd(end), 'c');

disp(num2str(sum(pB(:, 1)>0.5)))

[pbval, pbind] = max(pB(:, 2:3), [], 2);

pG = X.region.geomContext;

pgvs = [pG(:, 1) sum(pG(:, 2:4), 2) pG(:, 5)];
[pval, gvslab] = max(pgvs, [], 2);
ind = (gvslab==2);

[pval2, sublab] = max(pG(ind, 2:4), [], 2);

types(gvslab==1) = 1;
types(gvslab==3) = 2;
types(gvslab==2) = sublab+2;

newtypes = zeros(numel(splab));
newtypes(splab>0) = types(splab(splab>0));

spLR = bndinfo.edges.spLR;
spLR = [spLR ; spLR(:, [2 1])];
spLR2 = bndinfo1.edges.spLR;
%spLR2 = [spLR2 ; spLR2(:, [2 1])];
btype = zeros(bndinfo.ne*2, 1);
for k = 1:bndinfo.ne*2
    r1 = splab(spLR(k, 1));
    r2 = splab(spLR(k, 2));
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
lim = displayObjectLabels(im, newtypes, bndinfo.wseg);
lim = drawBoundaries(lim, bndinfo, btype);
figure(2), imagesc(lim), axis image

imwrite(dispim, ['./results/initial/' bn '_segs.jpg'], 'Quality', 90);

imwrite(lim, ['./results/initial/' bn '_result.jpg'], 'Quality', 90);
imwrite(im, ['./results/initial/' bndinfo.imname]);

lim2 = displayObjectLabels(im, bndinfo.labels, bndinfo.wseg);
lim2 = drawBoundaries(lim2, bndinfo, bndinfo.edges.boundaryType);
imwrite(lim2, ['./results/initial/' bn '_gt.jpg'], 'Quality', 90);

tmp = load([gcdir bn '.c.mat']);
gconf = tmp.cim_geom_7;
pv = cat(3, gconf(:, :, 1), sum(gconf(:, :, 2:6), 3), gconf(:, :, 7));
ph = gconf(:, :, 2:6) ./ repmat(sum(gconf(:, :, 2:6), 3), [1 1 5]);
limg = APPgetLabeledImage2(im, [], pv, ph);
imwrite(limg, ['./results/initial/' bn '_geom.jpg'], 'Quality', 90);

tmp = load('/usr1/projects/GeometricContext/data/ijcv06/multisegResults2.mat');
classifiers.vclassifier = tmp.vclassifier(cnum);
classifiers.hclassifier = tmp.hclassifier(cnum);

%% Recompute geometry
%classifiers = load('/usr1/projects/GeometricContext/data/ijcv06/ijcvClassifier');
imsegs.segimage = bndinfo.wseg;
imsegs.nseg = bndinfo.nseg;
stats = regionprops(imsegs.segimage, 'Area');
imsegs.npixels = vertcat(stats(:).Area);
imsegs.imsize = bndinfo.imsize;

[pg, data] = ijcvEstimateGeometry(im, imsegs, classifiers, splab);

pg = pg(splab, :);

[pv, ph] = splitpg(pg);
limg2 = APPgetLabeledImage2(im, imsegs, pv, ph);
imwrite(limg2, ['./results/initial/' bn '_geom2.jpg'], 'Quality', 90);

