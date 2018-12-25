function X = updateGeometricFeatures2(fn, X, segmaps, cvnum, gdatadir)

gclassdir = '/home/dhoiem/cmu/GeometricContext/data/ijcv06/';
%gclass1 = load([gclassdir 'multisegClassifierCv.mat']);
gclass2 = load([gclassdir 'perfectSegClassifierCv.mat']);

% gclass1.vclassifier = gclass1.vclassifier(cvnum);
% gclass1.hclassifier = gclass1.hclassifier(cvnum);
gclass2.vclassifier = gclass2.vclassifier(cvnum);
gclass2.hclassifier = gclass2.hclassifier(cvnum);

for f = 1:numel(fn)
    bn = strtok(fn{f}, '.');
    
    gdata = load([gdatadir bn '_gdata.mat']);
    
    if isempty(segmaps{f}), segmaps{f} = (1:gdata.imsegs.nseg); end
%    disp('changed geometry')    
%     pg = ijcvEstimateGeometry([], gdata.imsegs, gclass1, ...
%         segmaps{f}, gdata.spdata, gdata.imdata);
%     
%     X(f).region.pg2 = pg;    
    
    pg = ijcvEstimateGeometry([], gdata.imsegs, gclass2, ...
        segmaps{f}, gdata.spdata, gdata.imdata);
    
    X(f).region.pg2 = pg;
    
%     area = gdata.imsegs.npixels;
%     X(f).region.textureHist = zeros(max(segmaps{f}(:)), 15);
%     for k = 1:max(segmaps{f})
%         ind = find(segmaps{f}==k);
%         X(f).region.textureHist(k, :) = sum(gdata.spdata(ind, 30:44) .* ...
%             repmat(area(ind), [1 15])) / sum(area(ind));
%     end        
    
end