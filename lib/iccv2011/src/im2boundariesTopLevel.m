function [bndinfo, pbim, gconf, bndinfo_all] = im2boundariesTopLevel(im, thresh)
%
% [bndinfo, pbim, gconf, bndinfo_all] = im2boundariesTopLevel(im, thresh)
%
% High level function for estimating boundaries in an image using Hoiem et
% al. 2007 occlusion reasonining method.  
%
% Inputs:
%   im:       RGB double format image
%   thresh:   threshold values for segmentation hierarchies (default =
%             [0.105 0.25 0.6], threshold is lowest boundary confidence 
%             required not to merge to regions in hierarchical segmentation
%
% Outputs:
%   bndinfo.wseg: the superpixels
%   bndinfo.edges.indices: the indices for each edge
%   bndinfo.result.edgeProb: probability of edge being on
%   bndinfo.result.geomProb: probabilty of geometric label (gnd, vert, sky)
%   bndinfo.result.boundaries: most likely ege label (0=off, 1=on)
%   pbim: probability of boundary (Pb) image (imh, imw, 4 orient) 
%   gconf: surface likelihoods (imh, imw, [support, vert-planar/L/C/R,
%          non-planar solid/porous])
%   bndinfo_all: boundary results after each stage
%
% Notes:
%   - length of result.boundaries is twice length of edges.indices.  If
%   label in first half is "on", left side occludes; if label in second
%   half is "on", right side occludes.   


if ~exist('thresh', 'var') || isempty(thresh)
    thresh = [0.105 0.25 0.6];
end

% load classifiers
datadir = './data/';

load(fullfile(datadir, 'boundaryClassifiers.mat'));
load(fullfile(datadir, 'continuityClassifiers.mat'));
gclassifiers1 = load(fullfile(datadir, 'ijcvClassifier.mat'));
gclassifiers2 = load(fullfile(datadir, 'perfectSegClassifierCv.mat'));

% create pb confidences
disp('pb')
pbim = pbCGTG_nonmax(im);

% create geometric confidences
disp('geometry')
[pg, ~, imsegs] = ijcvTestImage(im, [], gclassifiers1);
gdata.imsegs = imsegs;
gconf = pg2confidenceImages(imsegs, {pg});
gconf = gconf{1}(:, :, 1:7);

% get occlusion boundary labels
disp('boundaries')
[bndinfo, bndinfo_all] = im2boundaries(im, pbim, gconf, dtBnd, dtBnd_fast, dtCont, ...
    gdata, gclassifiers2, thresh);
