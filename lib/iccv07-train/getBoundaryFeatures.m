function [X, Y] = getBoundaryFeatures(bndinfo, pbim, gconf)
% [X, Y] = getBoundaryFeatures(bndinfo, pbim, gconf)
% 
% Computes features for each edgelet based on boundary and geometry
% confidence images.
%
% Input:
%   bndinfo - structure of superpixel boundaries
%   pbim(imh, imw, norient) - probability of boundary image at each orientation
%   gconf(imh, imw, ngeoms) - confidence in each geometric label
%
% Output:
%   X(nboundaries, :) - features for the boundaries
%   Y(nboundaries, 2) - geometry on each size or 0 for no edge

edges = bndinfo.edges;
nbnd = bndinfo.ne*2; % number of directed edges

X = zeros(nbnd, 22);
Y = zeros(nbnd, 2);

% make gconf(:, :, [gnd planar porous solid sky])
gconf = cat(3, gconf(:, :, 1), sum(gconf(:, :, 2:4), 3), gconf(:, :, 5:7));
[imh, imw, ng] = size(gconf);

% compute average confidence (within half-disks at varying orientations) and
% difference confidence images
radius = 6; % note that this will be computed at quarter size
diskfil = double(getnhood(strel('disk', radius, 0)));
diskfil(radius+1:end, :) = 0; % make only upper half of disk ones


% create rotated versions of half-disk
nrot = 8;
diskfil = cat(3, diskfil, zeros([size(diskfil) nrot-1], 'double'));
for k = 2:nrot
    diskfil(:, :, k) = imrotate(diskfil(:, :, 1), (k-1)*45, 'bilinear', 'crop');
end
diskfil = diskfil ./ repmat(sum(sum(diskfil, 2), 1), [radius*2+1 radius*2+1 1]);


% compute average confidence over disk for different rotations
gcabs = zeros([imh imw ng nrot]);
gconf = imresize(single(gconf), 0.25, 'bilinear');
for k = 1:nrot
    gcabs(:, :, :, k) = imresize(imfilter(gconf, diskfil(:, :, k)), [imh imw]);    
end

% compute differences in confidence
gcdif = gcabs(:, :, :, 1:nrot/2) - gcabs(:, :, :, 1+nrot/2:nrot);

% get discretize orientation into 1=1|2, 2=1/2, 3=1_2, 4=2\1
theta = bndinfo.edges.thetaUndirected;
rels = (theta < -3*pi/8) + (theta < -pi/8) + (theta < pi/8) +  (theta < 3*pi/8); 
rels = mod(rels, 4) + 1;
rels = repmat(rels, [2 1]); % double to have directed edges

% gcabs(imh, imw, gtype, [T TL L BL B BR R TR])(Top, Bottom, Left, Right)
% gcdif(imh, imw, gtype, [ -, /, |, \ ]) (direction of edge)
% pbim(imh, imw, [ -, \, |, / ]) (direction of edge)
% rels(k): 1|2, 1/2, 1_2, 2\1 (relationship)

% map relation label in elabels to corresponding matric in gc and pb
gcabsmap1 = [3 2 1 8];
gcabsmap2 = [7 6 5 4];
gcdifmap1 = [3 2 1 -4];
gcdifmap2 = -gcdifmap1;
pbmap = [3 4 1 2];

% compute features
for k = 1:nbnd
    rel = elabels(k, 3);
    eind = edges.indices{mod(k, bndinfo.ne)};
    subi = abs([gcabsmap1(rel) gcdifmap1(rel) gcabsmap2(rel) gcdifmap2(rel) pbmap(rel)]);
    
    if ~isempty(eind)
        for y = 1:ng
            ind = eind + (y-1)*imh*imw + (subi(1)-1)*ng*imh*imw;
            X(k, y) = mean(gcabs(ind));
            ind = eind + (y-1)*imh*imw + (subi(2)-1)*ng*imh*imw;
            X(k, y+ng) = mean(gcdif(ind)) * sign(gcdifmap1(rel));
            ind = eind + (y-1)*imh*imw + (subi(3)-1)*ng*imh*imw;
            X(k, y+2*ng) = mean(gcabs(ind));
            ind = eind + (y-1)*imh*imw + (subi(4)-1)*ng*imh*imw;
            X(k, y+3*ng) = mean(gcdif(ind)) * sign(gcdifmap2(rel));        
        end
        ind = eind + (subi(5)-1)*imh*imw;
        X(k, 4*ng+1) = mean(pbim(ind));
        X(k, 4*ng+2) = rel;
        Y(k, :) = elabels(k, 1:2); 
    else
        Y(k, :) = -1;
    end
end
   
