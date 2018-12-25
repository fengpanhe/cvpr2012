function [fragments, neighbors, seg, polyfragments, poly_params] = seg2fragments(seg, img, min_area, order)
%
% [fragments, neighbors, new_seg, <polyfragments>, <poly_params>] = ...
%   seg2fragments(seg, img, min_area, <order>)
% 
%  Finds the boundary fragments for an oversegmentation of an image:
%
%  (1) Takes in an over-segmentation labeling (with labels 1:N),
%  (2) Fills any pixels with a label of zero with the neighboring label 
%      that has the most similar color value in the original RGB color 
%      image 'img', 
%  (3) Cleans up segments which have a single pixel jutting diagonally off
%      of a corner somewhere, by reassigning these offending pixels to the
%      neighboring segment with the most similar color (as in (2)).
%  (4) Merges segments that are too small with a neighbor with the most 
%      similar average color, and 
%  (5) Chains fragments along the boundaries between segments, stopping at 
%      junctions and image borders.  
%
%  Returns the fragments as a cell array of [Mx2] matrices which contain 
%  the M (x,y) coordinates of of each fragment's constituent elements, the 
%  fragment neighbor information (see below for details), and the new 
%  (filled-in / merged) segmentation.
%
%  'neighbors' is a struct with two fields:
%   - 'junction_fragmentlist': a list which contains for each junction an
%                              array of indices of the fragments that meet 
%                              there
%   - 'fragment_junctionlist': a list which contains for each fragment the
%                              indices of the junctions which bound it
%   - 'fragment_segments': a list which contains for each fragment the ID
%                          numbers of the segments found to the left and
%                          right of that fragment, stored [left right].
%
%  Thus, the indices of fragment i's neighbors on one end can be found by:
%    junction_fragmentlist{fragment_junctionlist{i}(1)}
%
%  If requested, also fits polynomials of degree 'order' to each fragment
%  and returns those as well (and their parameters).
%

% DEPENDENCIES:
%  - RGB2Lab.m
%  - cracks2fragments.m
%  - fill_in_segmentation.m
%  - fit_curves.m
%    - fit_poly_to_fragment.m
%  - remove_small_segments.m
%  - seg2cracks.m
% [- drawedges.m (will use it if available, not a problem if not)]
 
% Convert the image to Lab space so that color distances used below are
% more meaningful
if ~isempty(img)
    img_lab = RGB2Lab(img);
end

if(any(seg(:)==0))
    seg = fill_in_segmentation(img, seg);
end

% Find locations where there's a pixel jutting off diagonally from  a 
% corner of a segment.  We want to remove those and attach them to the 
% nearest-colored 4-connected neighbor
stragglers = seg~=seg(:,[2:end end]) & seg~=seg(:,[1 1:end-1]) & ...
    seg~=seg([2:end end],:) & seg~=seg([1 1:end-1],:);
if(any(stragglers(:)))
    seg(stragglers) = 0;
    while any(seg(:)==0)
        seg = fill_in_segmentation(img, seg, 0, 4);
        if any(seg(:)==0), disp('Trying again...'); end
    end
end

if(nargin < 3 || isempty(min_area))
    min_area = 0;
end

if(min_area > 0)
    seg = remove_small_segments(seg, img_lab, min_area, 1);
end

cracks = seg2cracks(seg);
[fragments, neighbors] = cracks2fragments(cracks, seg);

if(nargout>=4 || nargout==0)
    if(nargin < 4)
        order = 3;
    end
    
    [polyfragments, poly_params] = fit_curves(fragments, order);
end

if(nargout==0) 
    figure(gcf)
    subplot 121, hold off
    imagesc(img), axis image, hold on
    title('Segment Borders')
    
    subplot 122, hold off
    imagesc(img), axis image, hold on
    title(['Polynomial Fits to Fragments (order=' num2str(order)])
    
    if(exist('drawedges', 'file'))
        subplot 121, drawedges(fragments, 'rand');
        subplot 122, drawedges(polyfragments, 'rand');
    else
        subplot 121
        for(i=1:length(fragments))
            plot(fragments{i}(:,1), fragments{i}(:,2), 'r', 'LineWidth', 2);
        end
        subplot 122
        for(i=1:length(polyfragments))
            plot(polyfragments{i}(:,1), polyfragments{i}(:,2), 'r', 'LineWidth', 2);
        end
    end    
end
    