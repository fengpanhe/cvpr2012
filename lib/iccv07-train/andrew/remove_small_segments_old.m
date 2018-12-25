function [new_seg] = remove_small_segments(seg, img, min_area, reorder);
%
% [new_seg] = remove_small_segments(seg, img, min_area);
%
%  Remove segments below a minimum size, assigning them to a neighbor
%  segment with the most similar average color.
%

[nrows,ncols] = size(seg);
onborder = false(nrows,ncols);
onborder([1 end],:) = true;
onborder(:,[1 end]) = true;

stats = regionprops(seg, 'PixelIdxList');
all_indices = {stats(:).PixelIdxList};

new_seg = seg;

num_segments = max(seg(:));
areas = hist(seg(:), num_segments);
[sorted_areas, sort_index] = sort(areas);
too_small = find(sorted_areas<min_area);
neighbor_offsets = [-1 1 -nrows nrows -1-nrows -1+nrows 1-nrows 1+nrows];
segment_colors = cell(num_segments,1);
for(i=1:length(too_small))
    if(sorted_areas(i) > 0)
        segment_indices = all_indices{sort_index(i)};  %find(new_seg==sort_index(i));

        % make sure this segment hasn't already had a smaller segment merged
        % with it, making it big enough to pass the threshold
        if(length(segment_indices)<min_area)
            % Find this segment's neighbors, being careful at boundaries
            neighbor_indices = segment_indices*ones(1,8) + ones(length(segment_indices),1)*neighbor_offsets;
            [y,x] = ind2sub([nrows ncols], neighbor_indices);
            neighbor_indices(x<1 | x>ncols | y<1 | y>nrows) = [];
            neighbors = unique(new_seg(neighbor_indices));
            neighbors(neighbors==sort_index(i)) = [];

            % Compute this segment's average color if not already stored
            if(isempty(segment_colors{sort_index(i)}))
                color_indices = [segment_indices segment_indices+nrows*ncols segment_indices+2*nrows*ncols];
                segment_colors{sort_index(i)} = mean(img(color_indices), 1);
            end

            % Compute the neighbor segments average colors if not already
            % stored
            num_neighbors = length(neighbors);
            neighbor_indices = cell(1,num_neighbors);
            for(j=1:num_neighbors)
                if(isempty(segment_colors{neighbors(j)}))
                    neighbor_indices(j) = all_indices(neighbors(j)); % find(new_seg==neighbors(j));
                    color_indices = [neighbor_indices{j}(:) neighbor_indices{j}(:)+nrows*ncols neighbor_indices{j}(:)+2*nrows*ncols];
                    segment_colors{neighbors(j)} = mean(img(color_indices));
                end
            end

            % Lookup the neighbors' colors
            neighbor_colors = vertcat(segment_colors{neighbors});

            % Compute the distance between this segment's average color and the
            % neighbors' colors, and find the one with the most similar color.
            this_seg_color = ones(num_neighbors,1)*segment_colors{sort_index(i)};
            color_dist = sum( (this_seg_color-neighbor_colors).^2 , 2);
            [min_dist, which_neighbor] = min(color_dist);

            % Make this segment part of that neighbor:
            new_seg(segment_indices) = neighbors(which_neighbor);
            
            % Update the neighbor segment's avg. color now that it has grown:
            segment_area = length(segment_indices);
            neighbor_area = length(neighbor_indices{which_neighbor});
            total_area =  segment_area + neighbor_area;
            segment_colors{neighbors(which_neighbor)} = ...
                (segment_area*segment_colors{sort_index(i)} + ...
                 neighbor_area*segment_colors{neighbors(which_neighbor)}) / total_area;
        end
    end
end

if exist('reorder', 'var') && reorder
    nseg = max(new_seg(:));        
    diffs = setdiff((1:nseg), unique(new_seg));
    shifts = zeros(1, nseg);  shifts(diffs) = 1;  shifts = cumsum(shifts);
    new_seg = new_seg - shifts(new_seg);
end

if(nargout==0)
    figure
    subplot 131, imagesc(seg), axis image
    title([num2str(length(unique(seg(:)))) ' Segments Originally'])
    subplot 132, imagesc(new_seg), axis image
    title([num2str(length(unique(new_seg(:)))) ' Segments After Merging'])
    subplot 133, imagesc(seg~=new_seg), axis image
    title('Changes')
end
    
return;