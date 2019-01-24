%{
% Filename: plotTJunctions.m
% Project: display
% Created Date: Tuesday January 22nd 2019
% Author: Feng Panhe
% -----
% Last Modified:
% Modified By:
% -----
% Copyright (c) 2019 Feng Panhe
%}
function plotTJunctions(bndinfo)
%plotTJunctions - Description
%
% Syntax: plotTJunctions(bndinfo)
%
% Long description
    hold on
    
    im_size = bndinfo.imsize;
    annotation_radius = im_size(1) / 10;
    junction_positions = bndinfo.junctions.position;
    edge_endpoint_indexs = bndinfo.edges.junctions;
    
    for i = 1:size(junction_positions, 1)
        e_num = numel(find(edge_endpoint_indexs == i));
        if e_num == 3
            x = junction_positions(i, 1) - annotation_radius / 2;
            y = im_size(1) - junction_positions(i, 2)- annotation_radius / 2;
            if x < 0
                y = y - x;
                x = x - x;
            end
            if y < 0
                x = x - y;
                y = y - y;
            end
            x = x / im_size(2);
            y = y / im_size(1);
            
            dim = [x, y, annotation_radius / im_size(2), annotation_radius / im_size(1)];
            disp(dim);
            a = annotation('ellipse',dim);
            a.Color = 'red';
        end
    end
end