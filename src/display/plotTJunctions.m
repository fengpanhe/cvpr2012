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
    % annotation_radius = im_size(1) / 10;
    j_pos = bndinfo.junctions.position;
    edge_endpoint_indexs = bndinfo.edges.junctions;
    mapshow(j_pos(:, 1), j_pos(:, 2), 'Color', 'b', 'DisplayType', 'point', 'Marker', 'o');

    for i = 1:size(j_pos, 1)
        e_num = numel(find(edge_endpoint_indexs == i));

        if e_num == 3
            % x = j_pos(i, 1) - annotation_radius / 2;
            % y = im_size(1) - j_pos(i, 2)- annotation_radius / 2;
            % if x < 0
            %     y = y - x;
            %     x = x - x;
            % end
            % if y < 0
            %     x = x - y;
            %     y = y - y;
            % end
            % x = x / im_size(2);
            % y = y / im_size(1);

            % dim = [x, y, annotation_radius / im_size(2), annotation_radius / im_size(1)];
            % disp(dim);
            % a = annotation('ellipse',dim);
            ann_x = [abs(j_pos(i, 1) - 50), j_pos(i, 1)] / im_size(2);
            ann_y = 1 - [abs(j_pos(i, 2) - 50), j_pos(i, 2)] / im_size(1);
            a = annotation('textarrow', ann_x, ann_y, 'HeadStyle', 'vback3', 'String', num2str(i));
            a.Color = 'red';
        end

    end

end
