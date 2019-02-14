%{
% Filename: generateSyntheticDataset.m
% Project: syntheticDataset
% Created Date: Tuesday February 12th 2019
% Author: Feng Panhe
% -----
% Last Modified:
% Modified By:
% -----
% Copyright (c) 2019 Feng Panhe
%}
function res = generateSyntheticDataset(image_num, dest_dir_path)
    %generateSyntheticDataset - Description
    %
    % Syntax: res = generateSyntheticDataset(input)
    %
    % Long descriptions

    if ~exist('dest_dir_path', 'var') || isempty(dest_dir_path)
        dest_dir_path = fullfile('resources', 'SyntheticDataset');
    end

    if ~exist('image_num', 'var') || isempty(image_num)
        image_num = 1000;
    end

    gt_dir_path = fullfile(dest_dir_path, 'gtsave');
    im_dir_path = fullfile(dest_dir_path, 'images');

    max_shape_num = 10;

    image_size = [1200, 800];

    set(0,'DefaultFigureVisible', 'off');
    
    for i = 1:image_num
        im_name = sprintf('sd%d', i);
        im_file = fullfile(im_dir_path, strcat(im_name, '.jpg'));
        gt_file = fullfile(im_dir_path, strcat(im_name, '_gt.mat'));
        
        figure('Color', 'w', 'Position', [0, 0, image_size]);
        set(gca, 'position', [0 0 1 1]);
        shape_num = randi(max_shape_num);
        
        for j = 1:shape_num
            shape_type_id = randi(2);

            switch shape_type_id
                case 1
                    drawOneRandomTriangle([0, image_size(1)], [0, image_size(2)], rand(1, 3));
                case 2
                    drawOneRandomRectangle([0, image_size(1)], [0, image_size(2)], rand(1, 3));
            end

        end

        axis off;
        print(['-f' num2str(1)], '-djpeg', [im_file(1:end - 4) '.jpg']);
        close all;
    end

    set(0,'DefaultFigureVisible', 'on');

    % x1 = [100, 80, 1200];
    % y1 = [400, 300, 100];
    % x2 = [20, 20, 1200];
    % y2 = [40, 80, 600];
    % x1box = [100, 80, 1200, 100];
    % y1box = [400, 300, 100, 400];
    % x2box = [20, 20, 1200, 20];
    % y2box = [40, 80, 600, 40];
    % figure('Color', 'w', 'Position', [0, 0, image_size]);
    % set(gca, 'position', [0 0 1 1])
    % patch(x2, y2, 'green')
    % patch(x1, y1, 'green')
    % [xi, yi] = polyxpoly(x1box, y1box, x2box, y2box);
    % mapshow(xi, yi, 'DisplayType', 'point', 'Marker', 'o');
    % axis off;
    % %box off;
    % print(['-f' num2str(1)], '-djpeg', [im_file(1:end - 4) '.jpg']);

end

function drawOneRandomTriangle(x_range, y_range, shape_color)
    %drawOneRandomTriangle - Description
    %
    % Syntax:  drawOneRandomTriangle(x_range, y_range, shape_color)
    %
    % Long description
    shape_edge_num = 3;
    x = randi(x_range(2) - x_range(1), 1, shape_edge_num) + x_range(1);
    y = randi(y_range(2) - y_range(1), 1, shape_edge_num) + y_range(1);

    patch(x, y, shape_color);
end

function drawOneRandomRectangle(x_range, y_range, shape_color)
    %drawOneRandomRectangle - Description
    %
    % Syntax:  drawOneRandomRectangle(x_range, y_range, shape_color)
    %
    % Long description
    shape_edge_num = 4;
    x = zeros(1, shape_edge_num);
    y = zeros(1, shape_edge_num);

    x_edge_length = randi(x_range(2) - x_range(1));
    y_edge_length = randi(y_range(2) - y_range(1));

    x(1) = randi(x_range(2) - x_range(1) - x_edge_length) + x_range(1);
    y(1) = randi(y_range(2) - y_range(1) - y_edge_length) + y_range(1);

    x(2) = x(1);
    y(2) = y(1) + y_edge_length;

    x(3) = x(1) + x_edge_length;
    y(3) = y(1) + y_edge_length;

    x(4) = x(1) + x_edge_length;
    y(4) = y(1);

    patch(x, y, shape_color);
end
