%{
% Filename: writeContourDepthResults.m
% Project: display
% Created Date: Thursday February 21st 2019
% Author: Feng Panhe
% -----
% Last Modified:
% Modified By:
% -----
% Copyright (c) 2019 Feng Panhe
%}
function writeContourDepthResults(bndinfo_file, image_file, dest_dir_path, is_ex)
    %myFun - Draw a depth image on the original image
    %
    % Syntax:  writeContourDepthResults(fileNameList)
    %
    % 是 11年iccv的 iccvWriteContourDepthResults 的修改。
    % 加入了面的深度显示，深色在前，数字小。
    %
    % Inputs:
    %   bndinfo_file: The path of bndinfo file
    %   image_file: The path of image file
    %   dest_dir_path: Result directory path
    %   is_ex:

    if ~exist(bndinfo_file, 'file')
        error('File \"%s\" does not exist.', bndinfo_file);
    end

    if ~exist(image_file, 'file')
        error('File \"%s\" does not exist.', image_file);
    end

    if ~exist('is_ex', 'var') && isempty(is_ex)
        is_ex = false;
    end

    [~, image_name, image_ext_name] = fileparts(image_file);
    [~, bndinfo_name, ~] = fileparts(bndinfo_file);
    im = im2double(imread(image_file));
    im = rgb2gray(im);
    if is_ex

        load(bndinfo_file, 'bndinfo');
        lab = bndinfo.edges.boundaryType;
        lab = lab(1:end / 2) + 2 * lab(end / 2 + 1:end);
        out_name = fullfile(dest_dir_path, strcat(bndinfo_name, image_ext_name));
        printOcclusionResult(im, bndinfo, lab, out_name, 1);
    else
        load(bndinfo_file, 'bndinfo');
        out_name = fullfile(dest_dir_path, strcat(image_name, '_gcd', image_ext_name));

        [bndinfo2, ~] = updateBoundaryInfo2(bndinfo, bndinfo.labels);
        bndinfo2.names = bndinfo.names;
        bndinfo2.type = bndinfo.type;
        bndinfo2 = orderfields(bndinfo2);

        bndinfo2 = processGtBoundaryLabels(bndinfo2);

        lab = bndinfo2.edges.boundaryType;
        lab = (lab(1:end / 2) > 0) + 2 * (lab(end / 2 + 1:end) > 0);
        printOcclusionResult(im, bndinfo2, lab, out_name, 1);
    end

end
