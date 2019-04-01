%{
% Filename: generatePbFile.m
% Project: GeometricContextDataset
% Created Date: Wednesday March 27th 2019
% Author: Feng Panhe
% -----
% Last Modified:
% Modified By:
% -----
% Copyright (c) 2019 Feng Panhe
%}
function generatePbFile(image_names, imdir, outdir)

    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    file_names = strtok(image_names, '.');

    exist_file_names = dir(fullfile(outdir, '*_pb.mat'));
    exist_file_names = {exist_file_names.name};
    exist_file_names = split(exist_file_names, '_pb.mat');
    not_exist_file_names = setdiff(file_names, exist_file_names);

    not_exist_file_num = numel(not_exist_file_names);


    if not_exist_file_num > 0

        parfor i = 1:not_exist_file_num
            
            fprintf('\n\n image %d / %d .\n', i, not_exist_file_num);
            f_name = not_exist_file_names{i};
            image_f = imread(fullfile(imdir, strcat(f_name, '.jpg')));
            pbim_file = fullfile(outdir, strcat(f_name, '_pb.mat'));
            pbim = pbCGTG_nonmax(double(image_f) / 255);
            mat_file_save(pbim_file, pbim);
            
        end

    end

end

function mat_file_save(mat_file, pbim)
    save(mat_file, 'pbim');
end
