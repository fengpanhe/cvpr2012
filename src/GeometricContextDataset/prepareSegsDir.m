function prepareSegsDir(image_names, images_dir, pb_dir, save_dir)

    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    file_names = strtok(image_names, '.');

    exist_file_names = dir(fullfile(save_dir, '*_seg.mat'));
    exist_file_names = {exist_file_names.name};
    exist_file_names = split(exist_file_names, '_seg.mat');
    not_exist_file_names = setdiff(file_names, exist_file_names);

    not_exist_file_num = numel(not_exist_file_names);

    if not_exist_file_num > 0

        parfor i = 1:not_exist_file_num

            fprintf('\n\n prepareSegsDir: processing image %d / %d .\n', i, not_exist_file_num);
            f_name = not_exist_file_names{i};
            save_file = fullfile(save_dir, strcat(f_name, '_seg.mat'));

            image_f = im2double(imread(fullfile(images_dir, strcat(f_name, '.jpg'))));
            pbim = loadFile(fullfile(pb_dir, strcat(f_name, '_pb.mat')));

            wseg = pb2wseg(pbim, 100000);

            [edges, ~, neighbors, wseg, ~, ~, polyparams] = ...
                seg2fragments(double(wseg), image_f, 25);
            bndinfo = processBoundaryInfo3(wseg, edges, neighbors);
            bndinfo.imname = strcat(f_name, '.jpg');

            for k = 1:numel(polyparams)
                bndinfo.edges.polyparams{k} = single(polyparams{k});
            end

            imind = i;
            saveFile(save_file, bndinfo, imind);

        end

    end

end

function pbim = loadFile(file_name)
    load(file_name, 'pbim');
end

function saveFile(file_name, bndinfo, imind)
    save(file_name, 'bndinfo', 'imind');
end
