function prepareGdataDir(image_names, imsegs_list, outdir)
    %myFun - Description
    %
    % Syntax: prepareGdataDir(input)
    %
    % Long description
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    file_names = strtok(image_names, '.');

    exist_file_names = dir(fullfile(outdir, '*_gdata.mat'));
    exist_file_names = {exist_file_names.name};
    exist_file_names = split(exist_file_names, '_gdata.mat');
    not_exist_file_indexs = ~ismember(file_names, exist_file_names);

    if any(not_exist_file_indexs) > 0
        disp('prepareGdataDir');
        for imsegs = imsegs_list(not_exist_file_indexs)
            f_name = strtok(imsegs.imname, '.');
            save(fullfile(outdir, [f_name '_gdata.mat']), 'imsegs');
        end
    end

end
