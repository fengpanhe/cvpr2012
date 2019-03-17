function [ssvm_precision, mrf_precision, each_image_res] = multImageTJTest(dataset_path, file_name_list, model_file)
    file_num = numel(file_name_list);
    ssvm_precision_list = zeros(file_num, 1);
    mrf_precision_list = zeros(file_num, 1);
    ssvm_tj_errata_cell = cell(file_num, 1);
    mrf_tj_errata_cell = cell(file_num, 1);

    parfor i = 1:file_num
        [ssvm_precision_list(i), mrf_precision_list(i), ...
            ~, ~, ...
            ssvm_tj_errata, mrf_tj_errata] = ...
            oneImageTJTest(model_file, dataset_path, file_name_list{i}, 1);

        ssvm_tj_errata_cell{i} = ssvm_tj_errata;
        mrf_tj_errata_cell{i} = mrf_tj_errata;
    end

    image_tj_nums = cellfun(@numel, ssvm_tj_errata_cell) / 3;
    tj_num = sum(image_tj_nums);
    ssvm_correct_tj_num = sum(cellfun(@sum, ssvm_tj_errata_cell)) / 3;
    mrf_correct_tj_num = sum(cellfun(@sum, mrf_tj_errata_cell)) / 3;

    % 结果处理
    each_image_res = [file_name_list', num2cell(ssvm_precision_list), num2cell(mrf_precision_list), num2cell(image_tj_nums)];
    disp([{'image name', 'rank-svm precision', 'mrf precision', 'tj_num'}; each_image_res]);
    ssvm_precision = ssvm_correct_tj_num / tj_num;
    mrf_precision = mrf_correct_tj_num / tj_num;
end