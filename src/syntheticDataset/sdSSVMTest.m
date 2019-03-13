%{
% Filename: sdSSVMTest.m
% Project: syntheticDataset
% Created Date: Monday February 25th 2019
% Author: Feng Panhe
% -----
% Last Modified:
% Modified By:
% -----
% Copyright (c) 2019 Feng Panhe
%}
function [ssvm_precision, mrf_precision] = sdSSVMTest(model_file, off_plot)
    %sdSSVMTest - Description
    %
    % Syntax: [precision1, precision2] = sdSSVMTest(input)
    %
    % Long description

    dataset_path = fullfile('resources', 'SyntheticDataset');
    fileList = dir(fullfile(dataset_path, 'gtsave', '*_gt.mat'));
    fileList = {fileList.name};

    file_name_list = strtok(fileList, '_');

    tj_num = 0;
    ssvm_correct_tj_num = 0;
    mrf_correct_tj_num = 0;
    ssvm_precision_list = zeros(1, numel(file_name_list));
    mrf_precision_list = zeros(1, numel(file_name_list));

    for i = 1:numel(file_name_list)
        disp(file_name_list{i});
        [ssvm_precision, mrf_precision, ...
            ssvm_predict_score, mrf_predict_score, ...
            ssvm_tj_errata, mrf_tj_errata] = ...
            TJTest(model_file, dataset_path, file_name_list{i}, off_plot);
        ssvm_precision_list(i) = ssvm_precision;
        mrf_precision_list(i) = mrf_precision;
        tj_num = tj_num + numel(ssvm_tj_errata) / 3;
        ssvm_correct_tj_num = ssvm_correct_tj_num + sum(ssvm_tj_errata) / 3;
        mrf_correct_tj_num = mrf_correct_tj_num + sum(mrf_tj_errata) / 3;
    end

    disp(file_name_list);
    disp(ssvm_precision_list);
    disp(mrf_precision_list);
    ssvm_precision = ssvm_correct_tj_num / tj_num;
    mrf_precision = mrf_correct_tj_num / tj_num;

end
