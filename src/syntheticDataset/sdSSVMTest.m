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
function [ssvm_precision, mrf_precision, each_res] = sdSSVMTest(model_file)
    %sdSSVMTest - Description
    %
    % Syntax: [precision1, precision2] = sdSSVMTest(input)
    %
    % Long description

    dataset_path = fullfile('resources', 'SyntheticDataset');
    fileList = dir(fullfile(dataset_path, 'gtsave', '*_gt.mat'));
    fileList = {fileList.name};

    file_name_list = strtok(fileList, '_');

    [ssvm_precision, mrf_precision, each_res] = multImageTJTest(dataset_path, file_name_list, model_file);
end
