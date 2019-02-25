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
function [precision1, precision2, res1, res2] = sdSSVMTest(modelFilePath, is_new)
    %sdSSVMTest - Description
    %
    % Syntax: [precision1, precision2] = sdSSVMTest(input)
    %
    % Long description
    if ~exist('is_new', 'var') && isempty(is_new)
        is_new = false;
    end

    sdGtPath = fullfile('resources', 'SyntheticDataset', 'gtsave');
    fileList = dir(fullfile(sdGtPath, '*_gt.mat'));
    fileList = {fileList.name};

    file_name_list = strtok(fileList, '_');

    ssvmXYFile = 'result/tmp/sd_ssvmXY.mat';

    if exist(ssvmXYFile, 'file') &&~is_new
        load(ssvmXYFile, 'X', 'Y', 'infos');
    else
        X = [];
        Y = [];
        infos = [];

        for i = 1:1
            [X_, Y_, infos_] = sd2ssvmXY(file_name_list{i}, is_new);

            if size(X, 1) > 0
                X_(:, 1) = X_(:, 1) + X(end, 1);
            end

            X = [X; X_];
            Y = [Y; Y_];
            infos_(:, end + 1) = i;
            infos = [infos; infos_];
        end

        infos(:, [2, 1]) = infos(:, [1, 2]);
        save(ssvmXYFile, 'X', 'Y', 'infos');
    end

    disp(file_name_list{1:1});
    testDateIndex = 1:size(X, 1);
    predictScore = SSVMPredict(X(testDateIndex, :), Y(testDateIndex, :), modelFilePath);

    res1 = zeros(size(X, 1), 4);
    res1(:, 1:3) = [X(testDateIndex, 1), infos(:, 2), predictScore];
    res1(:, 4) = 1;
    indexs = (1:size(X, 1) / 3) * 3 - 2;
    tmp_ = res1(indexs, 3) < res1(indexs + 1, 3);
    res1(indexs(tmp_), 4) = 2;
    tmp_ = res1(indexs + 1, 3) < res1(indexs + 2, 3);
    res1(indexs(tmp_) + 1, 4) = 2;
    tmp_ = res1(indexs + 2, 3) < res1(indexs, 3);
    res1(indexs(tmp_) + 2, 4) = 2;

    precision1 = calcPrecision(X(testDateIndex, 1), Y(testDateIndex, 1), predictScore);

    mrfMatrix = [infos(testDateIndex, :), X(testDateIndex, 1), predictScore];
    mrfMatrix = MRFEnergy_mex(21, double(mrfMatrix), size(mrfMatrix));
    predictScore = mrfMatrix(:, 4);

    res2 = zeros(size(X, 1), 4);
    res2(:, 1:3) = [X(testDateIndex, 1), infos(:, 2), predictScore];
    res2(:, 4) = 1;
    indexs = (1:size(X, 1) / 3) * 3 - 2;
    tmp_ = res2(indexs, 3) < res2(indexs + 1, 3);
    res2(indexs(tmp_), 4) = 2;
    tmp_ = res2(indexs + 1, 3) < res2(indexs + 2, 3);
    res2(indexs(tmp_) + 1, 4) = 2;
    tmp_ = res2(indexs + 2, 3) < res2(indexs, 3);
    res2(indexs(tmp_) + 2, 4) = 2;

    precision2 = calcPrecision(X(testDateIndex, 1), Y(testDateIndex, 1), predictScore);

end

function precision = calcPrecision(tid, label, predictLabel)
    %myFun - Description
    %
    % Syntax: count = myFun(input)
    %
    % Long description
    res = [tid, label, predictLabel];
    count = 0;
    groupNum = size(res, 1) / 3;
    res = sortrows(res, [1, 3]);

    for i = 1:groupNum
        i3 = i * 3;
        rows = [i3 - 2; i3 - 1; i3];
        a = eq([1; 2; 3], res(rows, 2));

        if sum(a) == 3
            count = count + 1;
        end

    end

    precision = double(count) / double(groupNum);
end
