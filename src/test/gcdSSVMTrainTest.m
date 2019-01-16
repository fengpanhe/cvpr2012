function res = gcdSSVMTrainTest()
    %gcdSSVMTrainTest - Description
    %
    % Syntax: gcdSSVMTrainTest()
    %
    % Long description
    gcdGtPath = fullfile('resources', 'GeometricContextDataset', 'gtsave');
    fileList = dir(fullfile(gcdGtPath, '*_gt.mat'));
    fileList = {fileList.name};

    ssvmXYFile = 'result/tmp/ssvmXY.mat';
    if exist(ssvmXYFile, 'file')
        load(ssvmXYFile);
    else
        [X, Y, infos] = gcd2ssvmXY(fileList);
        save(ssvmXYFile, 'X', 'Y', 'infos');
    end

    dateGroupNum = size(X, 1) / 3;
    mid = ceil(dateGroupNum * 2 / 3) * 3;
    % randomOrder = randperm(dateGroupNum);
    trainDateIndex = 1:dateGroupNum * 3;
    testDateIndex = mid + 1:dateGroupNum * 3;

    % trainDateIndex = trainDateIndex * 3;
    % trainDateIndex = [trainDateIndex, trainDateIndex - 1, trainDateIndex - 2]
    % trainDateIndex = sort(trainDateIndex);

    % testDateIndex = testDateIndex * 3;
    % testDateIndex = [testDateIndex, testDateIndex - 1, testDateIndex - 2]
    % testDateIndex = sort(testDateIndex);

    modelFilePath = SSVMTrain(X(trainDateIndex, :), Y(trainDateIndex, :));
    predictScore = SSVMPredict(X(testDateIndex, :), Y(testDateIndex, :), modelFilePath);

    res = [X(testDateIndex, 1), Y(testDateIndex, 1), predictScore];
    res2 = sortrows(res, [1,3]);

    count = 0;
    groupNum = size(res2, 1) / 3;
    for i = 1:groupNum
        i3 = i * 3;
        rows = [i3-2; i3-1; i3];
        a = eq([1;2;3], res2(rows, 2));
        if sum(a) == 3
            count = count + 1;
        end
    end
    disp(count);
end
