function modelFilePath = gcdSSVMTrain(matFileNameList)
%train - gcdSSVMTrain
%
% Syntax: modelFilePath = gcdSSVMTrain(X, Y)
%
% Long description
    gcdGtPath = 'resources/GeometricContextDataset/gtsave/';
    gcdImPath = 'resources/GeometricContextDataset/images/';

    imNameList = strtok(matFileNameList, '_');
    matFileNum = numel(matFileNameList);
    X = [];
    Y = [];
    for index = 1:matFileNum
        matFile = [gcdGtPath, matFileNameList{index}];
        imFile = [gcdImPath, imNameList{index}, '.jpg'];
        if exist(matFile, 'file') && exist(imFile, 'file')
            im = imread(imFile);
            load(matFile);
            [bndinfo2, err] = updateBoundaryInfo2(bndinfo, bndinfo.labels);  
            combinedFeatures = getCombinedFeatures(bndinfo2, im);
            [features, lables] = getSSVMClassifierFeatures(bndinfo2, combinedFeatures, 'train');
            X = [X; features];
            Y = [Y; labels];
        end
    end
    modelFilePath = SSVMTrain(X, Y);
end


function modelFilePath = test_gcdSSVMTrain()
    matFileNameList = {'alley01_gt.mat', 'alley02_gt.mat'};
    modelFilePath = gcdSSVMTrain(matFileNameList);
end