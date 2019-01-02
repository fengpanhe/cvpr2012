function [X, Y] = gcd2ssvmXY(matFileNameList)
    %gcd2ssvmXY - 从gcd数据集生成ssvm的 X Y 形式的数组
    %
    % Syntax: [X, Y] = gcd2ssvmXY(matFileNameList)
    %
    % Long description
    gcdGtPath = 'resources/GeometricContextDataset/gtsave/';
    gcdImPath = 'resources/GeometricContextDataset/images/';
    gcdPbimPath = 'result/tmp/pbim/';

    imNameList = strtok(matFileNameList, '_');
    matFileNum = numel(matFileNameList);
    X = [];
    Y = [];
    for index = 1:matFileNum

        matFile = [gcdGtPath, matFileNameList{index}];
        imFile = [gcdImPath, imNameList{index}, '.jpg'];
        pbimFile = [gcdPbimPath, imNameList{index}, '_pbim.mat'];

        if ~exist(matFile, 'file')
            error('File \"%s\" does not exist.', matFile);
        end

        if ~exist(imFile, 'file')
            error('File \"%s\" does not exist.', imFile);
        end

        im = imread(imFile);
        load(matFile);

        if exist(pbimFile, 'file')
            load(pbimFile);
        else
            disp('pb');
            pbim = pbCGTG_nonmax(double(im) / 255);
            save(pbimFile, 'pbim');
        end

        [bndinfo2, err] = updateBoundaryInfo2(bndinfo, bndinfo.labels);
        bndinfo2.pbim = pbim;
        combinedFeatures = getCombinedFeatures(bndinfo2, im);
        [features, lables] = getSSVMClassifierFeatures(bndinfo2, combinedFeatures, 'train');

        if size(X, 1) > 0
            features(:, 1) = features(:, 1) + X(end, 1);
        end

        X = [X; features];
        Y = [Y; lables];
    end
end
