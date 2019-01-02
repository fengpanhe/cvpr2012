function modelFilePath = gcdSSVMTrain(matFileNameList)
    %train - gcdSSVMTrain
    %
    % Syntax: modelFilePath = gcdSSVMTrain(X, Y)
    %
    % Long description

    [X, Y] = gcd2ssvmXY(matFileNameList);

    modelFilePath = SSVMTrain(X, Y);
end

function modelFilePath = test_gcdSSVMTrain()
    matFileNameList = {'alley01_gt.mat', 'alley02_gt.mat'};
    modelFilePath = gcdSSVMTrain(matFileNameList);
end
