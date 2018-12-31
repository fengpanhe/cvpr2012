function predictFilePath = predictSSVM(testFilePath, modelFilePath)
%predictSSVM - Description
%
% Syntax: predictionFilePath = predictSSVM(input)
%
% Long description
    ssvmpredict = 'lib/svm_rank/build/svm_rank_classify';
    predictFilePath = 'resources/SSVMmodel/predict.txt';
    cmd = [ssvmpredict ' ' testFilePath ' ' modelFilePath ' ' predictFilePath];
    disp(cmd);
    system(cmd);
end