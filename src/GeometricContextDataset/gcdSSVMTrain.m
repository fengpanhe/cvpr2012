function model_file = gcdSSVMTrain(train_file_list, gcd_ssvm_model_file_path)
    %train - gcdSSVMTrain
    %
    % Syntax: modelFilePath = gcdSSVMTrain(X, Y)
    %
    % Long description
    if ~exist('gcd_ssvm_model_file_path', 'var') || isempty(gcd_ssvm_model_file_path)
        gcd_ssvm_model_file_path = fullfile('resources', 'SSVMmodel', 'gcd_ssvm.model');
    end

    X = [];
    Y = [];

    for i = 1:numel(train_file_list)
        [X_, Y_, ~] = gcd2ssvmXY(train_file_list{i});

        if size(X, 1) > 0
            X_(:, 1) = X_(:, 1) + X(end, 1);
        end

        X = [X; X_];
        Y = [Y; Y_];
    end

    model_file = SSVMTrain(X, Y, gcd_ssvm_model_file_path);
end
