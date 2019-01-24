function bndinfo = gcdGtPredict(image_name, model_file_path)
    %gcdSSVMPredict - Description
    %
    % Syntax: bndinfo = gcdGtPredict(image_name)
    %
    % Long description
    gcd_gt_bndinfo_path = fullfile('resources', 'GeometricContextDataset', 'gtsave');
    gt_bndinfo_file = fullfile(gcd_gt_bndinfo_path, [image_name '_gt.mat']);

    if ~exist(gt_bndinfo_file, 'file')
        error('File %s does not exist.', gt_bndinfo_file);
    end

    [X, Y, infos] = gcd2ssvmXY({sprintf('%s%s', image_name, '_gt.mat')});
    load(gt_bndinfo_file, 'bndinfo');

    predict_score = SSVMPredict(X, Y, model_file_path);

    mrf_matrix = [infos, X(:, 1), predict_score];
    mrf_matrix = MRFEnergy_mex(21, double(mrf_matrix), size(mrf_matrix));

    predict_score = mrf_matrix(:, 4);
    res = [X(:, 1), Y(:, 1), predict_score];
    bndinfo.res = res;
end
