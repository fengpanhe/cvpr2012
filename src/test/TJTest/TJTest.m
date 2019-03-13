function [ssvm_precision, mrf_precision, ssvm_predict_score, mrf_predict_score, ssvm_tj_errata, mrf_tj_errata] = TJTest(model_file, dataset_path, image_name, off_plot)
    %TJTest - Description
    %
    % Syntax: output = TJTest(input)
    %
    % example:
    % [ssvm_precision, mrf_precision] = TJTest('resources/SSVMmodel/ssvm.model', 'resources/SyntheticDataset', 'sd10', 1);

    if ~exist('off_plot', 'var') || isempty(off_plot)
        off_plot = false;
    end

    if ~exist(model_file, 'file')
        error('File \"%s\" does not exist.', model_file);
    end

    [gt_file, image_file, pbim_file] = getDatasetInfoFile(dataset_path, image_name);

    image_f = imread(image_file);
    load(gt_file, 'bndinfo');
    load(pbim_file, 'pbim');

    bndinfo.pbim = pbim;

    combined_features = getCombinedFeatures(bndinfo, image_f);

    bndinfo.combined_features = combined_features;

    item_infos = getSSVMClassifierFeatures(bndinfo, combined_features, 'train');

    Y = item_infos(:, 1);
    edge_id = item_infos(:, 2);
    X = item_infos(:, 3:end);

    % rank-svm 预测
    ssvm_predict_score = SSVMPredict(X, Y, model_file);

    image_id = zeros(numel(edge_id), 1);
    [ssvm_precision, ssvm_tj_errata] = calcTJPrecision(X(:, 1), Y(:, 1), ssvm_predict_score);

    tj_num = bndinfo.combined_features.TJnum;
    tj_infos = bndinfo.combined_features.TJInfo;
    edges_segid = [];

    for tj_i = 1:tj_num
        tj_edge_segid = cell2mat(tj_infos{tj_i}.edge_segid);
        edges_segid = [edges_segid; tj_edge_segid];
    end

    % mrf 优化
    penalty = 20;
    mrfMatrix = [image_id, edges_segid, X(:, 1), ssvm_predict_score];
    mrfMatrix = MRFEnergy_mex(penalty, double(mrfMatrix), size(mrfMatrix));
    mrf_predict_score = mrfMatrix(:, 4);

    % 有向图
    seg_num = max(max(bndinfo.wseg));

    tj_num = bndinfo.combined_features.TJnum;
    adjacency_matrix = zeros(seg_num);

    for tj_i = 1:tj_num
        
        index_tmp = tj_i * 3;
        tj_edge_score = mrf_predict_score(index_tmp - 2:index_tmp, 1);
        tj_edge_segid = edges_segid(index_tmp - 2:index_tmp, 1);
        segid_score = [tj_edge_segid, tj_edge_score];
        segid_score = sortrows(segid_score, 2);
        adjacency_matrix(segid_score(3, 1), segid_score(2, 1)) = adjacency_matrix(segid_score(3, 1), segid_score(2, 1)) + segid_score(3, 2) - segid_score(2, 2);
        adjacency_matrix(segid_score(3, 1), segid_score(1, 1)) = adjacency_matrix(segid_score(3, 1), segid_score(1, 1)) + segid_score(3, 2) - segid_score(1, 2);
        adjacency_matrix(segid_score(2, 1), segid_score(1, 1)) = adjacency_matrix(segid_score(2, 1), segid_score(1, 1)) + segid_score(2, 2) - segid_score(1, 2);
    end

    G = digraph(adjacency_matrix);
    plot(G, 'Layout', 'force', 'EdgeLabel', G.Edges.Weight);

    % 结果处理

    [mrf_precision, mrf_tj_errata] = calcTJPrecision(X(:, 1), Y(:, 1), mrf_predict_score);

    save(gt_file, 'bndinfo');

    if ~off_plot
        fprintf('正在画图片 %s 的 TJ...\n', image_name);
        plotTJTestResult(bndinfo, ssvm_predict_score);
    end

end

function plotTJTestResult(bndinfo, predict_score)
    %poltTJTestResult - Description
    %
    % Syntax: output = poltTJTestResult(bndinfo, tj_errata)
    %
    % Long description

    [~, image_name, ~] = fileparts(bndinfo.imname);

    image_save_dir = fullfile('result', 'tj_test');

    if ~exist(image_save_dir, 'dir')
        mkdir(image_save_dir);
    end

    edges_indices = bndinfo.edges.indices;
    edge_num = numel(edges_indices);
    blabels = bndinfo.edges.boundaryType;
    image_size = bndinfo.imsize;
    tj_infos = bndinfo.combined_features.TJInfo;
    tj_num = bndinfo.combined_features.TJnum;

    edges_lab = blabels(1:end / 2) + 2 * blabels(end / 2 + 1:end);
    edges_xy = cell(edge_num, 1);

    for edge_i = 1:edge_num
        [ey, ex] = ind2sub(image_size, edges_indices{edge_i});
        edges_xy{edge_i} = [ex, ey];
    end

    set(0, 'DefaultFigureVisible', 'off');

    for tj_i = 1:tj_num
        tj_edge_id = cell2mat(tj_infos{tj_i}.edgeId);
        tj_edge_flip = logical(cell2mat(tj_infos{tj_i}.edgeFlip));

        figure('Color', 'w', 'Position', [0, 0, image_size([2, 1])]);
        set(gca, 'position', [0 0 1 1]);
        set(gca, 'YDir', 'reverse');
        hold on;
        plot([0, 0, image_size(2), image_size(2)], [0, image_size(1), 0, image_size(1)], '.');

        predict_lab = zeros(3, 1);
        ps_i = tj_i * 3 - 2;
        predict_lab(1) = 1 + xor(predict_score(ps_i) < predict_score(ps_i + 1), tj_edge_flip(1));
        predict_lab(2) = 1 + xor(predict_score(ps_i + 1) < predict_score(ps_i + 2), tj_edge_flip(2));
        predict_lab(3) = 1 + xor(predict_score(ps_i + 2) < predict_score(ps_i), tj_edge_flip(3));
        plotOneTJ(edges_xy(tj_edge_id, :), predict_lab, 'r', image_size);
        plotOneTJ(edges_xy(tj_edge_id, :), edges_lab(tj_edge_id), 'b', image_size);
        save_file_name = strcat(image_name, '_tj', num2str(tj_i), '.jpg');

        if all(predict_lab == edges_lab(tj_edge_id))
            save_file_name = strcat('correct_', save_file_name);
        else
            save_file_name = strcat('incorrect_', save_file_name);
        end

        print('-djpeg', '-r0', fullfile(image_save_dir, save_file_name));
        hold off;
        close all;
    end

    set(0, 'DefaultFigureVisible', 'on');
end

function plotOneTJ(edges_xy, edges_lab, ann_color, image_size)
    %poltOneTJ - Description
    %
    % Syntax: output = poltOneTJ(edges, labs, edge_color)
    %
    % Long description
    arrowdist = ceil(sqrt(image_size(1).^2 + image_size(2).^2) / 10);

    for k = 1:size(edges_xy, 1)

        if edges_lab(k) > 0

            edge_xy = edges_xy{k};
            ex = edge_xy(:, 1);
            ey = edge_xy(:, 2);
            npix = numel(ex);

            narrows = ceil(npix / arrowdist);

            epos = ceil((1:narrows) / (narrows + 1) * npix);
            plot(ex, ey, 'Color', 'k', 'LineWidth', 1);

            for j = 1:numel(epos)

                a_xy = edge_xy(epos(j), :);

                if edges_lab(k) == 1
                    xy1 = edge_xy(max(epos(j) - 10, 1), :);
                    xy2 = edge_xy(min(epos(j), npix), :);
                else % blabels(k)==2;
                    xy1 = edge_xy(min(epos(j) + 10, npix), :);
                    xy2 = edge_xy(min(epos(j), npix), :);
                end

                theta = atan2(xy2(2) - xy1(2), xy2(1) - xy1(1));

                as_xy = a_xy - [cos(theta), sin(theta)];

                a_xy(1) = (a_xy(1) - 1) / (image_size(2) - 1);
                a_xy(2) = 1 - (a_xy(2) - 1) / (image_size(1) - 1);
                as_xy(1) = (as_xy(1) - 1) / (image_size(2) - 1);
                as_xy(2) = 1 - (as_xy(2) - 1) / (image_size(1) - 1);
                a_xy(1) = max(a_xy(1), 0);
                as_xy(1) = max(as_xy(1), 0);
                annotation('arrow', [as_xy(1) a_xy(1)], [as_xy(2) a_xy(2)], 'Color', ann_color, 'LineStyle', 'none', ...
                    'HeadStyle', 'vback3', 'HeadWidth', 17, 'HeadLength', 10);
            end

        end

    end

end
