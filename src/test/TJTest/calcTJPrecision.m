function [precision, errata] = calcTJPrecision(tid, label, predictScore)
    %calcTJPrecision - Description
    %
    % Syntax: count = calcTJPrecision(input)
    %
    % tid, label, predictScore 均是每三行一个 TJ
    %
    res = [tid, label, predictScore];
    count = 0;
    groupNum = size(res, 1) / 3;
    res = sortrows(res, [1, 3]);
    errata = zeros(numel(tid), 1);

    for i = 1:groupNum
        i3 = i * 3;
        rows_index = [i3 - 2; i3 - 1; i3];

        if all([1; 2; 3] == res(rows_index, 2))
            count = count + 1;
            errata(rows_index, :) = 1;
        end

    end

    precision = double(count) / double(groupNum);
end
