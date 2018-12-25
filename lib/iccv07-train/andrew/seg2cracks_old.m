function [crack_img] = seg2cracks(seg)
%
% [crack_img] = seg2cracks(seg)
% 
%  Converts a segmentation image into a crack-coded image.
%

% Bits:
%       |
%       | 1
% 4 ----+---- 2
%       | 
%       | 3

% Values:
% tic
UP    = 1;
RIGHT = 2;
DOWN  = 3;
LEFT  = 4;

% dx = uint8(seg ~= image_right(seg));
% dy = uint8(seg ~= image_down(seg) );
% crack_img = dx + bitshift(dy,LEFT-1) + ...
%     bitshift(image_right(dy),RIGHT-1) + bitshift(image_down(dx),DOWN-1);

% Removed dependency on image_down and image_right, to make this more
% easily packaged:
dx = uint8(seg ~= seg(:,[2:end end]));
dy = uint8(seg ~= seg([2:end end],:));
crack_img = dx + bitshift(dy,LEFT-1) + ...
    bitshift(dy(:,[2:end end]),RIGHT-1) + bitshift(dx([2:end end],:),DOWN-1);

% fprintf('New method: %.3f seconds\n', toc);
% 
% %% Old Method %%
% tic
% [nrows, ncols ] = size(seg);
% crack_img2 = zeros(nrows, ncols, 'uint8');
% % Bits:
% UP    = 1;
% RIGHT = 2;
% DOWN  = 3;
% LEFT  = 4;
% 
% index = find(seg ~= image_right(seg));
% crack_img2(index) = bitset(crack_img2(index), UP);
% 
% index = find(image_right(seg) ~= image_downright(seg));
% crack_img2(index) = bitset(crack_img2(index), RIGHT);
% 
% index = find(image_down(seg) ~= image_downright(seg));
% crack_img2(index) = bitset(crack_img2(index), DOWN);
% 
% index = find(seg ~= image_down(seg));
% crack_img2(index) = bitset(crack_img2(index), LEFT);
% 
% fprintf('Old method: %.3f seconds\n', toc);
% 
% if(all(crack_img(:)==crack_img2(:)))
%     disp('Results agree.')
% end
% 
