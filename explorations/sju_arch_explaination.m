clc

addpath(genpath('./functions'));

%% Set, Read and find difference of images Files 
previous_image_file =   "D:\GMIS_DataSet\MutliFocusV2\BedHeadBoard_Cups\seq_001.tiff";
next_image_file =       "D:\GMIS_DataSet\MutliFocusV2\BedHeadBoard_Cups\seq_002.tiff";
previous = imread(previous_image_file);
next = imread(next_image_file);

%% visualize difference
figure(1);
imagesc(previous(:,:,1));
colormap('gray');
colorbar;
title("Previous")

figure(2);
imagesc(next(:,:,1));
colormap('gray');
colorbar;
title("next")
%% Shifting and Indexing
[residue, idx1, idx2] = sju_residue_dpcm_pixels_forward(previous, next);

%%  Visualize Residue
figure(3);
image(residue(:,:,1));
colormap('gray');
colorbar;
title("residue")
%% Visualize Index Images

figure(4);
imagesc(idx1(:,:,1));
colormap([0 0 0; 1 1 1]);
% h = colorbar;
% set(h, 'Ticks', [0, 1], 'TickLabels', {'0', '1'});
title("index1")

figure(5);
imagesc(idx2(:,:,1));
colormap([0 0 0; 1 1 1]);
% h = colorbar;
% set(h, 'Ticks', [0, 1], 'TickLabels', {'0', '1'});
title("index2")

%% images histogram
d = int16(previous) - int16(next);
figure(6)
imhist(d(:,:,1), 65535)

figure(7)
imhist(d(:,:,1)+128, 65535)

%% histogram of residue
figure(8)
imhist(residue(:,:,1), 255)

