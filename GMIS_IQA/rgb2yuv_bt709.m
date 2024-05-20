function output_yuv = rgb2yuv_bt709(input_rgb)
%RGB2YUV_BT709 Summary of this function goes here
%   Detailed explanation goes here
% code is obtained and modified from https://www.mathworks.com/matlabcentral/fileexchange/47786-rgb-to-yuv-format?s_tid=answers_rc2-2_p5_MLT
% Santhana Raj (2023). RGB to YUV format (https://www.mathworks.com/matlabcentral/fileexchange/47786-rgb-to-yuv-format), MATLAB Central File Exchange. Retrieved March 11, 2023.

T = [0.1826,    0.6142,    0.0620;   -0.1006,   -0.3386,    0.4392;    0.4392,   -0.3989,   -0.0403];
yuvoffset = [16; 128; 128];
R = double(input_rgb(:,:,1));
G = double(input_rgb(:,:,2));
B = double(input_rgb(:,:,3));
Y = T(1,1) * R + T(1,2) * G + T(1,3) * B + yuvoffset(1);
U = T(2,1) * R + T(2,2) * G + T(2,3) * B + yuvoffset(2);
V = T(3,1) * R + T(3,2) * G + T(3,3) * B + yuvoffset(3);
Y = uint8(round(Y));
U = uint8(round(U));
V = uint8(round(V));
output_yuv = cat(3, Y, U, V);
end

