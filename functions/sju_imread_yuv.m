function image = sju_imread_yuv(filename, w, h)
%SJU_IMREAD_YUV Summary of this function goes here
%   Detailed explanation goes here
%   filename = filename of .yuv input file to read
                % the file should contain the pixels in planar format
%   w = width = number of columns of input image
%   h = height = number of rows of input image
%   image is the output of this fuction.
fid=fopen(filename,'r');
image= uint8(fread(fid));
image= reshape(image, [w h 3]);
fclose(fid);
image= permute(image, [2, 1, 3]);
end

