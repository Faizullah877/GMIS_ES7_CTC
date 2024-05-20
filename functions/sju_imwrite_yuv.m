function  sju_imwrite_yuv(Image, fileName)
% Image = a 3 channel image
% fileName = path + filename.yuv to write the .yuv file.
fid=fopen(fileName,'w');
Image = permute(Image,[2,1,3]);
fwrite(fid, Image, 'uint8'); % writes in 444 planar format
fclose(fid);
end

