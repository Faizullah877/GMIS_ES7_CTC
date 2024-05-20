% clc
clear all
% https://ieeexplore.ieee.org/document/9455944

set_name = "TestPlant_Teapot";
set_path = "D:\GMIS_EXPs\Test_CTC\dataset";

test_dir = "D:\GMIS_EXPs\Test_CTC2";

if(~exist(test_dir, 'dir')), mkdir(test_dir); end

codecs_exe_dir = "D:\GMIS_ES_7\codecs";
%% Conversion to YUV420 video 

tCstart = tic;
orig_ext = "*.tiff";
original_images_list = dir(fullfile(set_path,set_name, orig_ext));

no_of_images = numel(original_images_list);
codec_in_file_format = "yuv420p10le";
input_file = fullfile(original_images_list(1).folder, original_images_list(1).name);

[~, image_name, ~] = fileparts(original_images_list(1).name);
input =fullfile(set_path, set_name, erase(image_name,"001")+"%03d.tiff");

rgb_img = imread(input_file);
[he, we, ce] = size(rgb_img);
orig_raw_yuv_fname = sprintf("%s_%dx%dx%d_1Hz_prec10_P420.yuv", set_name, we,he, no_of_images);
orig_raw_yuv = fullfile(test_dir, orig_raw_yuv_fname);
cmd = sprintf("ffmpeg -hide_banner -i %s -pix_fmt %s -vf scale=in_range=full:in_color_matrix=bt709:out_range=full:out_color_matrix=bt709 -color_primaries bt709  -color_trc bt709 -colorspace bt709 -y %s", input, codec_in_file_format, orig_raw_yuv);
[status, out1] = system(cmd);

conversion_time = toc(tCstart);
fprintf(['   Conversion done in  ' num2str(conversion_time) ' sec\n']);

%% VVC Encoding.
config = "intra"; % intra or inter

tEncStart = tic;
encoder_exe = fullfile(codecs_exe_dir, "VVC_VTM/EncoderApp.exe");


q_v = 35;

compressed_file = fullfile(test_dir, sprintf("%s_q%03d.266", set_name, q_v));
de_compressed_file = fullfile(test_dir, sprintf("dec_%s", orig_raw_yuv_fname));

if strcmp(config, "inter")
    cfg_file = fullfile(codecs_exe_dir, "VVC_VTM","cfg", "encoder_randomaccess_vtm.cfg");
    command =sprintf('%s -c %s -i %s -wdt %d -hgt %d -b %s -f %d -fr 30 -q %d --InputBitDepth=10 --InputChromaFormat=420 --ChromaFormatIDC=420 --TemporalSubsampleRatio=1 --ConformanceWindowMode=1 --MinSearchWindow=1 -o %s', encoder_exe, cfg_file, orig_raw_yuv, we, he, compressed_file, no_of_images, q_v, de_compressed_file);
    % command =sprintf('%s -c %s -i %s -wdt %d -hgt %d -b %s -f %d -fr 30 -q %d --InputBitDepth=10 --InputChromaFormat=420 --ChromaFormatIDC=420 --TemporalSubsampleRatio=1 --ConformanceWindowMode=1 --MinSearchWindow=1', encoder_exe, cfg_file, orig_raw_yuv, we, he, compressed_file, no_of_images, q_v);
else
    cfg_file = fullfile(codecs_exe_dir, "VVC_VTM","cfg", "encoder_intra_vtm.cfg");
    command =sprintf('%s -c %s -i %s -wdt %d -hgt %d -b %s -f %d -fr 30 -q %d --InputBitDepth=10 --InputChromaFormat=420 --ChromaFormatIDC=420 --TemporalSubsampleRatio=1 --ConformanceWindowMode=1 --IntraPeriod=1 -o %s', encoder_exe, cfg_file, orig_raw_yuv, we, he, compressed_file, no_of_images, q_v, de_compressed_file);
    % command =sprintf('%s -c %s -i %s -wdt %d -hgt %d -b %s -f %d -fr 30 -q %d --InputBitDepth=10 --InputChromaFormat=420 --ChromaFormatIDC=420 --TemporalSubsampleRatio=1 --ConformanceWindowMode=1 --IntraPeriod=1', encoder_exe, cfg_file, orig_raw_yuv, we, he, compressed_file, no_of_images, q_v);
end

[status, out2] = system(command)

compression_time = toc(tEncStart);
fprintf(['   Compression done in  ' num2str(compression_time) ' sec\n']);


% %% VVC Decoding
% tDecStart = tic;
% dec_exe = fullfile(codecs_exe_dir, "VVC_VVenC/vvdecapp.exe");
% 
% de_compressed_file = fullfile(test_dir, sprintf("dec_%s", orig_raw_yuv_fname));
% 
% dec_cmd = sprintf("%s -b %s -o %s ", dec_exe,  compressed_file, de_compressed_file);
% 
% [status, out3] = system(dec_cmd)
% 
% decode_time = toc(tDecStart);
% fprintf(['   Decoding done in  ' num2str(decode_time) ' sec\n']);

%% back Conversion to images
tBackConvStart = tic;
decoded_images_dir = fullfile(test_dir, sprintf("Decoded_%d", q_v));
if(~exist(decoded_images_dir, 'dir')), mkdir(decoded_images_dir); end
decoded_output = fullfile(decoded_images_dir, sprintf("dec_%s_q%d_%%03d.png",set_name, q_v));

reconv_cmd = sprintf("ffmpeg -hide_banner -f rawvideo -s %dx%d -pix_fmt yuv420p10le -i %s   -y %s", we, he, de_compressed_file, decoded_output);
[status, out1] = system(reconv_cmd);

deconversion_time = toc(tBackConvStart);
fprintf(['   Re-Conversion done in  ' num2str(deconversion_time) ' sec\n']);


%% evaluation
