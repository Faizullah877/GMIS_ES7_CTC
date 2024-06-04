function report = encode_decode_DPCM_QDCT(input)
%ENCODE_DECODE_DPCM_QDCT Summary of this function goes here
%   Detailed explanation goes here

report = struct();
set_name = input.set_name;
q_v = input.q_value;
no_cols = input.width;
no_rows = input.height;
no_of_images = input.no_of_images;

encoding_scheme = "jpg1_dpcm_dct_seq";


encoder_exe = fullfile(input.codec_folder, "SJU_ARCH/gmis.exe");
decoder_exe = encoder_exe;  %fullfile(input.codec_folder, "gmis.exe");


ImagesResult_T1 = table('Size', [no_of_images 8], 'VariableTypes', {'uint16', 'string','uint16', 'uint16', 'double', 'single', 'double','double'});
ImagesResult_T1.Properties.VariableNames = {'ImageNo','ImageName', 'ImageWidth', 'ImageHeight', 'CompressedImageSize','bpp', 'enc_time', 'dec_time'};


gmis_file = fullfile(input.compressed_folder, sprintf("%s.jpg", set_name));
command = sprintf('%s -encode -i_yuv_dir %s -o %s -q %d -encoding_scheme %s -w %d -h %d -noi %d -pix_fmt YUV444P -en_cmg',...
    encoder_exe, input.input_raw  ,gmis_file, q_v, encoding_scheme, no_cols, no_rows, no_of_images);
% encoding to one file using dpcmed_jpeg
tCstartimg= tic;
[status, out1] = system(command);
report.set_enc_time = toc(tCstartimg);
% disp(['               => Set Compressed in ', num2str(set_enc_time), ' sec']);

command = sprintf('%s -decode %s -o_png_dir %s',decoder_exe, gmis_file, input.decompressed_folder);
tDstartimg = tic;
[status, out] = system(command);
report.set_dec_time = toc(tDstartimg);
% disp(['               => Set Decompressed in ', num2str(set_dec_time), ' sec']);

report.set_dec_folder = input.decompressed_folder;
dpcme_jpeg_file_info=dir(gmis_file);
compressed_set_size = dpcme_jpeg_file_info.bytes;
report.compressed_file = gmis_file;
report.compressed_file_size = compressed_set_size;
report.set_bpp = (compressed_set_size*8) / (no_cols*no_rows*no_of_images);
for imgNo = 1: no_of_images
    img_bpp = get_gmis_image_bpp(out1, imgNo);
    compressed_img_size = ceil((img_bpp*no_cols*no_rows)/8);
    img_enc_time = report.set_enc_time/no_of_images;
    img_dec_time = report.set_dec_time/no_of_images;
    image_name = sprintf("image%03d", imgNo);
    ImagesResult_T1(imgNo,:) ={imgNo, image_name,  uint16(no_cols), uint16(no_rows),  compressed_img_size, img_bpp, img_enc_time, img_dec_time};
end

report.images_result_table = ImagesResult_T1;
end

