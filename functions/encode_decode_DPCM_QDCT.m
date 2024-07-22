function report = encode_decode_DPCM_QDCT(input)
report = struct();

set_info = input.set_info;
set_name = set_info.set_name;
if ~set_info.same_resolution
    fprintf("The codec {gmis.exe} does not support sets {%s} with images having different resolution.", set_name);
    error("Compression Failed");
end

q_v = input.q_value;
no_of_images = set_info.no_of_images;

encoding_scheme = "jpeg1_dpcm_qdct_seq";


encoder_exe = fullfile(input.codec_folder, "SJU_ARCH/gmis.exe");
decoder_exe = encoder_exe;  %fullfile(input.codec_folder, "gmis.exe");

no_rows = set_info.images_details.ImageHeight(1);
no_cols = set_info.images_details.ImageWidth(1);

ImagesResult = table('Size', [no_of_images 8], 'VariableTypes', {'uint16', 'string','uint16', 'uint16', 'double', 'single', 'double','double'});
ImagesResult.Properties.VariableNames = {'ImageNo','ImageName', 'ImageWidth', 'ImageHeight', 'CompressedImageSize','bpp', 'enc_time', 'dec_time'};


gmis_file = fullfile(input.compressed_folder, sprintf("%s.jpg", set_name));
command = sprintf('%s -encode -i_yuv_dir %s -o %s -q %d -encoding_scheme %s -w %d -h %d -noi %d -pix_fmt YUV444P -en_cmg',...
    encoder_exe, input.input_raw  ,gmis_file, q_v, encoding_scheme, no_cols, no_rows, no_of_images);
% encoding to one file using dpcmed_jpeg
tCstartimg= tic;
[status, out1] = system(command);
report.set_enc_time = toc(tCstartimg);
if(status)
    fprintf("Failed to compress [%s] with [gmis.exe] codec.\n", input.input_raw);
    disp(out1);
    error('compression failed');
end
command = sprintf('%s -decode %s -o_ppm_dir %s',decoder_exe, gmis_file, input.decompressed_folder);
tDstartimg = tic;
[status, out2] = system(command);
report.set_dec_time = toc(tDstartimg);
if(status)
    fprintf("Failed to decompress [%s] with [gmis.exe] codec.\n", gmis_file);
    disp(out2);
    error('decompression failed');
end
report.set_dec_folder = input.decompressed_folder;
dpcme_jpeg_file_info=dir(gmis_file);
compressed_set_size = dpcme_jpeg_file_info.bytes;
report.compressed_file = gmis_file;
for imgNo = 1: no_of_images
    image_name_wo_ext = set_info.images_details.NameWoExt(imgNo);
    img_bpp = get_gmis_image_bpp(out1, imgNo);
    compressed_img_size = ceil((img_bpp*no_cols*no_rows)/8);
    img_enc_time = report.set_enc_time/no_of_images;
    img_dec_time = report.set_dec_time/no_of_images;
    ImagesResult(imgNo,:) ={imgNo, image_name_wo_ext,  uint16(no_cols), uint16(no_rows),  compressed_img_size, img_bpp, img_enc_time, img_dec_time};
end

report.set_enc_time = sum(ImagesResult.("enc_time"));
report.compressed_set_size = compressed_set_size;
report.set_bpp = (compressed_set_size*8) / double(set_info.num_pixels);
report.set_dec_time = sum(ImagesResult.("dec_time"));
report.set_dec_folder = input.decompressed_folder;
report.images_result_table = ImagesResult;
end

