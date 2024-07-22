function report = encode_decode_JPEG2000(input)

report = struct();
set_info = input.set_info;
no_of_images = set_info.no_of_images;

encoder_exe = fullfile(input.codec_folder, "JPEG_2000", "kdu_compress.exe");
decoder_exe = fullfile(input.codec_folder, "JPEG_2000", "kdu_expand.exe");

ImagesResult = table('Size', [no_of_images 8], 'VariableTypes', {'uint16', 'string','uint16', 'uint16', 'double', 'single', 'double','double'});
ImagesResult.Properties.VariableNames = {'ImageNo','ImageName', 'ImageWidth', 'ImageHeight', 'CompressedImageSize','bpp', 'enc_time', 'dec_time'};
q_v = input.q_value;

for imgNo = 1: no_of_images 
    image_name_wo_ext = set_info.images_details.NameWoExt(imgNo);
    input_raw = fullfile(input.input_raw, sprintf("%s.ppm", image_name_wo_ext));
    j2c_image_file = fullfile(input.compressed_folder, sprintf("%s.j2c", image_name_wo_ext));

    %encoding
    command = sprintf("%s -i %s -o %s Qfactor=%d", encoder_exe, input_raw, j2c_image_file, q_v);
                
    tCstartimg= tic;
    [status, ~] = system(command);
    img_enc_time  = toc(tCstartimg);
    if(status)
        fprintf("Failed to compress [%s] with [JPEG 2000] codec.\n", input_raw);
        error('compression failed');
    end
    % decoding
    decoded_image_ppm = fullfile(input.decompressed_folder, sprintf("%s.ppm", image_name_wo_ext));

    command = sprintf("%s -i %s -o %s", decoder_exe, j2c_image_file, decoded_image_ppm);
    tDstartimg = tic;
    [status, out2] = system(command);
    if(status)
        fprintf("Failed to decompress [%s] with [JPEG 2000] codec.\n", jpg_image_file);
        disp(out2);
        error('decompression failed');
    end
    img_dec_time = toc(tDstartimg);

    % get bpp and encoding time
    jpg_image_file_info=dir(j2c_image_file);
    compressed_img_size = jpg_image_file_info.bytes;
    no_rows = set_info.images_details.ImageHeight(imgNo);
    no_cols = set_info.images_details.ImageWidth(imgNo);
    img_bpp = (compressed_img_size*8) / (double(no_cols)*double(no_rows));
    ImagesResult(imgNo,:) ={imgNo, image_name_wo_ext, no_cols, no_rows, compressed_img_size, img_bpp, img_enc_time, img_dec_time};
end


report.set_enc_time = sum(ImagesResult.("enc_time"));
compressed_set_size = sum(ImagesResult.("CompressedImageSize"));
report.compressed_set_size = compressed_set_size;
report.set_bpp = (compressed_set_size*8) / double(set_info.num_pixels);
report.set_dec_time = sum(ImagesResult.("dec_time"));
report.set_dec_folder = input.decompressed_folder;
report.images_result_table = ImagesResult;

end

