function report = encode_decode_JPEG_XL(input)

report = struct();
no_of_images = input.no_of_images;

encoder_exe = fullfile(input.codec_folder, "JPEG_XL", "cjxl.exe");
decoder_exe = fullfile(input.codec_folder, "JPEG_XL", "djxl.exe");

raw_files = dir(fullfile(input.input_raw, input.input_raw_fileEXT));

ImagesResult_T1 = table('Size', [no_of_images 8], 'VariableTypes', {'uint16', 'string','uint16', 'uint16', 'double', 'single', 'double','double'});
ImagesResult_T1.Properties.VariableNames = {'ImageNo','ImageName', 'ImageWidth', 'ImageHeight', 'CompressedImageSize','bpp', 'enc_time', 'dec_time'};
q_v = input.q_value;

for imgNo = 1: no_of_images

    input_raw = fullfile(input.input_raw, raw_files(imgNo).name);
    Orig_IMG = imread(input_raw);
    [no_rows, no_cols, ~] = size(Orig_IMG);
    [~, image_name, ~] = fileparts(input_raw);
    jxl_image_file = fullfile(input.compressed_folder, sprintf("%s.jxl", image_name));

    %encoding
    command =sprintf('%s %s %s -q %d -e 6 --num_threads 0', encoder_exe, input_raw, jxl_image_file, q_v);
                
    tCstartimg= tic;
    [status, ~] = system(command);
    img_enc_time  = toc(tCstartimg);
    if(status)
        fprintf("Failed to compress [%s] with [JPEG XL] codec.\n", input_raw);
        error('compression failed');
    end
    % decoding
    decoded_image_ppm = fullfile(input.decompressed_folder, sprintf("dec_image_%03d.ppm", imgNo));
    command=sprintf('%s %s %s --num_threads 0',decoder_exe, jxl_image_file, decoded_image_ppm);
    tDstartimg = tic;
    [status, ~] = system(command);
    if(status)
        fprintf("Failed to decompress [%s] with [JPEG XL] codec.\n", jxl_image_file);
        error('decompression failed');
    end
    img_dec_time = toc(tDstartimg);

    % get bpp and encoding time
    jpg_image_file_info=dir(jxl_image_file);
    compressed_img_size = jpg_image_file_info.bytes;
    img_bpp = (compressed_img_size*8) / (no_cols*no_rows);
    ImagesResult_T1(imgNo,:) ={imgNo, image_name, no_cols, no_rows, compressed_img_size, img_bpp, img_enc_time, img_dec_time};
end


report.set_enc_time = sum(ImagesResult_T1.("enc_time"));
compressed_set_size = sum(ImagesResult_T1.("CompressedImageSize"));
report.compressed_file_size = compressed_set_size;
report.set_bpp = (compressed_set_size*8) / (input.width*input.height*no_of_images);
report.set_dec_time = sum(ImagesResult_T1.("dec_time"));
report.set_dec_folder = input.decompressed_folder;
report.images_result_table = ImagesResult_T1;

end

