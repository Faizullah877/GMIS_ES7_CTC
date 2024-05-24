function report = encode_decode_JPEG1(input)

report = struct();
no_of_images = input.no_of_images;
encoder_exe = fullfile(input.codec_folder, "JPEG1", "cjpeg.exe");
decoder_exe = fullfile(input.codec_folder, "JPEG1", "djpeg.exe");


raw_files = dir(fullfile(input.input_raw, "*.ppm"));

ImagesResult_T1 = table('Size', [no_of_images 8], 'VariableTypes', {'uint16', 'string','uint16', 'uint16', 'double', 'single', 'double','double'});
ImagesResult_T1.Properties.VariableNames = {'ImageNo','ImageName', 'ImageWidth', 'ImageHeight', 'CompressedImageSize','bpp', 'enc_time', 'dec_time'};
q_v = input.q_value;

for imgNo = 1: no_of_images

    input_raw = fullfile(input.input_raw, raw_files(imgNo).name);
    Orig_IMG = imread(input_raw);
    [no_rows, no_cols, ~] = size(Orig_IMG);
    [~, image_name, ~] = fileparts(input_raw);
    jpg_image_file = fullfile(input.compressed_folder, sprintf("%s.jpg", image_name));

    %encoding
    if strcmp(input.config  , "arith")
        command =sprintf('%s -quality %d -arithmetic -outfile %s %s', encoder_exe,  q_v, jpg_image_file, input_raw);
    else
        command =sprintf('%s -quality %d -outfile %s %s', encoder_exe,  q_v, jpg_image_file, input_raw);
    end
    tCstartimg= tic;
    [status, ~] = system(command);
    img_enc_time  = toc(tCstartimg);
    if(status)
        fprintf("Failed to compress [%s] with [JPEG 1] codec.\n", input_raw);
        error('compression failed');
    end
    % decoding
    decoded_image_ppm = fullfile(input.decompressed_folder, sprintf("dec_image_%03d.ppm", imgNo));
    command=sprintf('%s -outfile %s %s',decoder_exe, decoded_image_ppm, jpg_image_file);
    tDstartimg = tic;
    [status, ~] = system(command);
    if(status)
        fprintf("Failed to decompress [%s] with [JPEG 1] codec.\n", jpg_image_file);
        error('decompression failed');
    end
    img_dec_time = toc(tDstartimg);

    % get bpp and encoding time
    jpg_image_file_info=dir(jpg_image_file);
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

