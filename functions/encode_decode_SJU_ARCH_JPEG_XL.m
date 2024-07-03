function report = encode_decode_SJU_ARCH_JPEG_XL(input)
report = struct();
no_of_images = input.no_of_images;
no_cols = input.width;
no_rows = input.height;

encoder_exe = fullfile(input.codec_folder, "JPEG_XL", "cjxl.exe");
decoder_exe = fullfile(input.codec_folder, "JPEG_XL", "djxl.exe");
brotli_exe = fullfile(input.codec_folder, "JPEG_XL\brotli.exe");

raw_files = dir(fullfile(input.input_raw, "*.ppm"));

ImagesResult_T1 = table('Size', [no_of_images 8], 'VariableTypes', {'uint16', 'string','uint16', 'uint16', 'double', 'single', 'double','double'});
ImagesResult_T1.Properties.VariableNames = {'ImageNo','ImageName', 'ImageWidth', 'ImageHeight', 'CompressedImageSize','bpp', 'enc_time', 'dec_time'};
q_v = input.q_value;

tmp_folder = fullfile(input.set_q_folder, "tmp");
if(~exist(tmp_folder, 'dir')), mkdir(tmp_folder); end

% Mother Image
mother_image_file = fullfile(input.input_raw, raw_files(1).name);
[~, image_name, ~] = fileparts(mother_image_file);
mother_image_compressed_file = fullfile(input.compressed_folder, sprintf("%s.jxl", image_name));
%encoding
command =sprintf('%s %s %s -q %d -e 6 --num_threads 8', encoder_exe, mother_image_file, mother_image_compressed_file, q_v);

tCstartimg= tic;
[status, ~] = system(command);
img_enc_time  = toc(tCstartimg);
if(status)
    fprintf("Failed to compress [%s] with [JPEG XL] codec.\n", mother_image_file);
    error('compression failed');
end
% decoding
decoded_image_ppm = fullfile(input.decompressed_folder, sprintf("dec_image_%03d.ppm", 1));
command = sprintf('%s %s %s --num_threads 8',decoder_exe, mother_image_compressed_file, decoded_image_ppm);
    
tDstartimg = tic;
[status, ~] = system(command);
if(status)
    fprintf("Failed to decompress [%s] with [JPEG XL] codec.\n", mother_image_compressed_file);
    error('decompression failed');
end
img_dec_time = toc(tDstartimg);

% get bpp and encoding time
jpg_image_file_info=dir(mother_image_compressed_file);
compressed_img_size = jpg_image_file_info.bytes;
img_bpp = (compressed_img_size*8) / (no_cols*no_rows);
ImagesResult_T1(1,:) ={1, image_name, no_cols, no_rows, compressed_img_size, img_bpp, img_enc_time, img_dec_time};

Orig_IMG = imread(mother_image_file);
decoded_mother_img = imread(decoded_image_ppm);
psnrv = psnr(uint8(decoded_mother_img), Orig_IMG);
fprintf("%s  => BPP   %f  => PSNR = %f\n", image_name, img_bpp, psnrv);
previous_image = int16(decoded_mother_img);

for imgNo = 2: no_of_images

    input_raw = fullfile(input.input_raw, raw_files(imgNo).name);
    Orig_IMG = imread(input_raw);
    [no_rows, no_cols, ~] = size(Orig_IMG);
    [~, image_name, ~] = fileparts(input_raw);

    next_image = int16(Orig_IMG);

    [diff_img, idx1, idx2] = sju_residue_dpcm_pixels_forward(previous_image, next_image);

    % diff_img = previous_image - next_image;
    % diff_img = diff_img + 128;
    % idx1 = diff_img>255;
    % idx2 = diff_img<0;
    % diff_img(idx1) = 511 - diff_img(idx1);
    % diff_img(idx2) =   0 - diff_img(idx2);

    idx1_tmp_file = fullfile(tmp_folder, sprintf("%s_idx1.mat", image_name));
    idx2_tmp_file = fullfile(tmp_folder, sprintf("%s_idx2.mat", image_name));
    save(idx1_tmp_file, "idx1");
    save(idx2_tmp_file, "idx2");
    residue_tmp_file = fullfile(tmp_folder, sprintf("residue_%03d.ppm", imgNo));
    imwrite(uint8(diff_img), residue_tmp_file);

    idx1_brotli_file = fullfile(input.compressed_folder, sprintf("%s_idx1.br", image_name));
    idx2_brotli_file = fullfile(input.compressed_folder, sprintf("%s_idx2.br", image_name));
    brotli_cmd = sprintf("%s -f -o %s %s", brotli_exe, idx1_brotli_file, idx1_tmp_file);
    [status1 out1] = system(brotli_cmd);
    if(status1)
        fprintf(out1);
        error('compression failed');
    end
    brotli_cmd = sprintf("%s -f -o %s %s", brotli_exe, idx2_brotli_file, idx2_tmp_file);
    [status1 out1] = system(brotli_cmd);
    if(status1)
        fprintf(out1);
        error('compression failed');
    end

    jxl_image_file = fullfile(input.compressed_folder, sprintf("%s_residue.jxl", image_name));

    %encoding
    command = sprintf('%s %s %s -q %d -e 6 --num_threads 8', encoder_exe, residue_tmp_file, jxl_image_file, q_v);
                    
    tCstartimg= tic;
    [status, ~] = system(command);
    img_enc_time  = toc(tCstartimg);
    if(status)
        fprintf("Failed to compress [%s] with [JPEG XL] codec.\n", residue_tmp_file);
        error('compression failed');
    end
    % decoding
    decoded_residue_image_ppm = fullfile(tmp_folder, sprintf("dec_image_%03d.ppm", imgNo));
    
    command=sprintf('%s %s %s --num_threads 0',decoder_exe, jxl_image_file, decoded_residue_image_ppm);
    tDstartimg = tic;
    [status, ~] = system(command);
    if(status)
        fprintf("Failed to decompress [%s] with [JPEG XL] codec.\n", jxl_image_file);
        error('decompression failed');
    end
    img_dec_time = toc(tDstartimg);

    recon_residue = imread(decoded_residue_image_ppm);
    
    rec_idx1_tmp_file = fullfile(tmp_folder, sprintf("rec_%s_idx1.mat", image_name));
    rec_idx2_tmp_file = fullfile(tmp_folder, sprintf("rec_%s_idx2.mat", image_name));
    brotli_cmd = sprintf("%s -d -f -o %s %s", brotli_exe, rec_idx1_tmp_file, idx1_brotli_file);
    [status1 out1] = system(brotli_cmd);
    if(status1)
        fprintf(out1);
        error('brotli decompression failed');
    end
    brotli_cmd = sprintf("%s -d -f -o %s %s", brotli_exe, rec_idx2_tmp_file, idx2_brotli_file);
    [status1 out1] = system(brotli_cmd);
    if(status1)
        fprintf(out1);
        error('brotli decompression failed');
    end
    rec_idx1 = load(rec_idx1_tmp_file);
    rec_idx1 = rec_idx1.idx1;
    rec_idx2 = load(rec_idx2_tmp_file);
    rec_idx2 = rec_idx2.idx2;
    rec_next_image = sju_residue_dpcm_pixels_inverse(previous_image, int16(recon_residue), rec_idx1, rec_idx2);
    
    % rec_diff_image = int16(recon_residue);
    % rec_diff_image(rec_idx1) = 511 - rec_diff_image(rec_idx1);
    % rec_diff_image(rec_idx2) =   0 - rec_diff_image(rec_idx2);
    % rec_diff_image = rec_diff_image - 128;
    % rec_next_image = previous_image - rec_diff_image;

    previous_image = int16(rec_next_image);
    decoded_image_file = fullfile(input.decompressed_folder, sprintf("dec_image_%03d.ppm", imgNo));
    imwrite(uint8(rec_next_image), decoded_image_file);


    psnrv = psnr(uint8(rec_next_image), Orig_IMG);


    % get bpp and encoding time
    br1_file_info = dir(idx1_brotli_file);
    br1_bytes = br1_file_info.bytes;
    br2_file_info = dir(idx2_brotli_file);
    br2_bytes = br2_file_info.bytes;
    jpg_image_file_info=dir(jxl_image_file);
    compressed_img_size = jpg_image_file_info.bytes + br1_bytes + br2_bytes;
    img_bpp = (compressed_img_size*8) / (no_cols*no_rows);

    fprintf("%s  => BPP   %f  => PSNR = %f\n", image_name, img_bpp, psnrv);
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

