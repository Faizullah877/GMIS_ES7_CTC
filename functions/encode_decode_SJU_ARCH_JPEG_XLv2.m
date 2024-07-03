function report = encode_decode_SJU_ARCH_JPEG_XLv2(input)
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
    diff_img = previous_image - next_image;
    diff_img = diff_img + 256;
    diff_img = uint16(diff_img);
    diff_img = single(single(diff_img)/512);


    residue_tmp_file = fullfile(tmp_folder, sprintf("residue_%03d.pfm", imgNo));
    write_pfm(diff_img, residue_tmp_file);
    % imwrite(diff_img, residue_tmp_file);

    jxl_image_file = fullfile(input.compressed_folder, sprintf("%s_residue.jxl", image_name));

    %encoding
    command = sprintf('%s %s %s -q %d -e 6 --num_threads 8', encoder_exe, residue_tmp_file, jxl_image_file, q_v);
                    
    tCstartimg= tic;
    [status, out] = system(command);
    img_enc_time  = toc(tCstartimg);
    % if(status)
    %     fprintf(out);
    %     fprintf("Failed to compress [%s] with [JPEG XL] codec.\n", residue_tmp_file);
    %     error('compression failed');
    % end
    % decoding
    decoded_residue_image_ppm = fullfile(tmp_folder, sprintf("dec_image_%03d.pfm", imgNo));
    
    command=sprintf('%s %s %s --num_threads 0',decoder_exe, jxl_image_file, decoded_residue_image_ppm);
    tDstartimg = tic;
    [status, out2] = system(command);
    if(status)
        fprintf("Failed to decompress [%s] with [JPEG XL] codec.\n", jxl_image_file);
        error('decompression failed');
    end
    img_dec_time = toc(tDstartimg);

    recon_residue = imread(decoded_residue_image_ppm);
    
    rec_diff = int16(recon_residue) - 256;
    rec_next_image = previous_image - rec_diff;

    
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
    jpg_image_file_info=dir(jxl_image_file);
    compressed_img_size = jpg_image_file_info.bytes;
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

