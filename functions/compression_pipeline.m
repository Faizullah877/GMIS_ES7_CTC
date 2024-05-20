function M = compression_pipeline(exp_params, codec_params, M)
%COMPRESSION_PIPELINE Summary of this function goes here
%   Detailed explanation goes here
set_name = exp_params.set_name;
QA_metrics = exp_params.QA_Metrics;
do_discard = exp_params.do_discard;

codec_exp_time = tic;

codec = codec_params{1};
q_list = codec_params{2};
codec_in_file_format = codec_params{3};



quality_count = numel(q_list);
%% Data conversion
if(strcmp(codec, "VVC_VVenC_Inter") || strcmp(codec, "VVC_VVenC_Intra") ||codec == "VVC_VTM_Inter" || codec == "VVC_VTM_Intra")
    exp_params = convertDataset_to_video(exp_params, codec_in_file_format);
else
    exp_params = convert_dataset(exp_params, codec_in_file_format);
end

dir_struct = exp_params.dir_struct;

no_of_images = exp_params.no_of_images;

set_codec_results = struct();
codec_data_folder = fullfile(dir_struct.set_data_folder, codec);
if(~exist(codec_data_folder, 'dir')), mkdir(codec_data_folder); end

orig_files = dir(fullfile(exp_params.set_path, exp_params.images_files_ext));

%% Iterate over rate/q values
for q_itr = 1: quality_count
    q_v = q_list(q_itr);
    q_fieldName = ['q_' num2str(q_v)];
    q_fieldName = strrep(q_fieldName, ".", "p");
    if isfield(M, set_name)
        if isfield(M.(set_name), codec)
            if isfield(M.(set_name).(codec).q_fields, q_fieldName)
                str = sprintf("\n\tSet Name : [%s] => Codec : [%s] => Q field : [%s] already exists and the experiment is skipped.\n", set_name, codec, q_fieldName);
                fprintf(str);
                continue
            end
        end
    end

    set_q_folder = fullfile(codec_data_folder, q_fieldName);
    if(~exist(set_q_folder, 'dir')), mkdir(set_q_folder); end

    set_compressed_folder = fullfile(set_q_folder, "Compressed");
    if(~exist(set_compressed_folder, 'dir')), mkdir(set_compressed_folder); end
    set_dec_folder = fullfile(set_q_folder, "Decompressed");
    if(~exist(set_dec_folder, 'dir')), mkdir(set_dec_folder); end

    ImagesResult_T1 = table('Size', [no_of_images 12], 'VariableTypes', {'string', 'uint16', 'string','uint16', 'uint16', 'double', 'string', 'single', 'double', 'double','double','double'});
    ImagesResult_T1.Properties.VariableNames = {'set_name', 'ImageNo','ImageName', 'ImageWidth', 'ImageHeight', 'PixelsCount' ,'codec', 'rate', 'CompressedImageSize','bpp', 'enc_time', 'dec_time'};
    Images_IQA_T = table();

    print_current_time();
    fprintf("\n\tCompression Pipeline running for SET : [%s]   CODEC : [%s] and RATE : [%3.2f] => %d/%d\n",set_name, codec, q_v, q_itr, quality_count);
    q_start_time = tic();

    if(codec == "H264_AVC" || codec == "H264_AVCintra")
        %% H264 video codec
        input_orig = fullfile(orig_files(1).folder, orig_files(1).name);
        Orig_IMG = imread(input_orig);
        [no_rows, no_cols, ~] = size(Orig_IMG);
        encoder_exe=fullfile(dir_struct.codec_folder,'ffmpeg.exe');
        decoder_exe = fullfile(dir_struct.codec_folder, 'ffmpeg.exe');
        [~, image_name, ~] = fileparts(orig_files(1).name);
        ppm_input =fullfile(dir_struct.set_ppm_path, erase(image_name,"001")+"%03d.ppm");
        mp4_file = fullfile(set_compressed_folder, sprintf('%s_q%d.264',set_name,q_v ));
        ppm_output=fullfile(set_dec_folder,"dec_image_%03d.ppm");
        if(codec == "H264_AVC")
            enc_command = sprintf('%s -y -s %dx%d -r %d -i %s -bf 1 -c:v libx264 -crf %d  %s',encoder_exe, no_cols, no_rows, 1, ppm_input, q_v, mp4_file); % -y is to overwrite if output file already exist.
        else
            enc_command = sprintf('%s -y -s %dx%d -r %d -i %s -g 1 -c:v libx264 -crf %d  %s',encoder_exe, no_cols, no_rows, 1, ppm_input, q_v, mp4_file); % -y is to overwrite if output file already exist.
        end
        tCstartimg= tic;
        [status, out1] = system(enc_command);
        set_enc_time = toc(tCstartimg);

        fprintf('\t\t=> Set compressed.\n');
        % decoding to multiple .ppm files;
        dec_command = sprintf ('%s -y -i %s %s', decoder_exe, mp4_file, ppm_output );
        tDstartimg = tic;
        [status, out2] = system(dec_command);
        set_dec_time = toc(tDstartimg);
        fprintf('\t\t=> Set decompressed.\n');

        mp4_file_info=dir(mp4_file);
        compressed_set_size = mp4_file_info.bytes;
        set_bpp = (mp4_file_info.bytes*8) / (no_cols*no_rows*no_of_images);
        eval_set_start_time = tic;
        textprogressbar('            => Evaluating : ');
        for imgNo=1:no_of_images
            input_raw = fullfile(orig_files(imgNo).folder, orig_files(imgNo).name);
            Orig_IMG = imread(input_raw);
            [no_rows, no_cols, ~] = size(Orig_IMG);
            [~, image_name, ~] = fileparts(input_raw);
            decoded_ppm_file = fullfile(set_dec_folder, sprintf('dec_image_%03d.ppm', imgNo));
            compressed_img_size = compressed_set_size/no_of_images;
            img_bpp = set_bpp;
            img_enc_time = set_enc_time/no_of_images;
            img_dec_time = set_dec_time/no_of_images;
            %% Evaluation
            Recon_IMG = imread(decoded_ppm_file);
            %             image_IQA_values = evaluate_images(Orig_IMG, Recon_IMG, "all");
            image_IQA_values = evaluate_images_IQA(Orig_IMG, Recon_IMG, QA_metrics);
            ImagesResult_T1(imgNo,:) ={set_name, imgNo, image_name, no_cols, no_rows, no_rows*no_cols, codec, q_v, compressed_img_size, img_bpp, img_enc_time, img_dec_time};
            Images_IQA_T = [Images_IQA_T;struct2table(image_IQA_values)];
            textprogressbar(imgNo/no_of_images*100);
        end
        eval_set_time = toc(eval_set_start_time);
        textprogressbar(['   done in  ' num2str(set_enc_time+eval_set_time+set_dec_time) ' sec']);

    elseif(codec == "SJU_Arch_JPEG1_PIX_SEQ" || codec == "SJU_Arch_DPCM_QDCT" || codec == "SJU_Arch_JPEG1_only")
        %% SJU_Arch without RDOPT

        if(codec == "SJU_Arch_JPEG1_only")
            encoding_scheme = "jpg1_only";
        elseif(codec == "SJU_Arch_JPEG1_PIX_SEQ")
            encoding_scheme = "jpg1_dpcm_pixels_seq_v2";
        elseif(codec == "SJU_Arch_DPCM_QDCT")
            encoding_scheme = "jpg1_dpcm_dct_seq";
        end
        input_orig = fullfile(orig_files(1).folder, orig_files(1).name);
        Orig_IMG = imread(input_orig);
        [no_rows, no_cols, ~] = size(Orig_IMG);
        encoder_exe=fullfile(dir_struct.codec_folder,'gmis.exe');
        decoder_exe = fullfile(dir_struct.codec_folder, 'gmis.exe');
        index_images_folders = fullfile(set_q_folder, "IndexImages");
        if(~exist(index_images_folders, 'dir')), mkdir(index_images_folders); end
        gmis_file = fullfile(set_compressed_folder, sprintf("%s_%s.jpg", set_name, q_fieldName));
        command = sprintf('%s -encode -i_yuv_dir %s -o %s -index_images_dir %s -q %d -encoding_scheme %s -w %d -h %d -noi %d -pix_fmt YUV444P -en_cmg',...
            encoder_exe, dir_struct.set_jpeg_ycbcr_path ,gmis_file,index_images_folders, q_v, encoding_scheme, no_cols, no_rows, no_of_images);
        % encoding to one file using dpcmed_jpeg
        tCstartimg= tic;
        [status, out1] = system(command);
        set_enc_time = toc(tCstartimg);
        disp(['               => Set Compressed in ', num2str(set_enc_time), ' sec']);

        command = sprintf('%s -decode %s -o_png_dir %s -index_images_dir %s',decoder_exe, gmis_file, set_dec_folder, index_images_folders);
        tDstartimg = tic;
        [status, out] = system(command);
        set_dec_time = toc(tDstartimg);
        disp(['               => Set Decompressed in ', num2str(set_dec_time), ' sec']);

        dpcme_jpeg_file_info=dir(gmis_file);
        compressed_set_size = dpcme_jpeg_file_info.bytes;
        set_bpp = (compressed_set_size*8) / (no_cols*no_rows*no_of_images);
        eval_set_start_time = tic;
        textprogressbar('            => Evaluating : ');
        for imgNo=1:no_of_images
            input_orig = fullfile(orig_files(imgNo).folder, orig_files(imgNo).name);
            Orig_IMG = imread(input_orig);
            [no_rows, no_cols, ~] = size(Orig_IMG);
            [~, image_name, ~] = fileparts(input_orig);
            decoded_ppm_file = fullfile(set_dec_folder, sprintf('dec_image_%03d.png', imgNo));
            img_bpp = get_gmis_image_bpp(out1, imgNo);
            compressed_img_size = ceil((img_bpp*no_cols*no_rows)/8);
            img_enc_time = set_enc_time/no_of_images;
            img_dec_time = set_dec_time/no_of_images;
            %% Evaluation
            Recon_IMG = imread(decoded_ppm_file);
            image_IQA_values = evaluate_images_IQA(Orig_IMG, Recon_IMG, QA_metrics);
            ImagesResult_T1(imgNo,:) ={set_name, imgNo, image_name, no_cols, no_rows, no_rows*no_cols, codec, q_v, compressed_img_size, img_bpp, img_enc_time, img_dec_time};
            Images_IQA_T = [Images_IQA_T;struct2table(image_IQA_values)];
            textprogressbar(imgNo/no_of_images*100);
        end
        eval_set_time = toc(eval_set_start_time);
        textprogressbar(['   done in  ' num2str(set_enc_time+eval_set_time+set_dec_time) ' sec']);

    elseif(codec == "SJU_Arch_DPCM_PIXs_rdopt1")
        %% SJU_Arch with RDOPT (matlab)
        encoder_exe = fullfile(dir_struct.codec_folder, "cjpeg.exe");
        decoder_exe = fullfile(dir_struct.codec_folder, "djpeg.exe");
        rdopt_exe = fullfile(dir_struct.codec_folder, "RD_Opt.exe");
        img_qfile_name = "qtable_tmp.txt";
        orig_residue_dir = fullfile(set_q_folder, "Orig_Residue");
        if(~exist(orig_residue_dir, 'dir')), mkdir(orig_residue_dir); end
        dec_residue_dir = fullfile(set_q_folder, "Dec_Residue");
        if(~exist(dec_residue_dir, 'dir')), mkdir(dec_residue_dir); end
        jpgs_folder = fullfile(set_compressed_folder, "Compressed_jpgs");
        if(~exist(jpgs_folder, 'dir')), mkdir(jpgs_folder); end
        input_raw = fullfile(dir_struct.set_ppm_path, ppm_images_names(1).name);
        mother_rgb = imread(input_raw);
        [no_rows, no_cols, ~] = size(mother_rgb);

        [~, image_name, ~] = fileparts(input_raw);
        orig_ppm_fname = fullfile(orig_residue_dir, "mother_image.ppm");
        imwrite(mother_rgb, orig_ppm_fname);

        % compress mother image with jpeg-1
        jpeg_fname = fullfile(jpgs_folder, "mother_image_001.jpg");

        if(exist(img_qfile_name, "file"))
            delete(img_qfile_name);
        end
        rdopt_cmd = sprintf("%s -numtables 2 -bpp %f -qfile %s -im %s", rdopt_exe, q_v, img_qfile_name, input_raw);
        [status, out1] = system(rdopt_cmd);
        enc_command =sprintf('%s -qtables %s -arithmetic -outfile %s %s', encoder_exe,  img_qfile_name, jpeg_fname, input_raw);

        m_enc_stime = tic;
        [abc, out] = system(enc_command);
        if(abc)
            fprintf("Failed for [%s] to compress [%s] with [JPEG 1] codec.\n", codec, input_raw);
            error('compression failed');
        end
        img_enc_time = toc(m_enc_stime);
        % get bpp and encoding time
        jpg_image_file_info=dir(jpeg_fname);
        compressed_img_size = jpg_image_file_info.bytes;
        img_bpp = (compressed_img_size*8) / (no_cols*no_rows);

        decoded_ppm_fname = fullfile(set_dec_folder, "dec_mother_image.ppm");
        dec_command = sprintf("%s -outfile %s %s", decoder_exe, decoded_ppm_fname, jpeg_fname);
        m_dec_stime = tic;
        [xyz, out2] = system(dec_command);
        if(xyz)
            fprintf("Failed for [%s] to decompress [%s] with [JPEG 1] codec.\n", codec, jpg_image_file);
            error('decompression failed');
        end
        img_dec_time = toc(m_dec_stime);

        Recon_IMG = imread(decoded_ppm_fname);
        image_IQA_values = evaluate_images_IQA(mother_rgb, Recon_IMG, QA_metrics);
        ImagesResult_T1(1,:) ={set_name, 1, image_name, no_cols, no_rows, no_rows*no_cols, codec, q_v, compressed_img_size, img_bpp, img_enc_time, img_dec_time};
        Images_IQA_T = [Images_IQA_T;struct2table(image_IQA_values)];

        previous_recon_ycbcr = rgb2jpegycbcr(Recon_IMG);

        q_v = q_v  /3;
        q_start_time = tic;
        textprogressbar('        => Progress  : ');
        for imgNo =2:no_of_images
            next_raw_file = fullfile(dir_struct.set_ppm_path, ppm_images_names(imgNo).name);
            next_image_rgb = imread(next_raw_file);
            [h, w, c] = size(next_image_rgb);
            [no_rows, no_cols, ~] = size(next_image_rgb);
            [~, image_name, ~] = fileparts(next_raw_file);
            residue_gen_Tstart = tic;
            next_image_ycbcr = rgb2jpegycbcr(next_image_rgb);
            residue_9bit = int16(previous_recon_ycbcr) - int16(next_image_ycbcr);
            residue_9bit_ppm_fname = fullfile(orig_residue_dir, sprintf("residue9b_image_%03d.ppm", imgNo));
            residue_9bit_tmp = uint16(residue_9bit+255);
            imwrite(residue_9bit_tmp, residue_9bit_ppm_fname, 'MaxValue', 65535);

            % convert 16bit residue to 8 bit residue using SJU indexing method
            [residue8bit, idx, idx_bytes] = convert_10bit_residue_to_8bit(residue_9bit, set_compressed_folder, imgNo);
            residue_gen_time =  toc(residue_gen_Tstart);

            % compress 8 bit residue with JPEG 1 in arithmetic mode
            residue_ppm_fname = fullfile(orig_residue_dir, sprintf("residue8b_image_%03d.ppm", imgNo));
            imwrite(residue8bit, residue_ppm_fname);
            residue_jpeg_fname = fullfile(jpgs_folder, sprintf("residue8b_image_%03d.jpg", imgNo));

            if(exist(img_qfile_name, "file"))
                delete(img_qfile_name);
            end
            rdopt_cmd = sprintf("%s -numtables 2 -bpp %f -qfile %s -im %s", rdopt_exe, q_v, img_qfile_name, residue_ppm_fname);
            [status, ~] = system(rdopt_cmd);
            enc_command =sprintf('%s -qtables %s -rgb -arithmetic -outfile %s %s', encoder_exe,  img_qfile_name, residue_jpeg_fname, residue_ppm_fname);


            enc_startT = tic;
            [abc, out] = system(enc_command);
            if(abc)
                fprintf("Failed for [%s] to compress [%s] with [JPEG 1] codec.\n", codec, residue_ppm_fname);
                error('compression failed');
            end
            img_enc_time = toc(enc_startT) + residue_gen_time;

            jpg_info1= dir(residue_jpeg_fname);
            compressed_img_size_a = jpg_info1.bytes;
            jpg_bpp_a = (compressed_img_size_a*8)/(w*h);
            compressed_img_size_b = compressed_img_size_a + idx_bytes;
            compressed_img_size = compressed_img_size_b;
            jpg_bpp_b = (compressed_img_size_b*8)/(w*h);
            img_bpp = (compressed_img_size *8)/(w*h);

            decoded_residue8bit_ppm_fname = fullfile(dec_residue_dir, sprintf("dec_residue8b_image_%03d.ppm", imgNo));
            dec_command = sprintf("%s -outfile %s %s", decoder_exe, decoded_residue8bit_ppm_fname, residue_jpeg_fname);
            dec_startT = tic;
            [xyz, out2] = system(dec_command);
            if(xyz)
                fprintf("Failed for [%s] to decompress [%s] with [JPEG 1] codec.\n", codec, residue_jpeg_fname);
                error('decompression failed');
            end
            img_dec_time = toc(dec_startT);

            % read decoded residue and evaluate
            decoded_residue8bit = imread(decoded_residue8bit_ppm_fname);
            psnr_a = psnr(decoded_residue8bit, residue8bit);
            mse_a = immse(decoded_residue8bit, residue8bit);

            % inverse conversion of decoded 8bit residue to 16 bit residue
            % using SJU indexing method
            recon_residue = convert_8bit_residue_to_10bit(decoded_residue8bit, idx);
            psnr_b = psnr(uint16(residue_9bit_tmp), uint16(recon_residue+255), max(residue_9bit_tmp, [], "all"));
            %             mse_b = immse(residue_16bit, recon_residue);
            %             decoded_residue16bit_ppm_fname = fullfile(dec_residue_dir, sprintf("dec_residue16b_image_%03d.ppm", imgNo));
            %             imwrite(recon_residue, decoded_residue16bit_ppm_fname);

            % inverse residue
            residue_inv_tStart = tic;
            %get reconstructed image from 16 bit residue using USNW decoder
            %             recon_next_rgb = inverse_residue_regis_prediction_v3_based(set_name, set_compressed_folder, previous_recon_rgb, next_image_rgb, decoded_residue16bit_ppm_fname, ModelID, MAX_NUM_MODELS, QTree_Depth, SMALLEST_BLK_SIZE, imgNo);
            recon_next_ycbcr = int16(previous_recon_ycbcr) - int16(recon_residue);
            recon_next_ycbcr = uint8(recon_next_ycbcr);

            residue_inv_time = toc(residue_inv_tStart);
            img_dec_time = img_dec_time+residue_inv_time;

            recon_next_rgb = jpegycbcr2rgb(recon_next_ycbcr);

            % Evaluation the reconstructed image
            %             previous_recon_rgb = recon_next_rgb;
            %             previous_orig_rgb = next_image_rgb;
            previous_recon_ycbcr = recon_next_ycbcr;
            decoded_rgb_ppf_fname = fullfile(set_dec_folder, sprintf("dec_image_%03d.ppm", imgNo));
            imwrite(recon_next_rgb, decoded_rgb_ppf_fname);

            image_IQA_values = evaluate_images_IQA(next_image_rgb, recon_next_rgb, QA_metrics);
            ImagesResult_T1(imgNo,:) ={set_name, imgNo, image_name, no_cols, no_rows, no_rows*no_cols, codec, q_v, compressed_img_size, img_bpp, img_enc_time, img_dec_time};
            Images_IQA_T = [Images_IQA_T;struct2table(image_IQA_values)];
            textprogressbar(imgNo/no_of_images*100);
        end
        q_time_taken = toc(q_start_time);
        textprogressbar(['   done in  ' num2str(q_time_taken) ' sec']);
    elseif(codec == "UNSW_SJU_JPEG1huff" || codec == "UNSW_SJU_JPEG1arith" || codec == "UNSW_SJU_JPEG1huff_rdopt" || codec == "UNSW_SJU_JPEG1arith_rdopt" || codec == "UNSW_SJU_JPEG1arith_rdopt_22p")
        %% UNSW_SJU_Arch
        ModelID = 8;
        MAX_NUM_MODELS = 3;
        QTree_Depth = 4;
        SMALLEST_BLK_SIZE = 16;
        encoder_exe = fullfile(dir_struct.codec_folder, "cjpeg.exe");
        decoder_exe = fullfile(dir_struct.codec_folder, "djpeg.exe");
        rdopt_exe = fullfile(dir_struct.codec_folder, "RD_Opt.exe");
        img_qfile_name = "qtable_tmp.txt";
        orig_residue_dir = fullfile(set_q_folder, "Orig_Residue");
        if(~exist(orig_residue_dir, 'dir')), mkdir(orig_residue_dir); end
        dec_residue_dir = fullfile(set_q_folder, "Dec_Residue");
        if(~exist(dec_residue_dir, 'dir')), mkdir(dec_residue_dir); end
        jpgs_folder = fullfile(set_compressed_folder, "Compressed_jpgs");
        if(~exist(jpgs_folder, 'dir')), mkdir(jpgs_folder); end
        input_raw = fullfile(dir_struct.set_ppm_path, ppm_images_names(1).name);
        mother_rgb = imread(input_raw);
        [no_rows, no_cols, ~] = size(mother_rgb);

        [~, image_name, ~] = fileparts(input_raw);
        orig_ppm_fname = fullfile(orig_residue_dir, "mother_image.ppm");
        imwrite(mother_rgb, orig_ppm_fname);

        % compress mother image with jpeg-1
        jpeg_fname = fullfile(jpgs_folder, "mother_image_001.jpg");
        if(codec == "UNSW_SJU_JPEG1huff")
            enc_command = sprintf("%s -quality %d -outfile %s %s", encoder_exe, q_v, jpeg_fname, orig_ppm_fname); % -arithmetic
        elseif(codec == "UNSW_SJU_JPEG1arith")
            enc_command = sprintf("%s -quality %d -arithmetic -outfile %s %s", encoder_exe, q_v, jpeg_fname, orig_ppm_fname); % -arithmetic
        elseif(codec == "UNSW_SJU_JPEG1huff_rdopt")
            if(exist(img_qfile_name, "file"))
                delete(img_qfile_name);
            end
            rdopt_cmd = sprintf("%s -numtables 2 -bpp %f -qfile %s -im %s", rdopt_exe, q_v, img_qfile_name, input_raw);
            [status, ~] = system(rdopt_cmd);
            enc_command =sprintf('%s -qtables %s -outfile %s %s', encoder_exe,  img_qfile_name, jpeg_fname, input_raw);
        elseif(codec == "UNSW_SJU_JPEG1arith_rdopt" || codec == "UNSW_SJU_JPEG1arith_rdopt_22p")
            if(exist(img_qfile_name, "file"))
                delete(img_qfile_name);
            end
            rdopt_cmd = sprintf("%s -numtables 2 -bpp %f -qfile %s -im %s", rdopt_exe, q_v, img_qfile_name, input_raw);
            [status, ~] = system(rdopt_cmd);
            enc_command =sprintf('%s -qtables %s -arithmetic -outfile %s %s', encoder_exe,  img_qfile_name, jpeg_fname, input_raw);
        end
        m_enc_stime = tic;
        [abc, out] = system(enc_command);
        if(abc)
            fprintf("Failed for [%s] to compress [%s] with [JPEG 1] codec.\n", codec, input_raw);
            error('compression failed');
        end
        img_enc_time = toc(m_enc_stime);
        % get bpp and encoding time
        jpg_image_file_info=dir(jpeg_fname);
        compressed_img_size = jpg_image_file_info.bytes;
        img_bpp = (compressed_img_size*8) / (no_cols*no_rows);

        decoded_ppm_fname = fullfile(set_dec_folder, "dec_mother_image.ppm");
        dec_command = sprintf("%s -outfile %s %s", decoder_exe, decoded_ppm_fname, jpeg_fname);
        m_dec_stime = tic;
        [xyz, out2] = system(dec_command);
        if(xyz)
            fprintf("Failed for [%s] to decompress [%s] with [JPEG 1] codec.\n", codec, jpg_image_file);
            error('decompression failed');
        end
        img_dec_time = toc(m_dec_stime);

        Recon_IMG = imread(decoded_ppm_fname);
        image_IQA_values = evaluate_images_IQA(mother_rgb, Recon_IMG, QA_metrics);
        ImagesResult_T1(1,:) ={set_name, 1, image_name, no_cols, no_rows, no_rows*no_cols, codec, q_v, compressed_img_size, img_bpp, img_enc_time, img_dec_time};
        Images_IQA_T = [Images_IQA_T;struct2table(image_IQA_values)];

        previous_recon_rgb = Recon_IMG;
        previous_orig_rgb = mother_rgb;
        q_start_time = tic;
        textprogressbar('        => Progress  : ');
        for imgNo =2:no_of_images
            next_raw_file = fullfile(dir_struct.set_ppm_path, ppm_images_names(imgNo).name);
            next_image_rgb = imread(next_raw_file);
            [h, w, c] = size(next_image_rgb);
            [no_rows, no_cols, ~] = size(next_image_rgb);
            [~, image_name, ~] = fileparts(next_raw_file);
            residue_gen_Tstart = tic;
            if(codec == "UNSW_SJU_JPEG1arith_rdopt_22p")
                residue_16bit = get_residue_regis_prediction_v3_based(previous_orig_rgb, previous_orig_rgb, next_image_rgb, ModelID, MAX_NUM_MODELS, set_name, set_compressed_folder, imgNo);
            else
                residue_16bit = get_residue_regis_prediction_v3_based(previous_orig_rgb, previous_recon_rgb, next_image_rgb, ModelID, MAX_NUM_MODELS, set_name, set_compressed_folder, imgNo);
            end
            residue_16bit_ppm_fname = fullfile(orig_residue_dir, sprintf("residue16b_image_%03d.ppm", imgNo));
            imwrite(residue_16bit, residue_16bit_ppm_fname, 'MaxValue', 65535);

            % convert 16bit residue to 8 bit residue using SJU indexing method
            [residue8bit, idx, idx_bytes] = convert_16bit_residue_to_8bit(residue_16bit, set_compressed_folder, imgNo);
            residue_gen_time =  toc(residue_gen_Tstart);

            % compress 8 bit residue with JPEG 1 in arithmetic mode
            residue_ppm_fname = fullfile(orig_residue_dir, sprintf("residue8b_image_%03d.ppm", imgNo));
            imwrite(residue8bit, residue_ppm_fname);
            residue_jpeg_fname = fullfile(jpgs_folder, sprintf("residue8b_image_%03d.jpg", imgNo));

            if(codec == "UNSW_SJU_JPEG1huff")
                enc_command = sprintf("%s -rgb -quality %d -outfile %s %s", encoder_exe, q_v, residue_jpeg_fname, residue_ppm_fname); % -arithmetic
            elseif(codec == "UNSW_SJU_JPEG1arith")
                enc_command = sprintf("%s -quality %d -rgb -arithmetic -outfile %s %s", encoder_exe, q_v, residue_jpeg_fname, residue_ppm_fname); % -arithmetic
            elseif(codec == "UNSW_SJU_JPEG1huff_rdopt" || codec == "UNSW_SJU_JPEG1arith_rdopt_22p")
                if(exist(img_qfile_name, "file"))
                    delete(img_qfile_name);
                end
                rdopt_cmd = sprintf("%s -numtables 2 -bpp %f -qfile %s -im %s", rdopt_exe, q_v, img_qfile_name, residue_ppm_fname);
                [status, ~] = system(rdopt_cmd);
                enc_command =sprintf('%s -qtables %s -rgb -outfile %s %s', encoder_exe,  img_qfile_name, residue_jpeg_fname, residue_ppm_fname);
            elseif(codec == "UNSW_SJU_JPEG1arith_rdopt")
                if(exist(img_qfile_name, "file"))
                    delete(img_qfile_name);
                end
                rdopt_cmd = sprintf("%s -numtables 2 -bpp %f -qfile %s -im %s", rdopt_exe, q_v, img_qfile_name, residue_ppm_fname);
                [status, ~] = system(rdopt_cmd);
                enc_command =sprintf('%s -qtables %s -rgb -arithmetic -outfile %s %s', encoder_exe,  img_qfile_name, residue_jpeg_fname, residue_ppm_fname);
            end


            enc_startT = tic;
            [abc, out] = system(enc_command);
            if(abc)
                fprintf("Failed for [%s] to compress [%s] with [JPEG 1] codec.\n", codec, residue_ppm_fname);
                error('compression failed');
            end
            img_enc_time = toc(enc_startT) + residue_gen_time;

            jpg_info1= dir(residue_jpeg_fname);
            compressed_img_size_a = jpg_info1.bytes;
            jpg_bpp_a = (compressed_img_size_a*8)/(w*h);
            compressed_img_size_b = compressed_img_size_a + idx_bytes;
            compressed_img_size = compressed_img_size_b + 592 +768;
            jpg_bpp_b = (compressed_img_size_b*8)/(w*h);
            img_bpp = (compressed_img_size *8)/(w*h);

            decoded_residue8bit_ppm_fname = fullfile(dec_residue_dir, sprintf("dec_residue8b_image_%03d.ppm", imgNo));
            dec_command = sprintf("%s -outfile %s %s", decoder_exe, decoded_residue8bit_ppm_fname, residue_jpeg_fname);
            dec_startT = tic;
            [xyz, out2] = system(dec_command);
            if(xyz)
                fprintf("Failed for [%s] to decompress [%s] with [JPEG 1] codec.\n", codec, residue_jpeg_fname);
                error('decompression failed');
            end
            img_dec_time = toc(dec_startT);

            % read decoded residue and evaluate
            decoded_residue8bit = imread(decoded_residue8bit_ppm_fname);
            psnr_a = psnr(decoded_residue8bit, residue8bit);
            mse_a = immse(decoded_residue8bit, residue8bit);

            % inverse conversion of decoded 8bit residue to 16 bit residue
            % using SJU indexing method
            recon_residue = convert_8bit_residue_to_16bit(decoded_residue8bit, idx);
            psnr_b = psnr(residue_16bit, recon_residue);
            mse_b = immse(residue_16bit, recon_residue);
            decoded_residue16bit_ppm_fname = fullfile(dec_residue_dir, sprintf("dec_residue16b_image_%03d.ppm", imgNo));
            imwrite(recon_residue, decoded_residue16bit_ppm_fname);

            % inverse residue
            residue_inv_tStart = tic;
            %get reconstructed image from 16 bit residue using USNW

            if(codec == "UNSW_SJU_JPEG1arith_rdopt_22p")
                recon_next_rgb = inverse_residue_regis_prediction_v3_based(set_name, set_compressed_folder, previous_orig_rgb, next_image_rgb, decoded_residue16bit_ppm_fname, ModelID, MAX_NUM_MODELS, QTree_Depth, SMALLEST_BLK_SIZE, imgNo);
            else
                recon_next_rgb = inverse_residue_regis_prediction_v3_based(set_name, set_compressed_folder, previous_recon_rgb, next_image_rgb, decoded_residue16bit_ppm_fname, ModelID, MAX_NUM_MODELS, QTree_Depth, SMALLEST_BLK_SIZE, imgNo);
            end

            residue_inv_time = toc(residue_inv_tStart);
            img_dec_time = img_dec_time+residue_inv_time;

            % Evaluation the reconstructed image
            previous_recon_rgb = recon_next_rgb;
            previous_orig_rgb = next_image_rgb;
            decoded_rgb_ppf_fname = fullfile(set_dec_folder, sprintf("dec_image_%03d.ppm", imgNo));
            imwrite(recon_next_rgb, decoded_rgb_ppf_fname);

            image_IQA_values = evaluate_images_IQA(next_image_rgb, recon_next_rgb, QA_metrics);
            ImagesResult_T1(imgNo,:) ={set_name, imgNo, image_name, no_cols, no_rows, no_rows*no_cols, codec, q_v, compressed_img_size, img_bpp, img_enc_time, img_dec_time};
            Images_IQA_T = [Images_IQA_T;struct2table(image_IQA_values)];
            textprogressbar(imgNo/no_of_images*100);
        end
        q_time_taken = toc(q_start_time);
        textprogressbar(['   done in  ' num2str(q_time_taken) ' sec']);



    elseif(codec == "UNSW_SJU_JPEG1arith_rdopt_v4" || codec == "UNSW_SJU_JPEG1arith_v4")
        %% UNSW_SJU_Arch
        ModelID = 8;
        MAX_NUM_MODELS = 3;
        QTree_Depth = 4;
        SMALLEST_BLK_SIZE = 16;
        encoder_exe = fullfile(dir_struct.codec_folder, "cjpeg.exe");
        decoder_exe = fullfile(dir_struct.codec_folder, "djpeg.exe");
        rdopt_exe = fullfile(dir_struct.codec_folder, "RD_Opt.exe");
        img_qfile_name = "qtable_tmp.txt";
        orig_residue_dir = fullfile(set_q_folder, "Orig_Residue");
        if(~exist(orig_residue_dir, 'dir')), mkdir(orig_residue_dir); end
        dec_residue_dir = fullfile(set_q_folder, "Dec_Residue");
        if(~exist(dec_residue_dir, 'dir')), mkdir(dec_residue_dir); end
        jpgs_folder = fullfile(set_compressed_folder, "Compressed_jpgs");
        if(~exist(jpgs_folder, 'dir')), mkdir(jpgs_folder); end
        input_raw = fullfile(dir_struct.set_ppm_path, ppm_images_names(1).name);
        mother_rgb = imread(input_raw);
        [no_rows, no_cols, ~] = size(mother_rgb);

        [~, image_name, ~] = fileparts(input_raw);
        orig_ppm_fname = fullfile(orig_residue_dir, "mother_image.ppm");
        imwrite(mother_rgb, orig_ppm_fname);

        % compress mother image with jpeg-1
        jpeg_fname = fullfile(jpgs_folder, "mother_image_001.jpg");

        if(codec == "UNSW_SJU_JPEG1arith_v4")
            enc_command = sprintf("%s -quality 90 -arithmetic -outfile %s %s", encoder_exe, jpeg_fname, orig_ppm_fname); % -arithmetic
        elseif(codec == "UNSW_SJU_JPEG1arith_rdopt_v4")
            if(exist(img_qfile_name, "file"))
                delete(img_qfile_name);
            end
            rdopt_cmd = sprintf("%s -numtables 2 -bpp %f -qfile %s -im %s", rdopt_exe, q_v, img_qfile_name, input_raw);
            [status, ~] = system(rdopt_cmd);
            enc_command =sprintf('%s -qtables %s -outfile %s %s', encoder_exe,  img_qfile_name, jpeg_fname, input_raw);
        end
        m_enc_stime = tic;
        [abc, out] = system(enc_command);
        if(abc)
            fprintf("Failed for [%s] to compress [%s] with [JPEG 1] codec.\n", codec, input_raw);
            error('compression failed');
        end
        img_enc_time = toc(m_enc_stime);
        % get bpp and encoding time
        jpg_image_file_info=dir(jpeg_fname);
        compressed_img_size = jpg_image_file_info.bytes;
        img_bpp = (compressed_img_size*8) / (no_cols*no_rows);

        decoded_ppm_fname = fullfile(set_dec_folder, "dec_mother_image.ppm");
        dec_command = sprintf("%s -outfile %s %s", decoder_exe, decoded_ppm_fname, jpeg_fname);
        m_dec_stime = tic;
        [xyz, out2] = system(dec_command);
        if(xyz)
            fprintf("Failed for [%s] to decompress [%s] with [JPEG 1] codec.\n", codec, jpg_image_file);
            error('decompression failed');
        end
        img_dec_time = toc(m_dec_stime);

        Recon_IMG = imread(decoded_ppm_fname);
        image_IQA_values = evaluate_images_IQA(mother_rgb, Recon_IMG, QA_metrics);
        ImagesResult_T1(1,:) ={set_name, 1, image_name, no_cols, no_rows, no_rows*no_cols, codec, q_v, compressed_img_size, img_bpp, img_enc_time, img_dec_time};
        Images_IQA_T = [Images_IQA_T;struct2table(image_IQA_values)];

        previous_recon_rgb = Recon_IMG;
        previous_orig_rgb = mother_rgb;
        q_start_time = tic;
        textprogressbar('        => Progress  : ');
        for imgNo =2:no_of_images
            next_raw_file = fullfile(dir_struct.set_ppm_path, ppm_images_names(imgNo).name);
            next_image_rgb = imread(next_raw_file);
            [h, w, c] = size(next_image_rgb);
            [no_rows, no_cols, ~] = size(next_image_rgb);
            [~, image_name, ~] = fileparts(next_raw_file);
            residue_gen_Tstart = tic;

            %             residue = get_residue_unsw_v4_based(previous_orig_rgb, ref , tgt, ModelID, MAX_NUM_MODELS, set_name, residue_dir, imgNo)
            residue_16bit = get_residue_unsw_v4_based(previous_orig_rgb, previous_recon_rgb, next_image_rgb, ModelID, MAX_NUM_MODELS, set_name, set_compressed_folder, imgNo);
            residue_16bit_ppm_fname = fullfile(orig_residue_dir, sprintf("residue16b_image_%03d.ppm", imgNo));
            imwrite(residue_16bit, residue_16bit_ppm_fname, 'MaxValue', 65535);

            % convert 16bit residue to 8 bit residue using SJU indexing method
            [residue8bit, idx, idx_bytes] = convert_16bit_residue_2_8bit_v4(residue_16bit, set_compressed_folder, imgNo);
            residue_gen_time =  toc(residue_gen_Tstart);

            % compress 8 bit residue with JPEG 1 in arithmetic mode
            residue_ppm_fname = fullfile(orig_residue_dir, sprintf("residue8b_image_%03d.ppm", imgNo));
            imwrite(residue8bit, residue_ppm_fname);
            residue_jpeg_fname = fullfile(jpgs_folder, sprintf("residue8b_image_%03d.jpg", imgNo));

            if(codec == "UNSW_SJU_JPEG1arith_v4")
                enc_command = sprintf("%s -quality %d -qtables ./codecs/qfixedtable10.txt -arithmetic -outfile %s %s", encoder_exe, q_v, residue_jpeg_fname, residue_ppm_fname); % -arithmetic
            elseif(codec == "UNSW_SJU_JPEG1arith_rdopt_v4")
                if(exist(img_qfile_name, "file"))
                    delete(img_qfile_name);
                end
                rdopt_cmd = sprintf("%s -numtables 2 -bpp %f -qfile %s -im %s", rdopt_exe, q_v, img_qfile_name, residue_ppm_fname);
                [status, ~] = system(rdopt_cmd);
                enc_command =sprintf('%s -qtables %s -outfile %s %s', encoder_exe,  img_qfile_name, residue_jpeg_fname, residue_ppm_fname);
            end


            enc_startT = tic;
            [abc, out] = system(enc_command);
            if(abc)
                fprintf("Failed for [%s] to compress [%s] with [JPEG 1] codec.\n", codec, residue_ppm_fname);
                error('compression failed');
            end
            img_enc_time = toc(enc_startT) + residue_gen_time;

            jpg_info1= dir(residue_jpeg_fname);
            compressed_img_size_a = jpg_info1.bytes;
            jpg_bpp_a = (compressed_img_size_a*8)/(w*h);
            compressed_img_size_b = compressed_img_size_a + idx_bytes;
            compressed_img_size = compressed_img_size_b + 592 +768;
            jpg_bpp_b = (compressed_img_size_b*8)/(w*h);
            img_bpp = (compressed_img_size *8)/(w*h);

            decoded_residue8bit_ppm_fname = fullfile(dec_residue_dir, sprintf("dec_residue8b_image_%03d.ppm", imgNo));
            dec_command = sprintf("%s -outfile %s %s", decoder_exe, decoded_residue8bit_ppm_fname, residue_jpeg_fname);
            dec_startT = tic;
            [xyz, out2] = system(dec_command);
            if(xyz)
                fprintf("Failed for [%s] to decompress [%s] with [JPEG 1] codec.\n", codec, residue_jpeg_fname);
                error('decompression failed');
            end
            img_dec_time = toc(dec_startT);

            % read decoded residue and evaluate
            decoded_residue8bit = imread(decoded_residue8bit_ppm_fname);
            psnr_a = psnr(decoded_residue8bit, residue8bit);
            mse_a = immse(decoded_residue8bit, residue8bit);

            % inverse conversion of decoded 8bit residue to 16 bit residue
            % using SJU indexing method
            recon_residue = convert_8bit_residue_to_16bit_v4(decoded_residue8bit, idx);
            psnr_b = psnr(residue_16bit, recon_residue);
            mse_b = immse(residue_16bit, recon_residue);
            decoded_residue16bit_ppm_fname = fullfile(dec_residue_dir, sprintf("dec_residue16b_image_%03d.ppm", imgNo));
            imwrite(recon_residue, decoded_residue16bit_ppm_fname);

            % inverse residue
            residue_inv_tStart = tic;
            %get reconstructed image from 16 bit residue using USNW

            recon_next_rgb = inverse_residue_unsw_v4_based(set_name, set_compressed_folder, previous_recon_rgb, next_image_rgb, decoded_residue16bit_ppm_fname, ModelID, MAX_NUM_MODELS, QTree_Depth, SMALLEST_BLK_SIZE, imgNo);

            residue_inv_time = toc(residue_inv_tStart);
            img_dec_time = img_dec_time+residue_inv_time;

            % Evaluation the reconstructed image
            previous_recon_rgb = recon_next_rgb;
            previous_orig_rgb = next_image_rgb;
            decoded_rgb_ppf_fname = fullfile(set_dec_folder, sprintf("dec_image_%03d.ppm", imgNo));
            imwrite(recon_next_rgb, decoded_rgb_ppf_fname);

            image_IQA_values = evaluate_images_IQA(next_image_rgb, recon_next_rgb, QA_metrics);
            ImagesResult_T1(imgNo,:) ={set_name, imgNo, image_name, no_cols, no_rows, no_rows*no_cols, codec, q_v, compressed_img_size, img_bpp, img_enc_time, img_dec_time};
            Images_IQA_T = [Images_IQA_T;struct2table(image_IQA_values)];
            textprogressbar(imgNo/no_of_images*100);
        end
        q_time_taken = toc(q_start_time);
        textprogressbar(['   done in  ' num2str(q_time_taken) ' sec']);


    elseif strcmp(codec, "VVC_VTM_Inter") || strcmp(codec, "VVC_VTM_Intra")
        %% VVC-VTM Inter
        input_yuv_video = exp_params.yuv_video_file;
        input_raw = fullfile(orig_files(1).folder, orig_files(1).name);

        Orig_IMG = imread(input_raw);
        [no_rows, no_cols, ~] = size(Orig_IMG);
        [~, image_name, ~] = fileparts(input_raw);

        encoder_exe = fullfile(dir_struct.codec_folder, "VVC_VTM", "EncoderApp.exe");
        bin_video_file = fullfile(set_compressed_folder, sprintf("%s_%s.bin", set_name, q_fieldName));
        decoded_video_yuv = fullfile(set_dec_folder, sprintf("Rec_%s_%s.yuv", set_name, q_fieldName));

        if strcmp(codec, "VVC_VTM_Inter")
            cfg_file = fullfile(dir_struct.codec_folder, "VVC_VTM","cfg", "encoder_randomaccess_vtm.cfg");
            command =sprintf('%s -c %s -i %s -wdt %d -hgt %d -b %s -f %d -fr 30 -q %d --InputBitDepth=10 --InputChromaFormat=420 --ChromaFormatIDC=420 --TemporalSubsampleRatio=1 --ConformanceWindowMode=1 --MinSearchWindow 1 -o %s', encoder_exe, cfg_file, input_yuv_video, no_cols, no_rows, bin_video_file, no_of_images, q_v, decoded_video_yuv);
        else
            cfg_file = fullfile(dir_struct.codec_folder, "VVC_VTM","cfg", "encoder_intra_vtm.cfg");
            command =sprintf('%s -c %s -i %s -wdt %d -hgt %d -b %s -f %d -fr 30 -q %d --InputBitDepth=10 --InputChromaFormat=420 --ChromaFormatIDC=420 --TemporalSubsampleRatio=1 --ConformanceWindowMode=1 --IntraPeriod 1 -o %s', encoder_exe, cfg_file, input_yuv_video, no_cols, no_rows, bin_video_file, no_of_images, q_v, decoded_video_yuv);        
        end
        tCstartimg= tic;
        [status, outE] = system(command);
        img_enc_time  = toc(tCstartimg);
        if(status)
            fprintf("VVC codec=> SetName: %s, QP: %d, => Failed to compress [%s].\n", set_name,  q_v, input_yuv_video);
            error('compression failed');
        end

        fprintf("\t\tCompression/Decompression Done => VVC codec=> SetName: %s, QP: %d, =>  [%s].\n", set_name,  q_v, input_yuv_video);
        bin_video_file_info=dir(bin_video_file);
        compressed_video_size = bin_video_file_info.bytes;
        set_bpp = (compressed_video_size*8) / (no_cols*no_rows*no_of_images);

        num_of_bits_Table = get_images_bpp_vvc_out(outE);
        time_taken = extract_time(outE);


        decoded_rgb_files = fullfile(fullfile(set_dec_folder, sprintf("Recon_%s_ImgNO%%03d_%s.ppm", set_name, q_fieldName)));
        conv_cmd = sprintf('ffmpeg.exe -f rawvideo -vcodec rawvideo -s %dx%d -pix_fmt yuv420p10le -i %s -pix_fmt rgb24 -vf scale=in_range=full:in_color_matrix=bt709:out_range=full:out_color_matrix=bt709 -color_primaries bt709 -color_trc bt709 -colorspace bt709 -y %s', no_cols, no_rows, decoded_video_yuv, decoded_rgb_files);
        [status, ~] = system(conv_cmd);
        if(status)
            fprintf("Failed to convert [%s] with [ffmpeg] to [%s].\n", decoded_video_yuv, decoded_rgb_files);
            error('conversion failed');
        end
        tx = toc(tCstartimg);
        fprintf("\t\tCompression/Decompression and Conversion Done in %f sec.\n", tx);
        delete(decoded_video_yuv);
        decoded_images_list = dir(fullfile(set_dec_folder, '*.ppm'));

        textprogressbar('        => Evaluation Progress  : ');
        for imgNo = 1: no_of_images
            input_raw = fullfile(orig_files(imgNo).folder, orig_files(imgNo).name);

            Orig_IMG = imread(input_raw);
            [no_rows, no_cols, ~] = size(Orig_IMG);
            [~, image_name, ~] = fileparts(input_raw);

            decoded_rgb_file = fullfile(decoded_images_list(imgNo).folder, decoded_images_list(imgNo).name);

            Recon_IMG = imread(decoded_rgb_file);
            image_IQA_values = evaluate_images_IQA(Orig_IMG, Recon_IMG, QA_metrics);
            row = num_of_bits_Table.Image_No == imgNo;
            num_of_bits = num_of_bits_Table.Num_of_Bits(row);
            num_of_pixels = no_rows*no_cols;
            compressed_img_size = num_of_bits/8; % to convert to bytes from bits
            img_bpp = num_of_bits/num_of_pixels;
            img_enc_time = time_taken/2/no_of_images;
            img_dec_time = time_taken/2/no_of_images;
            ImagesResult_T1(imgNo,:) ={set_name, imgNo, image_name, no_cols, no_rows, num_of_pixels, codec, q_v, compressed_img_size, img_bpp, img_enc_time, img_dec_time};
            Images_IQA_T = [Images_IQA_T;struct2table(image_IQA_values)];
            textprogressbar(imgNo/no_of_images*100);
        end
        q_time_taken = toc(q_start_time);
        textprogressbar(['   done in  ' num2str(q_time_taken) ' sec']);
    elseif strcmp(codec, "VVC_VVenC_Inter") || strcmp(codec, "VVC_VVenC_Intra")
        %% VVC_VVenC_Inter
        input_yuv_video = exp_params.yuv_video_file;
        input_raw = fullfile(orig_files(1).folder, orig_files(1).name);

        Orig_IMG = imread(input_raw);
        [no_rows, no_cols, ~] = size(Orig_IMG);
        [~, image_name, ~] = fileparts(input_raw);

        encoder_exe = fullfile(dir_struct.codec_folder, "VVC_VVenC", "vvencFFapp.exe");
        bin_video_file = fullfile(set_compressed_folder, sprintf("%s_%s.266", set_name, q_fieldName));
        rnd_cfg_file = fullfile(dir_struct.codec_folder, "VVC_VVenC/cfg/randomaccess_slow.cfg");
        seq_cfg_file = fullfile(dir_struct.codec_folder, "VVC_VVenC/cfg/sequence.cfg");

        %command =sprintf('%s -c %s -i %s -wdt %d -hgt %d -b %s -f %d -fr 30 -q %d --InputBitDepth=10 --InputChromaFormat=444 --ChromaFormatIDC=444 --TemporalSubsampleRatio=1 --ConformanceWindowMode=1 -o %s', encoder_exe, cfg_file, input_yuv_video, no_cols, no_rows, bin_video_file, no_of_images, q_v, decoded_video_yuv);
        if(strcmp(codec, "VVC_VVenC_Inter"))
          enc_cmd = sprintf("%s --preset=slow -c %s -c %s -q %d --InputFile %s -s %dx%d --InputBitDepth 10 --FramesToBeEncoded %d --BitstreamFile %s --InputChromaFormat 420 --MinSearchWindow 1 --framerate 30", encoder_exe, rnd_cfg_file, seq_cfg_file, q_v, input_yuv_video, no_cols, no_rows, no_of_images, bin_video_file);
        else  
            enc_cmd = sprintf("%s --preset=slow -c %s -c %s -q %d --InputFile %s -s %dx%d --InputBitDepth 10 --FramesToBeEncoded %d --BitstreamFile %s --InputChromaFormat 420 --framerate 30 --IntraPeriod 1", encoder_exe, rnd_cfg_file, seq_cfg_file, q_v, input_yuv_video, no_cols, no_rows, no_of_images, bin_video_file);            
        end

        tCstartimg= tic;
        [status, outE] = system(enc_cmd);
        vid_enc_time  = toc(tCstartimg);
        if(status)
            fprintf("VVC VVC_VVenC => SetName: %s, QP: %d, => Failed to compress [%s].\n", set_name,  q_v, input_yuv_video);
            error('compression failed');
        end
        fprintf("\t\tCompression Done in %f => VVC codec=> SetName: %s, QP: %d, =>  [%s].\n", vid_enc_time, set_name,  q_v, input_yuv_video);
        bin_video_file_info=dir(bin_video_file);
        compressed_video_size = bin_video_file_info.bytes;
        set_bpp = (compressed_video_size*8) / (no_cols*no_rows*no_of_images);

        num_of_bits_Table = get_images_bpp_vvc_out(outE);
        time_taken = extract_time(outE);

        % decoding
        decoder_exe = fullfile(dir_struct.codec_folder, "VVC_VVenC/vvdecapp.exe");
        decoded_video_yuv = fullfile(set_dec_folder, sprintf("Rec_%s_%s.yuv", set_name, q_fieldName));

        dec_cmd = sprintf("%s -b %s -o %s ", decoder_exe,  bin_video_file, decoded_video_yuv);

        tDstartimg= tic;
        [status, outD] = system(dec_cmd);
        vid_dec_time  = toc(tDstartimg);
        if(status)
            fprintf("VVC vvdecapp => SetName: %s, QP: %d, => Failed to compress [%s].\n", set_name,  q_v, input_yuv_video);
            error('compression failed');
        end
        fprintf("\t\tDecompression Done in %f sec => VVC codec=> SetName: %s, QP: %d, =>  [%s].\n", vid_dec_time, set_name,  q_v, decoded_video_yuv);

        tConStartT = tic;
        decoded_rgb_files = fullfile(fullfile(set_dec_folder, sprintf("Recon_%s_ImgNO%%03d_%s.ppm", set_name, q_fieldName)));
        conv_cmd = sprintf('ffmpeg.exe -f rawvideo -vcodec rawvideo -s %dx%d -pix_fmt yuv420p10le -i %s -pix_fmt rgb24 -vf scale=in_range=full:in_color_matrix=bt709:out_range=full:out_color_matrix=bt709 -color_primaries bt709 -color_trc bt709 -colorspace bt709 -y %s', no_cols, no_rows, decoded_video_yuv, decoded_rgb_files);
        [status, ~] = system(conv_cmd);
        if(status)
            fprintf("Failed to convert [%s] with [ffmpeg] to [%s].\n", decoded_video_yuv, decoded_rgb_files);
            error('conversion failed');
        end
        tx = toc(tConStartT);
        fprintf("\t\tConversion Done in %f sec.\n", tx);
        delete(decoded_video_yuv);
        decoded_images_list = dir(fullfile(set_dec_folder, '*.ppm'));

        textprogressbar('        => Evaluation Progress  : ');
        for imgNo = 1: no_of_images
            input_raw = fullfile(orig_files(imgNo).folder, orig_files(imgNo).name);

            Orig_IMG = imread(input_raw);
            [no_rows, no_cols, ~] = size(Orig_IMG);
            [~, image_name, ~] = fileparts(input_raw);

            decoded_rgb_file = fullfile(decoded_images_list(imgNo).folder, decoded_images_list(imgNo).name);

            Recon_IMG = imread(decoded_rgb_file);
            image_IQA_values = evaluate_images_IQA(Orig_IMG, Recon_IMG, QA_metrics);
            row = num_of_bits_Table.Image_No == imgNo;
            num_of_bits = num_of_bits_Table.Num_of_Bits(row);
            num_of_pixels = no_rows*no_cols;
            compressed_img_size = num_of_bits/8; % to convert to bytes from bits
            img_bpp = num_of_bits/num_of_pixels;
            img_enc_time = vid_enc_time/no_of_images;
            img_dec_time = vid_dec_time/no_of_images;
            ImagesResult_T1(imgNo,:) ={set_name, imgNo, image_name, no_cols, no_rows, num_of_pixels, codec, q_v, compressed_img_size, img_bpp, img_enc_time, img_dec_time};
            Images_IQA_T = [Images_IQA_T;struct2table(image_IQA_values)];
            textprogressbar(imgNo/no_of_images*100);
        end
        q_time_taken = toc(q_start_time);
        textprogressbar(['   done in  ' num2str(q_time_taken) ' sec']);
   
    else
        %% JPEG1, JPEG 2000 and JPEG XL codecs
        textprogressbar('        => Progress  : ');
        for imgNo = 1: no_of_images
            ppm_images_names = dir(fullfile(dir_struct.set_ppm_path, "*.ppm"));
            input_raw = fullfile(dir_struct.set_ppm_path, ppm_images_names(imgNo).name);
            Orig_IMG = imread(input_raw);
            [no_rows, no_cols, ~] = size(Orig_IMG);
            [~, image_name, ~] = fileparts(input_raw);
            if(codec == "JPEG1" || codec == "JPEG1_Arithmetic")
                encoder_exe = fullfile(dir_struct.codec_folder, "JPEG1", "cjpeg.exe");
                decoder_exe = fullfile(dir_struct.codec_folder, "JPEG1", "djpeg.exe");
                jpg_image_file = fullfile(set_compressed_folder, sprintf("%s_ImgNO%03d_%s.jpg", image_name, imgNo, q_fieldName));
                decoded_image_ppm = fullfile(set_dec_folder, sprintf("Recon_%s_ImgNO%03d_%s.ppm", image_name, imgNo, q_fieldName));
                if(codec == "JPEG1")
                    command =sprintf('%s -quality %d -outfile %s %s', encoder_exe,  q_v, jpg_image_file, input_raw);
                else
                    command =sprintf('%s -quality %d -arithmetic -outfile %s %s', encoder_exe,  q_v, jpg_image_file, input_raw);
                end
                tCstartimg= tic;
                [status, ~] = system(command);
                img_enc_time  = toc(tCstartimg);
                if(status)
                    fprintf("Failed to compress [%s] with [JPEG 1] codec.\n", input_raw);
                    error('compression failed');
                end
                % get bpp and encoding time
                jpg_image_file_info=dir(jpg_image_file);
                compressed_img_size = jpg_image_file_info.bytes;
                img_bpp = (compressed_img_size*8) / (no_cols*no_rows);
                command=sprintf('%s -outfile %s %s',decoder_exe, decoded_image_ppm, jpg_image_file);
                tDstartimg = tic;
                [status, ~] = system(command);
                if(status)
                    fprintf("Failed to decompress [%s] with [JPEG 1] codec.\n", jpg_image_file);
                    error('decompression failed');
                end
                img_dec_time = toc(tDstartimg);
            elseif(codec == "JPEG1huff_rdopt" || codec == "JPEG1arith_rdopt" || codec == "SJU_Arch_JPEG1arith_rdopt1")
                encoder_exe = fullfile(dir_struct.codec_folder, "cjpeg.exe");
                decoder_exe = fullfile(dir_struct.codec_folder, "djpeg.exe");
                rdopt_exe = fullfile(dir_struct.codec_folder, "RD_Opt.exe");
                img_qfile_name = "qtable_tmp.txt";
                rdopt_cmd = sprintf("%s -numtables 2 -bpp %f -qfile %s -im %s", rdopt_exe, q_v, img_qfile_name, input_raw);
                [status, ~] = system(rdopt_cmd);

                jpg_image_file = fullfile(set_compressed_folder, sprintf("%s_ImgNO%03d_%s.jpg", image_name, imgNo, q_fieldName));
                decoded_image_ppm = fullfile(set_dec_folder, sprintf("Recon_%s_ImgNO%03d_%s.ppm", image_name, imgNo, q_fieldName));
                if(codec == "JPEG1huff_rdopt")
                    command =sprintf('%s -qtables %s -outfile %s %s', encoder_exe,  img_qfile_name, jpg_image_file, input_raw);
                elseif(codec == "JPEG1arith_rdopt")
                    command =sprintf('%s -qtables %s -arithmetic -outfile %s %s', encoder_exe,  img_qfile_name, jpg_image_file, input_raw);
                elseif(codec == "SJU_Arch_JPEG1arith_rdopt1")
                    command =sprintf('%s -qtables %s -arithmetic -dct float -outfile %s %s', encoder_exe,  img_qfile_name, jpg_image_file, input_raw);
                end
                tCstartimg= tic;
                [status, ~] = system(command);
                img_enc_time  = toc(tCstartimg);
                if(status)
                    fprintf("Failed to compress [%s] with [JPEG 1] codec.\n", input_raw);
                    error('compression failed');
                end
                % get bpp and encoding time
                jpg_image_file_info=dir(jpg_image_file);
                compressed_img_size = jpg_image_file_info.bytes;
                img_bpp = (compressed_img_size*8) / (no_cols*no_rows);
                command=sprintf('%s -outfile %s %s',decoder_exe, decoded_image_ppm, jpg_image_file);
                tDstartimg = tic;
                [status, ~] = system(command);
                if(status)
                    fprintf("Failed to decompress [%s] with [JPEG 1] codec.\n", jpg_image_file);
                    error('decompression failed');
                end
                img_dec_time = toc(tDstartimg);
                delete(img_qfile_name);
            elseif(codec == "JPEG2000")
                encoder_exe = fullfile(dir_struct.codec_folder, "JPEG_2000", "kdu_compress.exe");
                decoder_exe = fullfile(dir_struct.codec_folder, "JPEG_2000", "kdu_expand.exe");
                j2c_image_file = fullfile(set_compressed_folder, sprintf("%s_ImgNo%03d_%s.j2c", image_name, imgNo, q_fieldName));
                decoded_image_ppm = fullfile(set_dec_folder, sprintf("Recon_%s_ImgNo%03d_%s.ppm", image_name, imgNo, q_fieldName));
                command = sprintf("%s -i %s -o %s Qfactor=%d", encoder_exe, input_raw, j2c_image_file, q_v);
                tCstartimg= tic;
                [status, outTxt] = system(command);
                img_enc_time  = toc(tCstartimg);
                if(status)
                    fprintf("Failed to compress [%s] with [JPEG 2000] codec.\n", input_raw);
                    error('compression failed');
                end
                % get bpp and encoding time
                jpg_image_file_info=dir(j2c_image_file);
                compressed_img_size = jpg_image_file_info.bytes;
                img_bpp = (compressed_img_size*8) / (no_cols*no_rows);
                command = sprintf("%s -i %s -o %s", decoder_exe, j2c_image_file, decoded_image_ppm);
                tDstartimg = tic;
                [status, ~] = system(command);
                if(status)
                    fprintf("Failed to decompress [%s] with [JPEG 2000] codec.\n", j2c_image_file);
                    error('decompression failed');
                end
                img_dec_time = toc(tDstartimg);
            elseif(codec == "JPEGXL")
                encoder_exe = fullfile(dir_struct.codec_folder, "JPEG_XL", "cjxl.exe");
                decoder_exe = fullfile(dir_struct.codec_folder, "JPEG_XL", "djxl.exe");
                jxl_image_file = fullfile(set_compressed_folder, sprintf("%s_ImgNo%03d_%s.jxl", image_name, imgNo, q_fieldName));
                decoded_image_ppm = fullfile(set_dec_folder, sprintf("Recon_%s_ImgNo%03d_%s.ppm", image_name, imgNo, q_fieldName));
                command =sprintf('%s %s %s -q %d -e 6 --num_threads 0', encoder_exe, input_raw, jxl_image_file, q_v);

                tCstartimg= tic;
                [status, ~] = system(command);
                img_enc_time  = toc(tCstartimg);
                if(status)
                    fprintf("Failed to compress [%s] with [JPEG XL] codec.\n", input_raw);
                    error('compression failed');
                end
                % get bpp and encoding time
                jpg_image_file_info=dir(jxl_image_file);
                compressed_img_size = jpg_image_file_info.bytes;
                img_bpp = (compressed_img_size*8) / (no_cols*no_rows);
                command=sprintf('%s %s %s --num_threads 0',decoder_exe, jxl_image_file, decoded_image_ppm);
                tDstartimg = tic;
                [status, ~] = system(command);
                if(status)
                    fprintf("Failed to decompress [%s] with [JPEG XL] codec.\n", jxl_image_file);
                    error('decompression failed');
                end
                img_dec_time = toc(tDstartimg);

            else
                fprintf("Codec [%s] functionality not included yet", codec);
                error("Failed");
            end

            Recon_IMG = imread(decoded_image_ppm);
            image_IQA_values = evaluate_images_IQA(Orig_IMG, Recon_IMG, QA_metrics);
            ImagesResult_T1(imgNo,:) ={set_name, imgNo, image_name, no_cols, no_rows, no_rows*no_cols, codec, q_v, compressed_img_size, img_bpp, img_enc_time, img_dec_time};
            Images_IQA_T = [Images_IQA_T;struct2table(image_IQA_values)];
            textprogressbar(imgNo/no_of_images*100);
        end
        q_time_taken = toc(q_start_time);
        textprogressbar(['   done in  ' num2str(q_time_taken) ' sec']);
    end

    print_current_time();
    fprintf("\n");
    images_results_table = [ImagesResult_T1, Images_IQA_T];
    %     writetable(images_results_table, dir_struct.results_MS_excel_fname ,'Sheet', 'ImagesResults','WriteMode','Append')%, 'WriteVariableNames',false,'WriteRowNames',true)
    set_codec_results.(q_fieldName) = images_results_table;
    M.(set_name).(codec).q_fields.(q_fieldName) = images_results_table;
    save(fullfile(exp_params.dir_struct.results_folders, exp_params.results_mat_fname), 'M');

    % delete directorires and files
    if(do_discard)
        rmdir(set_compressed_folder, "s");
        rmdir(set_dec_folder, "s");
        rmdir(set_q_folder, "s");
    end
end
images_names = strings(1,no_of_images);
for i =1:no_of_images
    images_names(1, i) = orig_files(i).name;
end
M.(set_name).(codec).SET_AVG_RESULTS = calculate_set_results(M.(set_name).(codec).q_fields, QA_metrics);
% codec_results.SET_DETAILED_RESULTS = set_codec_results;
M.(set_name).(codec).IMAGES_RESULT = calculate_images_results(M.(set_name).(codec).q_fields,no_of_images, QA_metrics);
M.(set_name).(codec).RATE_VALUES = q_list;
M.(set_name).(codec).IMAGES_NAMES = images_names;

save(fullfile(exp_params.dir_struct.results_folders, exp_params.results_mat_fname), 'M');

codec_time_taken = toc(codec_exp_time);
fprintf("<strong>  => Experiment Done : => SET = [%s],  => CODEC = [%s]  => TIME TAKEN = [%f sec]</strong>\n\n", set_name, codec, codec_time_taken);
if(do_discard)
    rmdir(codec_data_folder, "s");
end
end

function print_current_time()
    nowTime = datetime('now', 'Format','dd-MMM-uuuu HH:mm:ss');
    timeStr = string(nowTime);
    fprintf("\tTime Now: %s", timeStr);
end

function codec_avg_results = calculate_set_results(set_codec_results, QA_metrics)
qfield_list = fieldnames(set_codec_results);
no_of_q_points = length(qfield_list);
avg_result_T1 = table('Size', [no_of_q_points 6], 'VariableTypes', {'uint16', 'single', 'uint32', 'double','double','double'});
avg_result_T1.Properties.VariableNames = {'QP_No','Qvalue', 'SET_CompressedSize','SET_BPP', 'SET_ENC_TIME', 'SET_DEC_TIME'};

for i = 1 : no_of_q_points
    q_p_table = set_codec_results.(qfield_list{i});
    q_p_value = q_p_table.rate(1);
    set_compressed_size = sum(q_p_table.CompressedImageSize);
    set_bpp = (set_compressed_size*8)/(sum(q_p_table.PixelsCount));
    set_enc_time = sum(q_p_table.enc_time);
    set_dec_time = sum(q_p_table.dec_time);
    avg_result_T1(i,:) = {i, q_p_value, set_compressed_size, set_bpp, set_enc_time, set_dec_time};

end

no_of_metrics = length(QA_metrics);
iqa_mat = zeros(no_of_q_points, no_of_metrics);
for i = 1: no_of_q_points
    q_p_table = set_codec_results.(qfield_list{i});
    for m = 1: no_of_metrics
        iqa_mat(i, m) = mean(q_p_table.(QA_metrics(m)));
    end
end
avg_IQA_table = array2table(iqa_mat);
avg_IQA_table.Properties.VariableNames = QA_metrics;
codec_avg_results = [avg_result_T1 avg_IQA_table];
end

function detailed_results = calculate_images_results(set_codec_results,no_of_images, QA_metrics)
detailed_results = struct();
qfield_list = fieldnames(set_codec_results);
no_of_q_points = length(qfield_list);
images_names = set_codec_results.(qfield_list{1}).ImageName;
extra_columns = 7; % S No, Rate, MEAN, MAX, MIN, STD, VAR
images_comp_sizes_mat = uint32(zeros(no_of_q_points, no_of_images + extra_columns));
images_bpp_mat = single(zeros(no_of_q_points, no_of_images + extra_columns));
images_enc_time_mat = single(zeros(no_of_q_points, no_of_images + extra_columns));
images_dec_time_mat = single(zeros(no_of_q_points, no_of_images + extra_columns));

for qf = 1: no_of_q_points
    images_bpp_mat(qf, 1) = qf;
    images_comp_sizes_mat(qf, 1) = qf;
    images_enc_time_mat(qf, 1) = qf;
    images_dec_time_mat(qf, 1) = qf;
    q_value = set_codec_results.(qfield_list{qf}).rate(1);
    images_bpp_mat(qf, 2) = q_value;
    images_comp_sizes_mat(qf, 2) = q_value;
    images_enc_time_mat(qf, 2) = q_value;
    images_dec_time_mat(qf, 2) = q_value;
    for img =1 : no_of_images
        images_comp_sizes_mat(qf, img+extra_columns) = set_codec_results.(qfield_list{qf}).CompressedImageSize(img);
        images_bpp_mat(qf, img+extra_columns) = set_codec_results.(qfield_list{qf}).bpp(img);
        images_enc_time_mat(qf, img+extra_columns) = set_codec_results.(qfield_list{qf}).enc_time(img);
        images_dec_time_mat(qf, img+extra_columns) = set_codec_results.(qfield_list{qf}).dec_time(img);
    end
    images_comp_sizes_mat(qf, 3) = mean(images_comp_sizes_mat(qf, extra_columns+1:end));
    images_comp_sizes_mat(qf, 4) = max(images_comp_sizes_mat(qf, extra_columns+1:end));
    images_comp_sizes_mat(qf, 5) = min(images_comp_sizes_mat(qf, extra_columns+1:end));
    images_comp_sizes_mat(qf, 6) = std(double(images_comp_sizes_mat(qf, extra_columns+1:end)));
    images_comp_sizes_mat(qf, 7) = var(double(images_comp_sizes_mat(qf, extra_columns+1:end)));
    images_bpp_mat(qf, 3) = mean(images_bpp_mat(qf, extra_columns+1:end));
    images_bpp_mat(qf, 4) = max(images_bpp_mat(qf, extra_columns+1:end));
    images_bpp_mat(qf, 5) = min(images_bpp_mat(qf, extra_columns+1:end));
    images_bpp_mat(qf, 6) = std(images_bpp_mat(qf, extra_columns+1:end));
    images_bpp_mat(qf, 7) = var(images_bpp_mat(qf, extra_columns+1:end));
    images_enc_time_mat(qf, 3) = mean(images_enc_time_mat(qf, extra_columns+1:end));
    images_enc_time_mat(qf, 4) = max(images_enc_time_mat(qf, extra_columns+1:end));
    images_enc_time_mat(qf, 5) = min(images_enc_time_mat(qf, extra_columns+1:end));
    images_enc_time_mat(qf, 6) = std(images_enc_time_mat(qf, extra_columns+1:end));
    images_enc_time_mat(qf, 7) = var(images_enc_time_mat(qf, extra_columns+1:end));
    images_dec_time_mat(qf, 3) = mean(images_dec_time_mat(qf, extra_columns+1:end));
    images_dec_time_mat(qf, 4) = max(images_dec_time_mat(qf, extra_columns+1:end));
    images_dec_time_mat(qf, 5) = min(images_dec_time_mat(qf, extra_columns+1:end));
    images_dec_time_mat(qf, 6) = std(images_dec_time_mat(qf, extra_columns+1:end));
    images_dec_time_mat(qf, 7) = var(images_dec_time_mat(qf, extra_columns+1:end));
end

table_header = {'S No', 'rate', 'MEAN', 'MAX', 'MIN', 'STD', 'VAR'};
Compressed_SizesT = array2table(images_comp_sizes_mat);
Compressed_SizesT.Properties.VariableNames(1:extra_columns) = table_header;
Compressed_SizesT.Properties.VariableNames(extra_columns+1:end) = images_names;
bpp_T = array2table(images_bpp_mat);
bpp_T.Properties.VariableNames(1:extra_columns) = table_header;
bpp_T.Properties.VariableNames(extra_columns+1:end) = images_names;
enc_time_T = array2table(images_enc_time_mat);
enc_time_T.Properties.VariableNames(1:extra_columns) = table_header;
enc_time_T.Properties.VariableNames(extra_columns+1:end) = images_names;
dec_time_T = array2table(images_dec_time_mat);
dec_time_T.Properties.VariableNames(1:extra_columns) = table_header;
dec_time_T.Properties.VariableNames(extra_columns+1:end) = images_names;
detailed_results.CompressedSizes = Compressed_SizesT;
detailed_results.BPPs = bpp_T;
detailed_results.ENC_TIMES = enc_time_T;
detailed_results.DEC_TIMES = dec_time_T;

metric_results = struct();
no_of_metrics = length(QA_metrics);
for m = 1: no_of_metrics
    metric_v = QA_metrics(m);
    metric_mat = single(zeros(no_of_q_points, no_of_images + extra_columns));
    for qf = 1: no_of_q_points
        metric_mat(qf, 1) = qf;
        q_value = set_codec_results.(qfield_list{qf}).rate(1);
        metric_mat(qf, 2) = q_value;
        for img = 1: no_of_images
            metric_mat(qf, img+extra_columns) = set_codec_results.(qfield_list{qf}).(metric_v)(img);
        end
        metric_mat(qf, 3) = mean(metric_mat(qf, extra_columns+1:end));
        metric_mat(qf, 4) = max(metric_mat(qf, extra_columns+1:end));
        metric_mat(qf, 5) = min(metric_mat(qf, extra_columns+1:end));
        metric_mat(qf, 6) = std(metric_mat(qf, extra_columns+1:end));
        metric_mat(qf, 7) = var(metric_mat(qf, extra_columns+1:end));
    end
    metric_T = array2table(metric_mat);
    metric_T.Properties.VariableNames(1:extra_columns) = table_header;
    metric_T.Properties.VariableNames(extra_columns+1:end) = images_names;
    metric_results.(metric_v) = metric_T;
end

detailed_results.METRICS = metric_results;

end

function exp_params = convert_dataset(exp_params, codec_in_file_format)
orig_ext = exp_params.images_files_ext;
original_images_list = dir(fullfile(exp_params.set_path, orig_ext));
no_of_images = numel(original_images_list);
exp_params.no_of_images = no_of_images;
orig_ext = extractAfter(orig_ext,1);
dir_struct = exp_params.dir_struct;
tCstart = tic;
if strcmp(codec_in_file_format, "PPM")
    ppm_folder = fullfile(dir_struct.ground_truth, codec_in_file_format);
    dir_struct.set_ppm_path = ppm_folder;
    if(~exist(ppm_folder, 'dir')), mkdir(ppm_folder); end
    textprogressbar('=> Conversion to PPM : ');
    for img=1:no_of_images
        rgb_image = original_images_list(img).name;
        ppm_name = strrep(rgb_image, orig_ext, ".ppm");
        IMG = imread(fullfile(original_images_list(img).folder, original_images_list(img).name));
        imwrite(IMG, fullfile(ppm_folder, ppm_name));
        textprogressbar(img/no_of_images*100);
    end

elseif strcmp(codec_in_file_format, "yuv444p10le")
    orig_yuv_folder = fullfile(dir_struct.ground_truth, codec_in_file_format);
    if(~exist(orig_yuv_folder, 'dir')), mkdir(orig_yuv_folder); end
    dir_struct.set_yuv444p_path = orig_yuv_folder;
    textprogressbar('  => Converting to yuv444p10le : ');
    for img=1:no_of_images
        input_file = fullfile(original_images_list(img).folder, original_images_list(img).name);
        rgb_img = imread(input_file);
        [he, we, ce] = size(rgb_img);
        output_file = fullfile(orig_yuv_folder, sprintf("image_%03d_%dx%d_1Hz_prec10_P444.yuv", img, we,he));
        cmd = sprintf("ffmpeg -hide_banner -i %s -pix_fmt yuv444p10le -vf scale=in_range=full:in_color_matrix=bt709:out_range=full:out_color_matrix=bt709 -color_primaries bt709  -color_trc bt709 -colorspace bt709 -y %s", input_file, output_file);
        [status, out1] = system(cmd);
        textprogressbar(img/no_of_images*100);
    end

elseif strcmp(codec_in_file_format, "JPEG-YCbCr")
    orig_ycbcr_folder = fullfile(dir_struct.ground_truth, "JPEG_YCbCr");
    if(~exist(orig_ycbcr_folder, 'dir')), mkdir(orig_ycbcr_folder); end
    dir_struct.set_jpeg_ycbcr_path = orig_ycbcr_folder;
    textprogressbar('  => Converting to JPEG-YCbCr : ');
    for img=1:no_of_images
        rgb_img = imread(fullfile(original_images_list(img).folder, original_images_list(img).name));
        [he, we, ce] = size(rgb_img);
        ycbcr_img = rgb2jpegycbcr(rgb_img);
        sju_imwrite_yuv(ycbcr_img, fullfile(orig_ycbcr_folder, sprintf("image_%03d_%dx%d_1Hz_prec8_P444.yuv", img, we,he)));
        textprogressbar(img/no_of_images*100);
    end
end

conversion_time = toc(tCstart);
textprogressbar(['   Conversion done in  ' num2str(conversion_time) ' sec']);
exp_params.dir_struct = dir_struct;

end


function exp_params = convertDataset_to_video(exp_params, codec_in_file_format)

tCstart = tic;
orig_ext = exp_params.images_files_ext;
original_images_list = dir(fullfile(exp_params.set_path, orig_ext));
no_of_images = numel(original_images_list);
exp_params.no_of_images = no_of_images;
dir_struct = exp_params.dir_struct;

orig_yuv_folder = fullfile(dir_struct.ground_truth, codec_in_file_format);
if(~exist(orig_yuv_folder, 'dir')), mkdir(orig_yuv_folder); end
dir_struct.video_yuv444p_path = orig_yuv_folder;

input_file = fullfile(original_images_list(1).folder, original_images_list(1).name);

[~, image_name, ~] = fileparts(original_images_list(1).name);
input =fullfile(exp_params.set_path, erase(image_name,"001")+"%03d.tiff");

% input =erase(input_file,"001")+"%03d.ppm";

rgb_img = imread(input_file);
[he, we, ce] = size(rgb_img);
if codec_in_file_format == "yuv444p10le"
    output_file = fullfile(orig_yuv_folder, sprintf("%s_%dx%dx%d_1Hz_prec10_P444.yuv", exp_params.set_name, we,he, no_of_images));
    cmd = sprintf("ffmpeg -hide_banner -i %s -pix_fmt yuv444p10le -vf scale=in_range=full:in_color_matrix=bt709:out_range=full:out_color_matrix=bt709 -color_primaries bt709  -color_trc bt709 -colorspace bt709 -y %s", input, output_file);
elseif codec_in_file_format == "yuv420p10le"
    output_file = fullfile(orig_yuv_folder, sprintf("%s_%dx%dx%d_1Hz_prec10_P420.yuv", exp_params.set_name, we,he, no_of_images));
    cmd = sprintf("ffmpeg -hide_banner -i %s -pix_fmt yuv420p10le -vf scale=in_range=full:in_color_matrix=bt709:out_range=full:out_color_matrix=bt709 -color_primaries bt709  -color_trc bt709 -colorspace bt709 -y %s", input, output_file);
else
    error("Unable to convert to video file for Pixel Format %s", pix_fmt);
end

[status, out1] = system(cmd);


exp_params.dir_struct = dir_struct;
exp_params.yuv_video_file = output_file;

conversion_time = toc(tCstart);
fprintf(['   Conversion done in  ' num2str(conversion_time) ' sec']);
end

function time = extract_time(charVector)
% Split the char vector into lines
lines = strsplit(charVector, '\n');

time = [];  % Initialize the time vector

for i = 1:length(lines)
    line = strtrim(lines{i});  % Remove leading and trailing white spaces
    if startsWith(line, 'Total Time:')  % Check if the line starts with 'Total Time:'
        tokens = strsplit(line);  % Split the line into words
        time = [time, str2double(tokens{3})];  % Add the number after 'Total Time:' to time
    end
end
end



function img_bpp = get_gmis_image_bpp(out1, imgNo)
% Assuming 'text' is your multiline text
lines = splitlines(out1);
if imgNo ==1
    pattern1 = "Compressed Image No :   1";
else
    if imgNo<11
        pattern1 = sprintf("ImageNo : %d  Residue of ", imgNo-1);
    else
        pattern1 = sprintf("ImageNo :%d  Residue of ", imgNo-1);
    end
end

for i = 1:length(lines)
    if contains(lines{i}, pattern1)
        % Extract the line
        line = lines{i};

        % Split the line into words
        words = split(line);

        % The last word is the desired number
        numStr = words{end};
        img_bpp = str2double(numStr);
        break;

    end
end

end