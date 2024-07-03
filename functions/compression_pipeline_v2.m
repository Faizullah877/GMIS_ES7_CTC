function M = compression_pipeline_v2(exp_params, codec_params, M)
    %COMPRESSION_PIPELINE Summary of this function goes here
    %   Detailed explanation goes here

    %general parameters
    set_name = exp_params.set_name;
    QA_metrics = exp_params.QA_Metrics;
    do_discard = exp_params.do_discard;   
    codec_exp_time = tic;    
    codec = codec_params{1};
    q_list = codec_params{2};
    data_format = codec_params{3};   
    no_of_QAmetrics = numel(QA_metrics);    
    quality_count = numel(q_list);
    dir_struct = exp_params.dir_struct;  
    no_of_images = exp_params.set_info.no_of_images;    
    set_codec_results = struct();

    % codec data folder
    codec_data_folder = fullfile(dir_struct.set_data_folder, codec);
    if(~exist(codec_data_folder, 'dir')), mkdir(codec_data_folder); end    
    
    % check for existing average results table for a set.
    if isfield(M, set_name) && isfield(M.(set_name), codec) &&  isfield(M.(set_name).(codec), "SET_AVG_RESULTS")
        set_avg_results = M.(set_name).(codec).SET_AVG_RESULTS;
        append_row = true;
    else
        set_avg_results = table('Size', [quality_count 9], 'VariableTypes', ["single", "uint32", "double", "double","double", "double", "double", "double", "double"]);
        variablenames = [["Qvalue", "SET_CompressedSize","SET_BPP", "SET_ENC_TIME", "SET_DEC_TIME"], QA_metrics];
        set_avg_results.Properties.VariableNames = cellstr(variablenames);
        append_row = false;
    end
    

    %% Iterate over rate/q values
    for q_itr = 1: quality_count
        q_v = q_list(q_itr);
        q_fieldName = ['q_' num2str(q_v)];
        q_fieldName = strrep(q_fieldName, ".", "p");

        % check and skip if the data for the q value is already there.
        if isfield(M, set_name)
            if isfield(M.(set_name), codec)
                if isfield(M.(set_name).(codec).q_fields, q_fieldName)
                    str = sprintf("\n\tSet Name : [%s] => Codec : [%s] => Q field : [%s] already exists and the experiment is skipped.\n", set_name, codec, q_fieldName);
                    fprintf(str);
                    continue
                end
            end
        end
    
        % create the folders for the compressed and decompressed data.
        set_q_folder = fullfile(codec_data_folder, q_fieldName);
        if(~exist(set_q_folder, 'dir')), mkdir(set_q_folder); end
        set_compressed_folder = fullfile(set_q_folder, "Compressed");
        if(~exist(set_compressed_folder, 'dir')), mkdir(set_compressed_folder); end
        set_dec_folder = fullfile(set_q_folder, "Decompressed");
        if(~exist(set_dec_folder, 'dir')), mkdir(set_dec_folder); end

        % necessary variables for codec and results
        no_rows = exp_params.set_info.images_details.img_height(1);
        no_cols = exp_params.set_info.images_details.img_width(1);
        ImagesResult_T1 = table('Size', [no_of_images 8], 'VariableTypes', {'uint16', 'string','uint16', 'uint16', 'double', 'single', 'double','double'});
        ImagesResult_T1.Properties.VariableNames = {'ImageNo','ImageName', 'ImageWidth', 'ImageHeight', 'CompressedImageSize','bpp', 'enc_time', 'dec_time'};    
        Images_IQA_T = double(zeros([no_of_images no_of_QAmetrics])); 
    
        % input structure contain data for codec
        input = struct();
        input.set_name = set_name;
        input.codec_folder = dir_struct.codec_folder;
        if strcmp(codec, "JPEG_AI_VM") | strcmp(codec, "SJU_Arch_JPEG_AI")
            input.input_raw = exp_params.set_path;
        else
            input.input_raw = fullfile(dir_struct.ground_truth, data_format);
        end
        input.set_q_folder = set_q_folder;
        input.compressed_folder = set_compressed_folder;
        input.decompressed_folder = set_dec_folder;
        input.q_value = q_v;
        input.width = exp_params.set_info.images_details.img_width(1);
        input.height = exp_params.set_info.images_details.img_height(1);
        input.no_of_images = no_of_images;
    
        fprintf("\tStart Time: %s.\n", get_current_time_str());
        fprintf("\tCompression Pipeline running for SET : [%s {%dx%dx%d}]   CODEC : [%s] and RATE : [%3.2f] => %d/%d\n",set_name,no_cols, no_rows, no_of_images, codec, q_v, q_itr, quality_count);
    
    
        q_start_time = tic();
    
    
    
        if strcmp(codec, "H264_AVC")
            input.config = "inter";
            codec_report = encode_decode_H264_AVC(input);
        elseif strcmp(codec, "H264_AVCintra")
            input.config = "intra";
            codec_report = encode_decode_H264_AVC(input);
        elseif strcmp(codec, "JPEG1")
            input.config = "huff";
            codec_report = encode_decode_JPEG1(input);
        elseif strcmp(codec, "JPEG1_Arithmetic")
            input.config = "arith";
            codec_report = encode_decode_JPEG1(input);
        elseif strcmp(codec, "JPEG2000")
            codec_report = encode_decode_JPEG2000(input);
        elseif strcmp(codec, "JPEGXL")
            codec_report = encode_decode_JPEG_XL(input);
        elseif strcmp(codec, "VVC_VVenC_Inter")
            input.config = "inter";
            codec_report = encode_decode_VVC_VVenC(input);
        elseif strcmp(codec, "VVC_VVenC_Inter_1_thread")
            input.config = "inter";
            input.no_of_threads = 1; 
            codec_report = encode_decode_VVC_VVenC(input);
        elseif strcmp(codec, "VVC_VVenC_Inter_8_thread")
            input.config = "inter";
            input.no_of_threads = 8; 
            codec_report = encode_decode_VVC_VVenC(input);
        elseif strcmp(codec, "VVC_VVenC_Intra")
            input.config = "intra";
            codec_report = encode_decode_VVC_VVenC(input);
        elseif strcmp(codec, "VVC_VTM_Inter")
            input.config = "inter";
            codec_report = encode_decode_VVC_VTM(input);
        elseif strcmp(codec, "VVC_VTM_Intra")
            input.config = "intra";
            codec_report = encode_decode_VVC_VTM(input);
        elseif strcmp(codec, "SJU_Arch_DPCM_QDCT")
            input.encoding_scheme = "jpeg1_dpcm_qdct_seq";
            codec_report = encode_decode_SJU_ARCH(input);
            % codec_report = encode_decode_DPCM_QDCT(input);
        elseif strcmp(codec, "SJU_Arch_DPCM_PIXELSv2")
            input.config = "arith";
            input.codec = "JPEG1";
            codec_report = encode_decode_SJU_ARCHs(input);
            % input.encoding_scheme = "jpeg1_dpcm_pixels_seq_idx";
            % codec_report = encode_decode_SJU_ARCH(input);
        elseif strcmp(codec, "JPEG_AI_VM")
            input.data_folder = codec_data_folder;
            codec_report = organize_and_get_report_JPEG_AI_VM(input);
        elseif strcmp(codec, "SJU_Arch_JPEG_AI")
            input.data_folder = codec_data_folder;
            codec_report = organize_and_get_report_SJU_Arch_JPEG_AI(input);
        elseif strcmp(codec, "SJU_Arch_JPEG2000v2")
            input.codec = "JPEG2000";
            codec_report = encode_decode_SJU_ARCHs(input);
        elseif strcmp(codec, "SJU_Arch_JPEG_XLv2")
            input.codec = "JPEG_XL";
            codec_report = encode_decode_SJU_ARCHs(input);
        else
            fprintf("Codec [%s] functionality not included yet", codec);
            error("Failed");
        end
    
        codec_time = toc(q_start_time);
    
        fprintf('\t\t=> Set Compressed and Decompressed in %f sec.\n', codec_time);
    
        eval_set_start_time = tic;
        textprogressbar('        => Evaluating : ');
        for imgNo=1:no_of_images
            input_raw = fullfile(exp_params.set_info.set_path, exp_params.set_info.images_details.img_names(imgNo));
            Orig_IMG = imread(input_raw);
            [no_rows, no_cols, ~] = size(Orig_IMG);
            [~, image_name, ~] = fileparts(input_raw);
            decoded_ppm_file = fullfile(codec_report.set_dec_folder, sprintf('dec_image_%03d.ppm', imgNo));
            dec_file_exists = isfile(decoded_ppm_file);
            if dec_file_exists
                Recon_IMG = imread(decoded_ppm_file);
                image_IQA_values = evaluate_images_IQA(Orig_IMG, Recon_IMG, QA_metrics);
                if ~isfield(codec_report, "images_result_table")
                    compressed_img_size = codec_report.compressed_file_size/no_of_images;
                    img_bpp = codec_report.set_bpp;
                    img_enc_time = codec_report.set_enc_time/no_of_images;
                    img_dec_time = codec_report.set_dec_time/no_of_images;
                    ImagesResult_T1(imgNo,:) ={imgNo, image_name, no_cols, no_rows, compressed_img_size, img_bpp, img_enc_time, img_dec_time};
                end
                tmp = struct2cell(image_IQA_values);
                v = [tmp{:}];
                Images_IQA_T(imgNo,:) = v;
            else             
                Images_IQA_T(imgNo,:) = NaN;
            end
            textprogressbar(imgNo/no_of_images*100);
        end
        eval_set_time = toc(eval_set_start_time);
        textprogressbar(['   done in  ' num2str(eval_set_time) ' sec']);
    
        if isfield(codec_report, "images_result_table")
            ImagesResult_T1 = codec_report.images_result_table;
        end
    
        fprintf("\tEnd Time: %s.\n\n", get_current_time_str());
    
        m_table = array2table(Images_IQA_T, 'VariableNames',QA_metrics);

        columnNames = m_table.Properties.VariableNames;
        for i = 1:length(columnNames)
            columnName = columnNames{i};
            if isnumeric(m_table.(columnName))
                meanValue = mean(m_table.(columnName), 'omitnan');
                nanIndices = isnan(m_table.(columnName));
                m_table.(columnName)(nanIndices) = meanValue;
            end
        end

        images_results_table = [ImagesResult_T1, m_table];
        if append_row
            existing_rows_count = height(set_avg_results);
            qp_no = existing_rows_count+1;
        else
            qp_no = q_itr;
        end
        set_avg_results(qp_no, :) = {q_v, codec_report.compressed_file_size, codec_report.set_bpp, codec_report.set_enc_time, codec_report.set_dec_time, ...
            mean(images_results_table.(QA_metrics(1))), ...
            mean(images_results_table.(QA_metrics(2))), ...
            mean(images_results_table.(QA_metrics(3))), ...
            mean(images_results_table.(QA_metrics(4)))};
    
        set_codec_results.(q_fieldName) = images_results_table;
        M.(set_name).(codec).q_fields.(q_fieldName) = images_results_table;
    
        M.(set_name).(codec).SET_AVG_RESULTS = sortrows(set_avg_results, "SET_BPP");
        M.(set_name).(codec).RATE_VALUES = q_list;
        M.(set_name).(codec).IMAGES_NAMES = exp_params.set_info.images_details.img_names;
    
        save(fullfile(exp_params.dir_struct.results_folders, exp_params.results_mat_fname), 'M');
    
        % delete directorires and files
        if(do_discard)
            rmdir(set_compressed_folder, "s");
            rmdir(set_dec_folder, "s");
            rmdir(set_q_folder, "s");
        end
    
        M.(set_name).(codec).IMAGES_RESULT = calculate_images_results(M.(set_name).(codec).q_fields, no_of_images, QA_metrics);
    
        save(fullfile(exp_params.dir_struct.results_folders, exp_params.results_mat_fname), 'M');
    end
    codec_time_taken = toc(codec_exp_time);
    fprintf("<strong>  => Experiment Done : => SET = [%s],  => CODEC = [%s]  => TIME TAKEN = [%f sec]</strong>\n\n", set_name, codec, codec_time_taken);
    if(do_discard)
        rmdir(codec_data_folder, "s");
    end
end

function timestr = get_current_time_str()
nowTime = datetime('now', 'Format','dd-MMM-uuuu HH:mm:ss');
timestr = string(nowTime);
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
    q_field_name = qfield_list{qf};
    q_field_name = strrep(q_field_name, 'p', '.');
    numbers = regexp(q_field_name, '\d+', 'match');
    q_value = str2double(numbers);
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
detailed_results.CompressedSizes = sortrows(Compressed_SizesT, "rate");
detailed_results.BPPs = sortrows(bpp_T, "rate");
detailed_results.ENC_TIMES = sortrows(enc_time_T, "rate");
detailed_results.DEC_TIMES = sortrows(dec_time_T, "rate");

metric_results = struct();
no_of_metrics = length(QA_metrics);
for m = 1: no_of_metrics
    metric_v = QA_metrics(m);
    metric_mat = single(zeros(no_of_q_points, no_of_images + extra_columns));
    for qf = 1: no_of_q_points
        metric_mat(qf, 1) = qf;
        q_field_name = qfield_list{qf};
        q_field_name = strrep(q_field_name, 'p', '.');
        numbers = regexp(q_field_name, '\d+', 'match');
        q_value = str2double(numbers);
        % q_value = set_codec_results.(qfield_list{qf}).rate(1);
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
    metric_results.(metric_v) = sortrows(metric_T, "rate");
end

detailed_results.METRICS = metric_results;

end

