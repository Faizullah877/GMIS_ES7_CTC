clc
clear all

date_time_str = datestr(datetime('now', 'Format','yyyy-MM-dd HH:mm:ss.SSS'));
date_time_str = erase(date_time_str, ["-", ":", " "]);

diary_name = strcat("diary_", date_time_str,".txt");
diary(diary_name);
diary on

%% Set Experiment parameters
exp_params = struct();
exp_params.set_name = "TestPlant_Teapot";
exp_params.images_files_ext = "*.tiff"; % * must be provided.  
exp_params.set_path = fullfile("D:\GMIS_EXPs\Test_CTC\dataset", exp_params.set_name); 
exp_params.experiment_folder = "D:\GMIS_EXPs\Test_CTC";
exp_params.results_mat_fname = 'Metrics.mat';
exp_params.QA_Metrics = ["PSNR_Y", "PSNR_YUV", "PSNR_HVS_M_Y", "MSSSIM_Y"];

convert_dataset = true;

exp_params.Anchors_parameters = {...
    % {'JPEG1', [15, 20, 30, 40], "PPM", true}... %[7, 20, 40, 60, 80, 85]
    % {'JPEG1_Arithmetic', [10, 25, 15, 35, 45], "PPM", convert_dataset}... %[7, 20, 40, 60, 80, 85]
    % {'SJU_Arch_DPCM_QDCT', [10, 25, 50, 75, 90, 95], "JPEG-YCbCr"}... %[10, 30, 50, 70, 85]
    % {'JPEG2000', [10, 25, 35, 50], "PPM", convert_dataset}... %[7, 20, 40, 60, 80, 85]
    % {'JPEGXL', [10, 25, 35, 50], "PPM", convert_dataset}.... %[10, 30, 50, 70, 85, 90]
    % {'VVC_VVenC_Inter', [50, 35, 25, 20], "yuv420p10le", true}... %[50, 35, 30, 25, 20, 15, 10, 5],
    % {'VVC_VVenC_Intra', [50, 35, 30, 25], "yuv420p10le", true}... % [50, 35, 30, 25, 20, 15, 10, 5]
    % {'VVC_VTM_Inter', [50, 35], "yuv420p10le", true}... % [50, 35, 25, 20, 15, 12, 10, 8]
    {'VVC_VTM_Intra', [50, 35], "yuv420p10le", true}... % [50, 35, 25, 20, 15, 13, 11, 9]
%         {'JPEG1arith_rdopt', [0.15, 0.30, 0.50, 0.75, 0.80, 0.90, 1.0]}...
%         {'JPEG1huff_rdopt', [0.15, 0.30, 0.50, 0.75, 0.80, 0.90, 1.0]}...
%         {'SJU_Arch_DPCM_PIXs_rdopt1', [0.15, 0.30, 0.50, 0.8, 1.0], "JPEG-YCbCr"}...
%         {'SJU_Arch_JPEG1_DCT_SEQ', [10, 30, 50, 70, 85], "JPEG-YCbCr"}...
%         {'UNSW_SJU_JPEG1arith_rdopt_v4', [0.15, 0.30, 0.50], "JPEG-YCbCr"}...
%         {'UNSW_SJU_JPEG1arith_v4', [7, 20, 40, 60, 80, 85], "JPEG-YCbCr"}...
%           {'SJU_Arch_JPEG1arith_rdopt1', [0.15, 0.30, 0.50], "PPM"}...
    % {'H264_AVC', [25, 30], "*.tiff"};
%     {'H264_AVCintra', [25, 30], 'PPM'}
    };

exp_params.do_discard = false;

%% Add paths to functions etc
addpath(genpath('./functions'));
addpath(genpath('./GMIS_IQA'));
addpath(genpath('./UNSW_v4_functions'));

exp_params.dir_struct = prepare_directories(exp_params);


if exist(fullfile(exp_params.dir_struct.results_folders, exp_params.results_mat_fname), 'file')
    load(fullfile(exp_params.dir_struct.results_folders, exp_params.results_mat_fname));
else
    M=struct();
end



%% loop through each anchor
no_of_anchors = numel(exp_params.Anchors_parameters);

for ancNo = 1:no_of_anchors
    anch_P = exp_params.Anchors_parameters{ancNo}; 
    M = compression_pipeline_v2(exp_params, anch_P, M);
    save(fullfile(exp_params.dir_struct.results_folders, exp_params.results_mat_fname), 'M');
end

if exp_params.do_discard
    rmdir(exp_params.dir_struct.set_data_folder, "s");
end
diary off

%% create the necessary directories. 
function dir_struct = prepare_directories(exp_params)
    
    current_dir = pwd;
    dir_struct.set_path = exp_params.set_path;
    dir_struct.codec_folder = fullfile(current_dir, "codecs");
    dir_struct.gmis_iqa = fullfile(current_dir, "GMIS_IQA");
    
    dir_struct.experiment_folder = exp_params.experiment_folder;
    if(~exist(dir_struct.experiment_folder, 'dir')), mkdir(dir_struct.experiment_folder); end

    dir_struct.results_folders = fullfile(dir_struct.experiment_folder, "Results");
    if(~exist(dir_struct.results_folders, 'dir')), mkdir(dir_struct.results_folders); end
      
    dir_struct.set_data_folder = fullfile(dir_struct.experiment_folder, exp_params.set_name);
    if(~exist(dir_struct.set_data_folder, 'dir')), mkdir(dir_struct.set_data_folder); end
    
    dir_struct.ground_truth = fullfile(dir_struct.set_data_folder, "GroundTruth");
    if(~exist(dir_struct.ground_truth, 'dir')), mkdir(dir_struct.ground_truth); end
end