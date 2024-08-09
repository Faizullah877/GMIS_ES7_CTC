clc
% clear all

M_File_Path = "D:\GMIS_EXPs\EXP8_test\Results";

M_file = fullfile(M_File_Path, "Metrics.mat");
plot_all_codecs = false;

figures_folder = fullfile(M_File_Path, "PlotsPCHIP");
if(~exist(figures_folder, 'dir')), mkdir(figures_folder); end

load(M_file);
set_name = "EPFL_Set4_8b";

set_results_present = isfield(M, set_name);
if(set_results_present ~= 1)
    fprintf("\t\tResults for Set = %s does not exist in the Metrics file.\n", set_name);
    return;
end

MaxBPP = 1.0; 

set_fig_folder = fullfile(figures_folder, set_name);
if(~exist(set_fig_folder, 'dir')), mkdir(set_fig_folder); end

plot_metrics = [ "PSNR_Y", "PSNR_YUV", "PSNR_HVS_M_Y", "MSSSIM_Y"];
% plot_metrics = ["PSNR_RGB", "PSNR_Y",  "PSNR_HVS_M_Y", "MSSSIM_RGB"];

if(plot_all_codecs)
    plot_codecs = string(fieldnames(M.(set_name)));
else
     plot_codecs = ["JPEG1", "JPEG1_Arithmetic", "SJU_Arch_DPCM_QDCT", "SJU_Arch_DPCM_PIXELSv2", "JPEG2000", "SJU_Arch_JPEG2000v2", "JPEGXL", "SJU_Arch_JPEG_XLv2", "JPEG_AI_VM", "SJU_Arch_JPEG_AI", "VVC_VTM_Intra", "VVC_VVenC_Inter"];
    
    % plot_codecs = ["JPEG1", "JPEG1_Arithmetic", "SJU_Arch_DPCM_PIXELSv2", "JPEG2000", "SJU_Arch_JPEG2000v2", "JPEGXL", "SJU_Arch_JPEG_XLv2", "JPEG_AI_VM", "SJU_Arch_JPEG_AI"];  
    % plot_codecs = ["JPEG2000", "SJU_Arch_JPEG2000v2", "JPEGXL", "SJU_Arch_JPEG_XLv2"];    
    % plot_codecs = ["JPEG1", "JPEG1_Arithmetic", "SJU_Arch_DPCM_PIXELSv2", "JPEG2000", "SJU_Arch_JPEG2000v2", "JPEGXL", "SJU_Arch_JPEG_XLv2", "JPEG_AI_VM", "SJU_Arch_JPEG_AI", "VVC_VTM_Intra", "VVC_VVenC_Inter"];
    % plot_codecs = ["JPEG1", "JPEG1_Arithmetic", "JPEG2000", "JPEGXL", "VVC_VVenC_Inter", "VVC_VTM_Intra", "SJU_Arch_DPCM_QDCT", "SJU_Arch_JPEG2000v1", "SJU_Arch_JPEG_XLv1", "SJU_Arch_DPCM_PIXELSv2", "SJU_Arch_JPEG2000v2", "SJU_Arch_JPEG_XLv2", "JPEG_AI_VM"];
    % plot_codecs = ["JPEG1", "JPEG1_Arithmetic", "JPEG2000", "JPEGXL", "VVC_VVenC_Inter", "VVC_VTM_Intra", "SJU_Arch_DPCM_QDCT", "SJU_Arch_JPEG2000", "SJU_Arch_JPEG_XL", "JPEG_AI_VM"];
    % plot_codecs = ["JPEG1", "JPEG1_Arithmetic", "JPEG2000", "JPEGXL", "VVC_VVenC_Inter", "VVC_VVenC_Intra"];
    % plot_codecs = ["JPEG1", "JPEG1_Arithmetic","JPEG2000", "JPEGXL", "VVC_VTM_Intra", "VVC_VVenC_Inter", "JPEG_AI_VM", "SJU_Arch_DPCM_QDCT","SJU_Arch_DPCM_PIXELS"];
    % plot_codecs = ["VVC_VTM_Inter", "VVC_VVenC_Inter",   "VVC_VTM_Intra", "VVC_VVenC_Intra"];
    % plot_codecs = ["JPEG2000", "JPEGXL", "SJU_Arch_DPCM_PIXELSv1", "SJU_Arch_DPCM_PIXELSv2", "SJU_Arch_JPEG2000v1", "SJU_Arch_JPEG2000v2", "SJU_Arch_JPEG_XLv1", "SJU_Arch_JPEG_XLv2"];
    % plot_codecs = ["JPEG1", "JPEG1_Arithmetic", "JPEG2000", "JPEGXL", "SJU_Arch_DPCM_PIXELSv2", "SJU_Arch_JPEG1_DPCMPIXELS_RDOPT"];
end
plot_styles = { ...
    % Codec Name ,          lineStyle ,LineColor, Marker     , plot_name
    {"JPEG1"                 , '-'    , 'k'     , 'o'        , "JPEG 1"               } ...
    {"JPEG1_Arithmetic"      , "-"    , 'y'     , '*'        , "JPEG 1 Arithmethic"   } ...
    {"JPEG2000"              , "-"    , 'g'     , 'square'   , "JPEG 2000"            } ...
    {"JPEGXL"                , "-"    , 'b'     , 'diamond'  , "JPEG XL"              } ...
    {"VVC_VVenC_Inter"       , ":"    , 'm'     , '^'        , "VVC inter"            } ...
    {"VVC_VVenC_Intra"       , "--"   , 'm'     , '>'        , "VVC intra"            } ...
    {"VVC_VTM_Inter"         , "-."   , 'm'     , 'v'        , "VVC inter"            } ...
    {"VVC_VTM_Intra"         , "-"    , 'm'     , 'x'        , "VVC intra"            } ...
    {"SJU_Arch_DPCM_PIXELSv1", "-."   , 'k'     , 'pentagram', "JPEG 1 DPCM-PIXELS-YCbCr"      } ...
    {"SJU_Arch_DPCM_PIXELSv2", ":"    , 'k'     , 'o'        , "JPEG 1 DPCM-PIXELS-RGB"      } ...
    {"SJU_Arch_DPCM_QDCT"    , "-."    , 'k'     , 'o'        , "JPEG 1 DPCM-QDCT"      } ...
    {"SJU_Arch_JPEG2000v1"   , "-."   , 'g'     , 'o'        , "SJU-Arch JPEG 2000"   } ...
    {"SJU_Arch_JPEG2000v2"   , ":"    , 'g'     , '>'        , "SJU-Arch JPEG 2000"   } ...
    {"SJU_Arch_JPEG_XLv1"    , "-."   , 'b'     , '|'        , "SJU-Arch JPEG XL"     } ...
    {"SJU_Arch_JPEG_XLv2"    , ":"    , 'b'     , 'o'        , "SJU-Arch JPEG XL"     } ...
    {"JPEG_AI_VM"            , "-"    , 'r'     , '*'        , "JPEG AI"              } ...
    {"SJU_Arch_JPEG_AI"      , ":"    , 'r'     , '*'        , "SJU-Arch JPEG AI"     } ...
    {"SJU_Arch_JPEG1_DPCMPIXELS_RDOPT", "-.", 'k', '+', "SJU-Arch JPEG1 Rdopt"}
    };
% plot_styles = { ...
%     % Codec Name ,          lineStyle ,LineColor, Marker     , plot_name
%     {"JPEG1"                 , '-'    , 'k'     , 'o'        , "JPEG 1"               } ...
%     {"JPEG1_Arithmetic"      , "-"    , 'y'     , '*'        , "JPEG 1 Arithmethic"   } ...
%     {"JPEG2000"              , "-"    , 'k'     , 'square'   , "JPEG 2000"            } ...
%     {"JPEGXL"                , ":"    , 'k'     , 'diamond'  , "JPEG XL"              } ...
%     {"VVC_VVenC_Inter"       , "-"    , 'g'     , 'o'        , "VVC VVenC inter"            } ...
%     {"VVC_VTM_Inter"         , "-."   , 'b'     , '*'        , "VVC VTM inter"            } ...
%     {"VVC_VVenC_Intra"       , "--"   , 'g'     , '|'        , "VVC VVenC intra"            } ...
%     {"VVC_VTM_Intra"         , ":"    , 'b'     , 'x'        , "VVC VTM intra"            } ...
%     {"SJU_Arch_DPCM_PIXELSv1", "-."   , 'k'     , 'pentagram', "SJU-Arch JPEG 1"      } ...
%     {"SJU_Arch_DPCM_PIXELSv2", ":"    , 'k'     , 'o'        , "SJU-Arch JPEG 1"      } ...
%     {"SJU_Arch_DPCM_QDCT"    , "-"    , 'y'     , 'o'        , "SJU-Arch JPEG 1"      } ...
%     {"SJU_Arch_JPEG2000v1"   , "-."   , 'g'     , 'o'        , "SJU-Arch JPEG 2000"   } ...
%     {"SJU_Arch_JPEG2000v2"   , ":"    , 'g'     , '>'        , "SJU-Arch JPEG 2000"   } ...
%     {"SJU_Arch_JPEG_XLv1"    , "-."   , 'b'     , '|'        , "SJU-Arch JPEG XL"     } ...
%     {"SJU_Arch_JPEG_XLv2"    , ":"    , 'b'     , 'o'        , "SJU-Arch JPEG XL"     } ...
%     {"JPEG_AI_VM"            , "-"    , 'r'     , '*'        , "JPEG AI"              } ...
%     {"SJU_Arch_JPEG_AI"      , ":"    , 'r'     , '*'        , "SJU-Arch JPEG AI"     } ...
%     {"SJU_Arch_JPEG1_DPCMPIXELS_RDOPT", "-.", 'k', '+', "SJU-Arch JPEG1 Rdopt"}
%     };

plot_set_avg_results_interpolated(set_fig_folder, set_name, M.(set_name), plot_metrics, plot_codecs, MaxBPP, plot_styles);
% plot_each_image_results_interpolated(set_fig_folder, set_name, M.(set_name), plot_metrics, plot_codecs,  MaxBPP, plot_styles)

% plot_each_image_results(set_fig_folder, set_name, M.(set_name), plot_metrics, plot_codecs, Markers, MaxBPP);

% plot_metric_distribution(set_fig_folder, set_name, M.(set_name), plot_metrics, plot_codecs, markers);
% 
% plot_bpp_distribution(set_fig_folder, set_name, M.(set_name), plot_codecs, markers);



function plot_set_avg_results_interpolated(set_fig_folder, set_name, M_set, plot_metrics, plot_codecs, MaxBPP, plot_styles)
set_avg_fig_dir = fullfile(set_fig_folder, "Averaged_Overall");
if(~exist(set_avg_fig_dir, 'dir')), mkdir(set_avg_fig_dir); end

no_of_metrics = length(plot_metrics);
for figureNo = 1: no_of_metrics
    metric_name = plot_metrics(figureNo);
    figure_file_name = fullfile(set_avg_fig_dir, sprintf("%s_%s.png",set_name, metric_name));
    no_of_codecs_to_plot = numel(plot_codecs);
    min_values = nan(1, no_of_codecs_to_plot);
    max_values = nan(1, no_of_codecs_to_plot);
    fig = figure('Name',metric_name);
    fig.Position = [200 200 720 512];
    hold on
    for i = 1: no_of_codecs_to_plot
        codec_name = plot_codecs(i);
        matches = cellfun(@(x) strcmp(x{1}, codec_name), plot_styles);
        plot_style =  plot_styles(matches);
        plot_style = plot_style{1};
        linestyle = plot_style{2};
        linecolor = plot_style{3};
        marker = plot_style{4};
        plot_name = plot_style{5};
        if isfield(M_set, codec_name)
            Exist_Column = strcmp(metric_name, M_set.(codec_name).SET_AVG_RESULTS.Properties.VariableNames);
            val = any(Exist_Column);
            if val
                set_avg_bpp = M_set.(codec_name).SET_AVG_RESULTS.SET_BPP;
                bpp_array = set_avg_bpp;
                set_avg_psnr_hvs_m_y = M_set.(codec_name).SET_AVG_RESULTS.(metric_name);
                metric_array = set_avg_psnr_hvs_m_y;
                max_bpp = max(bpp_array);
                min_bpp = min(bpp_array);
                xx1 = min_bpp:0.01:max_bpp;
                xx = [bpp_array' xx1];
                xx = sort(xx);
                yy = pchip(bpp_array, metric_array, xx);
                xx = xx(xx<MaxBPP);
                no_ele = length(xx);
                yy = yy(1:no_ele);
                bpp_array = bpp_array(bpp_array<MaxBPP);
                no_ele = length(bpp_array);
                metric_array = metric_array(1:no_ele);
                min_values(i) = min(metric_array(:));
                max_values(i) = max(metric_array(:));
                plot( xx, yy,'LineStyle', linestyle, 'Color', linecolor, 'LineWidth', 1.5, 'DisplayName',plot_name)
                plot(bpp_array, metric_array, 'Marker', marker, 'MarkerEdgeColor', linecolor, 'LineStyle', 'none', 'HandleVisibility', 'off', 'DisplayName', '');
 
            end
        end
    end
    min_value = mean(min_values, 'omitnan');
    max_value = max(max_values);
    if strcmp(metric_name, "MSSSIM_Y")
        ylim([min_value 1.0]);
    elseif strcmp(metric_name, "PSNR_HVS_M_Y")
        ylim([min_value max_value]);
    else
        ylim([min_value max_value]);
    end           
    title(set_name, 'FontSize', 14, 'Interpreter', 'none');
    xlabel('bpp', 'FontSize',12);
    ylabel(metric_name, 'FontSize',12, 'Interpreter', 'none');
    legend('Location', 'best', 'Interpreter', 'none', 'FontSize',12)
    legend('boxoff');
    grid on
    grid minor
    hold off
    saveas(fig, figure_file_name);
end
close all
end
function plot_each_image_results_interpolated(set_fig_folder, set_name, M_set, plot_metrics, plot_codecs,  MaxBPP, plot_styles)
images_fig_dir = fullfile(set_fig_folder, "Images_Result");
if(~exist(images_fig_dir, 'dir')), mkdir(images_fig_dir); end

codec_list = fieldnames(M_set);
no_of_codecs = numel(codec_list);

no_of_metrics = length(plot_metrics);
for cod = 1: no_of_metrics
    metric_name = plot_metrics(cod);
    image_metric_fig_folder = fullfile(images_fig_dir, metric_name);
    if(~exist(image_metric_fig_folder, 'dir')), mkdir(image_metric_fig_folder); end
end

no_of_images = numel(M_set.(codec_list{1}).IMAGES_NAMES);

for img =1:no_of_images
    % image_name = M_set.(codec_list{1}).IMAGES_NAMES(img);
    [~, image_name, ~] = fileparts(M_set.(codec_list{1}).IMAGES_NAMES(img));
    % image_name = erase(image_name, ".ppm");
    img_col_idx = 7+img;
    for figureNo = 1: no_of_metrics
        metric_name = plot_metrics(figureNo);
        figure_name = fullfile(images_fig_dir, metric_name, sprintf("%s_%s_%s.png",set_name, image_name, metric_name));
        fig = figure('Name',sprintf("%s_%s_%s", set_name, image_name, metric_name));
        fig.Position = [200 200 720 512];
        hold on
        for i = 1: no_of_codecs
            codec_name = codec_list{i};
            matches = cellfun(@(x) strcmp(x{1}, codec_name), plot_styles);
            plot_style =  plot_styles(matches);
            plot_style = plot_style{1};
            linestyle = plot_style{2};
            linecolor = plot_style{3};
            marker = plot_style{4};
            if(codec_is_member(codec_name, plot_codecs))
                headers = M_set.(codec_name).IMAGES_RESULT.BPPs.Properties.VariableNames;
                image_name1 = headers{img_col_idx};
                image_bpp = M_set.(codec_name).IMAGES_RESULT.BPPs.(image_name1);
                bpp_array = image_bpp(image_bpp<MaxBPP);
                no_ele = length(bpp_array);
                image_metrics = M_set.(codec_name).IMAGES_RESULT.METRICS.(metric_name).(image_name1);
                metric_array = image_metrics(1:no_ele);
                plot(bpp_array, metric_array, 'LineStyle', linestyle, 'Color', linecolor, 'LineWidth', 1.5, 'Marker', marker, 'MarkerEdgeColor', linecolor,   'DisplayName',codec_name);
            end
        end
        figure_title = sprintf("%s_%s", set_name, image_name);
        title(figure_title, 'FontSize', 14, 'Interpreter', 'none');
        xlabel('bpp', 'FontSize',12);
        ylabel(metric_name, 'FontSize',12, 'Interpreter', 'none');
        legend('Location', 'best', 'Interpreter', 'none', 'FontSize',12)
        legend('boxoff');
        grid on
        grid minor
        hold off
        saveas(fig, figure_name);
        %         close(fig);
    end
    close all
end
end



function is_member = codec_is_member(codec, codec_list)
no_of_codecs = numel(codec_list);
is_member = false;
for i=1: no_of_codecs
    ab = strcmp(string(codec), string(codec_list{i}));
    if(ab)
        is_member = true;
        return;
    end
end
end