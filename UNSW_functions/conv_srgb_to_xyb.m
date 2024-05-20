% convert from rgb to xyb using reference C code
function [outxyb] = conv_srgb_to_xyb(inrgb)
    addpath('pfm');

    ppm_folder = 'data\';
    prefix_name = 'inRGB';
    in_ppm_ref_name = [ppm_folder prefix_name '.ppm'];
    imwrite(uint8(inrgb), in_ppm_ref_name, 'MaxValue', 255);
    
    
    % read ppm and convert to XYB and save as pfm
    prefix_name = 'outXYB';
    out_pfm_ref_name = [ppm_folder prefix_name '.pfm'];
    
    EXE_DIR = 'C:\Users\reji_\MatlabProjects\xyb_colorspace\xyb_converter-main\xybConversion\x64\Debug\';
    CONV_XYB = ['xybConversion.exe ' in_ppm_ref_name  ' ' out_pfm_ref_name];
    EXE_COMMAND = [EXE_DIR CONV_XYB];
    [status_conv,cmdout_conv] = system(EXE_COMMAND);

    outxyb = read_pfm(out_pfm_ref_name, 0);
end