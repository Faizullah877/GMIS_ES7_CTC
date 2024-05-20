% function to convert rgb -> JPEG-YCbCr
% code provided by FAIZULLAH <faiz@sju.ac.kr>
function ycbcr_image = conv_rgb2jpegycbcr(rgb_image)
    %rgb_image : uint8 
    %ycbcr_image : uint8
    T = [0.299,0.587,0.114;-0.168736,-0.331264,0.5;0.5,-0.418688,-0.081312]*255;
    temp_img = double(rgb_image)/255;
    ycbcr_image = zeros(size(temp_img));
    ycbcr_image(:,:,1) = T(1)*temp_img(:,:,1)+T(4)*temp_img(:,:,2) + T(7)*temp_img(:,:,3);
    ycbcr_image(:,:,2) = T(2)*temp_img(:,:,1)+T(5)*temp_img(:,:,2) + T(8)*temp_img(:,:,3)+127.5;
    ycbcr_image(:,:,3) = T(3)*temp_img(:,:,1)+T(6)*temp_img(:,:,2) + T(9)*temp_img(:,:,3)+127.5;
    ycbcr_image = uint8(ycbcr_image);
end
