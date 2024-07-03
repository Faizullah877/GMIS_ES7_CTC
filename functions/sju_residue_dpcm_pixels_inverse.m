function recon_img = sju_residue_dpcm_pixels_inverse(previous, residue8b, idx1, idx2)
in8_C1 = residue8b(:,:,1);
in8_C2 = residue8b(:,:,2);
in8_C3 = residue8b(:,:,3);
i1_C1 = idx1(:,:,1);
i1_C2 = idx1(:,:,2);
i1_C3 = idx1(:,:,3);
i2_C1 = idx2(:,:,1);
i2_C2 = idx2(:,:,2);
i2_C3 = idx2(:,:,3);

ss = size(in8_C1);
out16_C1 = int16(zeros(ss));
[h, w, ~] = size(residue8b);
for m=1:h
    for n=1:w
        if(i1_C1(m,n) == 1)
            out16_C1(m,n)= 510 - int16(in8_C1(m,n)) -128;
        elseif(i2_C1(m,n) == 1)
            out16_C1(m,n) = 0 - int16(in8_C1(m,n)) - 128;
        else
            out16_C1(m,n) = int16(in8_C1(m,n)) - 128;
        end
    end
end
out16_C2 = int16(zeros(ss));
for m=1:h
    for n=1:w
        if(i1_C2(m,n) == 1)
            out16_C2(m,n)= 510 - int16(in8_C2(m,n)) - 128;
        elseif(i2_C2(m,n) == 1)
            out16_C2(m,n) = 0 - int16(in8_C2(m,n)) - 128;
        else
            out16_C2(m,n) = int16(in8_C2(m,n)) - 128;
        end
    end
end
out16_C3 = int16(zeros(ss));
for m=1:h
    for n=1:w
        if(i1_C3(m,n) == 1)
            out16_C3(m,n)=510 - int16(in8_C3(m,n)) - 128;
        elseif(i2_C3(m,n) == 1)
            out16_C3(m,n) = 0 - int16(in8_C3(m,n)) - 128;
        else
            out16_C3(m,n) = int16(in8_C3(m,n)) - 128;
        end
    end
end
recon_C1 = out16_C1 ;
recon_C2 = out16_C2;
recon_C3 = out16_C3 ;
recon_residue = cat(3, int16(recon_C1), int16(recon_C2), int16(recon_C3));

recon_img = uint8(previous - recon_residue);
end