function recon_residue = convert_8bit_residue_to_16bit_v4(residue8bit, idx)
in8_C1  = residue8bit(:,:,1);
in8_C2 = residue8bit(:,:,2);
in8_C3 = residue8bit(:,:,3);
ss = size(in8_C1);
out16_C1 = int16(zeros(ss));
[h, w, ~] = size(residue8bit);
for m=1:h
    for n=1:w
        if(idx.i1_C1(m,n) == 1)
            out16_C1(m,n)= 384 - int16(in8_C1(m,n));
        elseif(idx.i2_C1(m,n) == 1)
            out16_C1(m,n) = 0 - int16(in8_C1(m,n));
        else
            out16_C1(m,n) = int16(in8_C1(m,n));
        end
    end
end
out16_C2 = int16(zeros(ss));
for m=1:h
    for n=1:w
        if(idx.i1_C2(m,n) == 1)
            out16_C2(m,n)= 384 - int16(in8_C2(m,n));
        elseif(idx.i2_C2(m,n) == 1)
            out16_C2(m,n) = 0 - int16(in8_C2(m,n));
        else
            out16_C2(m,n) = int16(in8_C2(m,n));
        end
    end
end
out16_C3 = int16(zeros(ss));
for m=1:h
    for n=1:w
        if(idx.i1_C3(m,n) == 1)
            out16_C3(m,n)=384 - int16(in8_C3(m,n));
        elseif(idx.i2_C3(m,n) == 1)
            out16_C3(m,n) = 0 - int16(in8_C3(m,n));
        else
            out16_C3(m,n) = int16(in8_C3(m,n));
        end
    end
end
%ch_size = size(y_cap_Y);
recon_C1 = (out16_C1 * 2) + 255; % recon_C1 = (out16_C1 + 128)*2;
recon_C2 = (out16_C2 * 2) + 255; % recon_C2 = (out16_C2 + 128)*2;
recon_C3 = (out16_C3 * 2) + 255; % recon_C3 = (out16_C3 + 128)*2;
recon_residue = cat(3, uint16(recon_C1), uint16(recon_C2), uint16(recon_C3));

end

