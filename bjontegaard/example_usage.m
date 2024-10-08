%% Example file for calculation of the BD metric

% BSD 3-Clause License
%
% Copyright (c) 2022-2023, Friedrich-Alexander-Universität Erlangen-Nürnberg.
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
%
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution.
%
% * Neither the name of the copyright holder nor the names of its
%   contributors may be used to endorse or promote products derived from
%   this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


% addpath('bjontegaard');


%% Test data for test case 2 (low bitrates)
rate_anchor = [9487.76, 4593.60, 2486.44, 1358.24];
psnr_anchor = [40.037, 38.615, 36.845, 34.851];
rate_test = [9787.80, 4469.00, 2451.52, 1356.24];
psnr_test = [40.121, 38.651, 36.970, 34.987];

bdr_pchip = bd_pchip(rate_anchor, psnr_anchor, rate_test, psnr_test);
bdr_akima = bd_akima(rate_anchor, psnr_anchor, rate_test, psnr_test);

disp('Low bitrates: ');
disp(['PCHIP BDR: ' num2str(bdr_pchip*100,'%4.16f') '%']);
disp(['Akima BDR: ' num2str(bdr_akima*100,'%4.16f') '%']);

%% Test data for test case 2 (high bitrates)
rate_anchor = [82472.76, 26875.24, 9487.76, 4593.6];
psnr_anchor = [43.825, 41.404, 40.037, 38.615];
rate_test = [99896.96, 33952.04, 9787.8, 4469.0];
psnr_test = [44.731, 41.875, 40.121, 38.651];

bdr_pchip = bd_pchip(rate_anchor, psnr_anchor, rate_test, psnr_test);
bdr_akima = bd_akima(rate_anchor, psnr_anchor, rate_test, psnr_test);

disp('High bitrates: ');
disp(['PCHIP BDR: ' num2str(bdr_pchip*100,'%4.16f') '%']);
disp(['Akima BDR: ' num2str(bdr_akima*100,'%4.16f') '%']);

%% Test data for test case 4 (two supporting points)
rate_anchor = [82472.76, 26875.24];
psnr_anchor = [43.825, 41.404];
rate_test = [99896.96, 33952.04];
psnr_test = [44.731, 41.875];

bdr_pchip = bd_pchip(rate_anchor, psnr_anchor, rate_test, psnr_test);
bdr_akima = bd_akima(rate_anchor, psnr_anchor, rate_test, psnr_test);

disp('Two supporting points: ');
disp(['PCHIP BDR: ' num2str(bdr_pchip*100,'%4.16f') '%']);
disp(['Akima BDR: ' num2str(bdr_akima*100,'%4.16f') '%']);



%% Test data for test case 4 (six supporting points + RCD curve visualization)
rate_anchor = [82472.76, 26875.24, 9487.76, 4593.6, 2486.44, 1358.24];
psnr_anchor = [43.825, 41.404, 40.037, 38.615, 36.845, 34.851];
rate_test = [99896.96, 33952.04, 9787.8, 4469.0, 2451.52, 1356.24];
psnr_test = [44.731, 41.875, 40.121, 38.651, 36.970, 34.987];

bdr_pchip = bd_pchip(rate_anchor, psnr_anchor, rate_test, psnr_test, true);
bdr_akima = bd_akima(rate_anchor, psnr_anchor, rate_test, psnr_test, true);

disp('Six supporting points: ');
disp(['PCHIP BDR: ' num2str(bdr_pchip*100,'%4.16f') '%']);
disp(['Akima BDR: ' num2str(bdr_akima*100,'%4.16f') '%']);
