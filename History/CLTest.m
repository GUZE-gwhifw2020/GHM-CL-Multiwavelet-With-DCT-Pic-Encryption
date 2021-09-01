%% Birth Certificate
% ===================================== %
% DATE OF BIRTH:    2021.03.13
% NAME OF FILE:     CLTest
% FILE OF PATH:     /History
% FUNC:
%    验证上级文件夹CL变换与CLT逆变换函数。
% ===================================== %
addpath('../')
%% 读取图片
AImg = imread('1.png');

%% CL正变换，CL逆变换
E = CLT(AImg);
RA = CLTInv(E);

%% 结果验证
RAImg = uint8(RA);
fprintf('\tMSE = %.6e\n', mean((AImg - RAImg).^2, 'all'));
fprintf('\tPixel Diff: %d/%d\n', nnz(AImg ~= RAImg), numel(AImg));
figure(1)
imshow(RAImg)
