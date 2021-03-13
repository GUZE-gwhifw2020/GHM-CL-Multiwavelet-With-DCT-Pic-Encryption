%% Birth Certificate
% ===================================== %
% DATE OF BIRTH:    2021.03.13
% NAME OF FILE:     CLTest
% FILE OF PATH:     /..
% FUNC:
%    验证CL变换与CLT逆变换。
% ===================================== %
addpath('../')
%% 读取图片
AImg = imread('1.png');
%% 图片延拓
A = zeros([size(AImg,[1 2]) + 4, size(A, 3)]);
A(5:end,5:end,:) = AImg;

%% CL正变换，CL逆变换
E = CLT(A);
RA = CLTInv(E);

%% 结果验证
RAImg = uint8(RA(5:end,5:end,:));
fprintf('\tMSE = %.6e\n', mean((AImg - RAImg).^2, 'all'));
fprintf('\tPixel Diff: %d/%d\n', nnz(AImg ~= RAImg), numel(AImg));
figure
imshow(RAImg)