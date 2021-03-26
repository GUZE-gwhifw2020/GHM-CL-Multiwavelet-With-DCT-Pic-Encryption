%% Birth Certificate
% ===================================== %
% DATE OF BIRTH:    2021.03.10
% NAME OF FILE:     Main
% FILE OF PATH:     /ForestAndTree
% FUNC:
%   加密图片主工程文件。
% ===================================== %

addpath('../')
clc;clearvars

forestImg = double(imread('Forest1.png'));

%% Step1 载体多小波变换，提取LL1图四个分量
CLResult = CLT(forestImg);
LL1 = CLResult(1:end/2, 1:end/2, 1);

LL2 = LL1(1:end/2, 1:end/2);
LH2 = LL1(1:end/2, end/2+1:end);
HL2 = LL1(end/2+1:end, 1:end/2);
HH2 = LL1(end/2+1:end, end/2+1:end);
clearvars LL1

%% Step2 LL2分量DCT变换，提取高频系数
LL2D = dct(LL2);

%% Step3 HL2,LH2分量DCT变换，提取低频系数
HL2D = dct(HL2);
LH2D = dct(LH2);
CLL12 = double(oddEven(HL2D(1:end/2, 1:end/2)));
CLL13 = double(oddEven(LH2D(1:end/2, 1:end/2)));

C = reshape([CLL12 CLL13], [], 1);

%% Step4 隐藏信息混沌置乱
treeImg = imread('Tree1.bmp');
% CodePre = imbinarize(rgb2gray(treeImg));
CodePre = double(treeImg);
CodePre = CodePre(:);

[s,miu] = chaosSequence(length(CodePre));

neta = 0.2;
xn = s(1); 
miu = miu{1};

CodeH = s > neta;
CInX = xor(CodePre,CodeH)';
% clearvars s
%% Step5 应用优化算法调整CinX
% 暂时无法执行
F = nnz(CInX == C) / numel(C);
CInY = CInX;
%% Step6 最优隐藏结果CInY存入LH2,HL2
% 低bit重写
HL2D(1:end/2, 1:end/2) = oddEvenInsert(HL2D(1:end/2, 1:end/2), reshape(CInY(1:end), size(HL2D(1:end/2, 1:end/2))));
LH2D(1:end/2, 1:end/2) = oddEvenInsert(LH2D(1:end/2, 1:end/2), reshape(CInY(1:end), size(LH2D(1:end/2, 1:end/2))));

%% Step7 校验数据，最优置乱参数，信息HASH存入LL2
% bit图片转换为HASH序列
RL = hashEncode(CodePre);

% LL2分量低bit重写
Temp = LL2D(end/2+1:end,end/2+1:end);
Temp(1:128) = oddEvenInsert(Temp(1:128), RL');
LL2D(end/2+1:end,end/2+1:end) = Temp;
%% Step8 HH2部分lab灰度位平面保存信息HASH
% % beta分量
% HH2Lab = rgb2lab(HH2);
% LbetaG = HH2Lab(:,:,3);
% LbetaG = rgb2gray(repmat(LbetaG(:,:,3),[1 1 3]));
% 
% % Lab中b分量改写
% RH = RL;
% LbetaG = LbetaG - double(oddEven(LH2D(1:N^2/2))) + RH;
% 
% % 恢复
% HH2Lab = cat(HH2Lab(:,:,1:2), LbetaG);

%% Step9 重组图像
LH2 = idct(LH2D);
HL2 = idct(HL2D);
LL2 = idct(LL2D);
LL1 = [LL2 LH2; HL2 HH2];
CLResult(1:end/2, 1:end/2, 1) = LL1;
ImageEncode = CLTInv(CLResult);
T = mean((ImageEncode-forestImg).^2, 'all');

ImageEncode = ImageEncode + 2*randn(size(ImageEncode));
%% Step1 载体多小波变换，提取LL1图四个分量
CLResult_R = CLT(ImageEncode);
LL1_R = CLResult_R(1:end/2, 1:end/2, 1);

LL2_R = LL1_R(1:end/2, 1:end/2);
LH2_R = LL1_R(1:end/2, end/2+1:end);
HL2_R = LL1_R(end/2+1:end, 1:end/2);
HH2_R = LL1_R(end/2+1:end, end/2+1:end);
clearvars LL1_R

%% Step3 HL2,LH2分量DCT变换，提取低频系数
LH2D_R = dct(LH2_R);
HL2D_R = dct(HL2_R);
Code1 = xor(reshape(CodeH,[64 64]), oddEven(LH2D_R(1:end/2, 1:end/2)));
Code2 = xor(reshape(CodeH,[64 64]), oddEven(HL2D_R(1:end/2, 1:end/2)));

%% Step4 是否被攻击
% 提取图片转换为HASH序列
RL_R = hashEncode(Code2);

% LL2分量提取原有HASHI
LL2D_R = dct(LL2_R);
Temp = LL2D_R(end/2+1:end,end/2+1:end);
RL_Ori = oddEven(Temp(1:128))';

if(~isequal(RL_R, RL_Ori))
    fprintf('\t被攻击。\n')
else
    fprintf('\t没有被攻击。\n')
end
%% 绘图
figure(1)
subplot(1,4,1);
imshow(Code1);
title(sprintf('Code1 Error:%.1f%%', mean(Code1 ~= treeImg, 'all') * 100));
subplot(1,4,2);
imshow(Code2);
title(sprintf('Code2 Error:%.1f%%', mean(Code2 ~= treeImg, 'all') * 100));
subplot(1,4,3);
imshow(treeImg);
title('原始信息图像')
subplot(1,4,4);
imshow(uint8(ImageEncode));
title(sprintf('加密载体图像 Error:%.1f', T));