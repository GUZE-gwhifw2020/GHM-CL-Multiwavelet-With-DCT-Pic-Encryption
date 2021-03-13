%% Birth Certificate
% ===================================== %
% DATE OF BIRTH:    2021.03.10
% NAME OF FILE:     CodeInsert
% FILE OF PATH:     /ForestAndTree
% FUNC:
%   DCT系数奇偶嵌入实例。
% ===================================== %
A = imread('../1.png');
addpath('../')
%% Encode

% CL变换
CLA = CLT(A); N = size(CLA,1);

% 暂时只取第一通道LL1-HL2分量
indX = 1:N/4; indY = N/4+1:N/2;
CLA1 = squeeze(CLA(indX, indY, 1));
N = size(CLA1,1);

% DCT变换
CLA1DCT = dct(CLA1);

% 待嵌入二进制文件
codeInsert = randi([0,1],size(CLA1)/2);

% 仅在DCT高频分量添加，下标索引
indS = N/2+1:N;

% 嵌入比特(有待优化，目前的嵌入方式过于暴力)
CLA1DCT(indS,indS) = round(CLA1DCT(indS,indS)/2)*2 + codeInsert;

% DCT逆变换
CLA1 = idct(CLA1DCT);

% CL第一通道重写
CLA(indX, indY, 1) = CLA1;

% CL逆变换
RA = CLTInv(CLA);
%% Decode
CLRA = CLT(RA);
CLRA1 = squeeze(CLRA(indX, indY, 1));
CLRA1DCT = dct(CLRA1);
%% 破坏图片Decode
RADamg = RA;
RADamg(40:60,50:200,:) = 0; RADamg(100:250,100:130,:) = 0;
CLRA = CLT(RADamg);
CLRA2 = squeeze(CLRA(indX, indY, 1));
CLRA2DCT = dct(CLRA2);
%% 原始图片与重写后图片
figure
subplot(1,3,1); imshow(A);
subplot(1,3,2); imshow(uint8(RA));
subplot(1,3,3); imshow(uint8(RADamg));
%% 原始二进制文件与解密二进制文件
figure
subplot(1,4,1); imshow(codeInsert)
subplot(1,4,2); imshow(double(oddEven(CLRA1DCT(indS,indS))));
subplot(1,4,3); imshow(double(oddEven(CLRA2DCT(indS,indS))));
subplot(1,4,4); imshow(codeInsert == double(oddEven(CLRA2DCT(indS,indS))));
%%
% 比特恢复率
nnz(codeInsert == oddEven(CLRA1DCT(indS,indS))) / numel(codeInsert)
nnz(codeInsert == oddEven(CLRA2DCT(indS,indS))) / numel(codeInsert)