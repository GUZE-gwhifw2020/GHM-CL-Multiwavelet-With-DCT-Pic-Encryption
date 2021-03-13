%% Birth Certificate
% ===================================== %
% DATE OF BIRTH:    2021.03.03
% NAME OF FILE:     CL
% FILE OF PATH:     /..
% FUNC:
%   CL变换测试文件。
% ===================================== %

clc;clearvars;close all;
%% 图片
aP = double(imread('1.png')) / 256;
tic

% 图片大小N
N = size(aP,1) + 4;
nChannel = size(aP,3);

%% 滤波系数
% 低通滤波器
L0 = [1/2 -1/2;sqrt(7)/4 -sqrt(7)/4] / sqrt(2);
L1 = [1 0;0 1/2] / sqrt(2);
L2 = [1/2 1/2;-sqrt(7)/4 -sqrt(7)/4] / sqrt(2);

% 高通滤波器
H0 = [1/2 -1/2;-1/4 1/4] / sqrt(2);
H1 = [-1 0;0 sqrt(7)/2] / sqrt(2);
H2 = [1/2 1/2;1/4 1/4] / sqrt(2);

L = {L0;L1;L2};
H = {H0;H1;H2};
clearvars L0 L1 L2 H0 H1 H2
% 前置滤波器
Pre = [1/4 1/4;1/(1+sqrt(7)) -1/(1+sqrt(7))];
% 后置滤波器
Post = [2 (1+sqrt(7))/2;2 -(1+sqrt(7))/2];

%% Step0 通道提取与延拓
A = zeros(size(aP,[1 2]) + 4);
A(5:end,5:end) = squeeze(aP(:,:,1));

%% Step1-1 行前置滤波
B = zeros(N);
for ii = 1:N
    for nn = 1:N/2
        Arow = [A(ii,2*nn-1);A(ii,2*nn)];
        % 公式中k=0，滤波只有一项
        Brow = Pre * Arow;
        B(ii,nn) = Brow(1);
        B(ii,N/2+nn) = Brow(2);
    end
end
clearvars Arow Brow
%% Step1-2 列前置滤波
C = zeros(N);
for ii = 1:N
    for nn = 1:N/2
        Bcol = [B(2*nn-1,ii);B(2*nn,ii)];
        Ccol = Pre * Bcol;
        C(nn,ii) = Ccol(1);
        C(N/2+nn,ii) = Ccol(2);
    end
end
clearvars Bcol Ccol
%% Step 2-1 行多小波正变换
D = zeros(N);
for mm = 1:N/4
    for ii = 1:N
        DrowL = [0;0];
        DrowH = [0;0];
        
        for nn = (2*mm):(2*mm+2)
            if(N/2+nn <= N)
                Crow = [C(ii,nn);C(ii,N/2+nn)];
            else
                % Crow = [C(ii,nn);C(ii,mod(N/2+nn-1,N)+1)];
                Crow = [C(ii,nn);0];
            end
            DrowL = DrowL + L{nn-2*mm+1}*Crow;
            DrowH = DrowH + H{nn-2*mm+1}*Crow;
        end
        D(ii,mm) = DrowL(1);
        D(ii,N/4+mm) = DrowL(2);
        D(ii,N/2+mm) = DrowH(1);
        D(ii,N/2+N/4+mm) = DrowH(2);
    end
end
clearvars Crow DrowL DrowH
%% Step 2-2 列多小波正变换
E = zeros(N);
for mm = 1:N/4
    for ii = 1:N
        ErowL = [0;0];
        ErowH = [0;0];
        
        for nn = (2*mm):(2*mm+2)
            if(N/2+nn <= N)
                Dcol = [D(nn,ii);D(N/2+nn,ii)];
            else
                % Dcol = [D(nn,ii);D(mod(N/2+nn-1,N)+1,ii)];
                Dcol = [D(nn,ii);0];
            end
            ErowL = ErowL + L{nn-2*mm+1}*Dcol;
            ErowH = ErowH + H{nn-2*mm+1}*Dcol;
        end
        E(mm,ii) = ErowL(1);
        E(N/4+mm,ii) = ErowL(2);
        E(N/2+mm,ii) = ErowH(1);
        E(N/2+N/4+mm,ii) = ErowH(2);
    end
end
clearvars ErowL ErowH Dcol 

%% Step 1-1 列多小波逆变换
RD = zeros(N);
for nn = 1:N/2
    for ii = 1:N
        RDcolL = [0;0];
        RDcolH = [0;0];
        for mm = 1:N/4
            if(nn-2*mm >= 0 && nn-2*mm <=2)
                
                EcolL = [E(mm,ii);E(N/4+mm,ii)];
                EcolH = [E(N/2+mm,ii);E(3*N/4+mm,ii)];
                
                RDcolL = RDcolL + L{nn-2*mm+1}' * EcolL;
                RDcolH = RDcolH + H{nn-2*mm+1}' * EcolH;
            end
        end
        RD(nn,ii) = RDcolL(1) + RDcolH(1);
        RD(N/2+nn,ii) = RDcolL(2) + RDcolH(2);
    end
end
clearvars RDcolL RDcolH EcolL EcolH
% RD(1,:) = RD(2,:);
%% Step 1-2 行多小波逆变换
RC = zeros(N);
for nn = 1:N/2
    for ii = 1:N
        RCrowL = [0;0];
        RCrowH = [0;0];
        
        for mm = 1:N/4
            if(nn-2*mm >= 0 && nn-2*mm <= 2)
                RDrowL = [RD(ii,mm);RD(ii,N/4+mm)];
                RDrowH = [RD(ii,N/2+mm);RD(ii,3*N/4+mm)];
        
                RCrowL = RCrowL + L{nn-2*mm+1}' * RDrowL;
                RCrowH = RCrowH + H{nn-2*mm+1}' * RDrowH;
            end
        end
        RC(ii,nn) = RCrowL(1) + RCrowH(1);
        RC(ii,nn+N/2) = RCrowL(2) + RCrowH(2);
    end
end
clearvars RCrowL RCrowH RDrowL RDrowH
% RC(:,1) = RC(:,2);
%% Step 2-1 列后置滤波
RB = zeros(N);
for ii = 1:N
    for nn = 1:N/2
        RCcol = [RC(nn,ii);RC(nn+N/2,ii)];
        RBcol = Post * RCcol;
        RB(2*nn-1,ii) = RBcol(1);
        RB(2*nn,ii) = RBcol(2);
    end
end
clearvars RCcol RBcol
%% Step 2-2 行后置滤波
RA = zeros(N);
for ii = 1:N
    for nn = 1:N/2
        RBrow = [RB(ii,nn);RB(ii,nn+N/2)];
        RArow = Post * RBrow;
        RA(ii,2*nn-1) = RArow(1);
        RA(ii,2*nn) = RArow(2);
    end
end
clearvars RBrow RArow

%%
toc
figure('Name','结果');
subplot(3,3,1);
imshow(A);xlabel('原始图像 - A')
subplot(3,3,2);
imshow(B);xlabel('行前置滤波 - B')
subplot(3,3,3);
imshow(C);xlabel('列前置滤波 - C')
subplot(3,3,4);
imshow(D);xlabel('行多小波正变换 - D')
subplot(3,3,5);
imshow(E);xlabel('列多小波正变换 - E')
subplot(3,3,6);
imshow(RD);xlabel('列多小波逆变换 - RD')
subplot(3,3,7);
imshow(RC);xlabel('行多小波逆变换 - RC')
subplot(3,3,8);
imshow(RB);xlabel('列后置滤波 - RB')
subplot(3,3,9);
imshow(RA);xlabel('行后置滤波 - RA')
%%
A = A(5:end,5:end);
RA = RA(5:end,5:end);
fprintf('\tMSE = %.6e\n', mse(RA,A))
fprintf('\tPixel Diff: %d/%d\n', nnz(round(A*256) ~= round(RA*256)), numel(RA));
figure
imshow(RA)

