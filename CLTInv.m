function X = CLTInv(Y)
%CLTINV CL Multi-Wavelet Inverse Transform
%   此处显示详细说明

X = zeros(size(Y));
N = size(Y,1);

% 滤波器系数
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
% 前置滤波器
% Pre = [1/4 1/4;1/(1+sqrt(7)) -1/(1+sqrt(7))];
% 后置滤波器
Post = [2 (1+sqrt(7))/2;2 -(1+sqrt(7))/2];

% 主通道处理
for iChannel = 1:size(Y, 3)
    E = squeeze(Y(:,:,iChannel));
    RD = zeros(N);
    for nn = 1:N/2
        for ii = 1:N
            RDcolL = [0;0];
            RDcolH = [0;0];
            for kk = 0:2
                mm = (nn - kk)/2;
                if(mod(mm, 1) ~= 0)
                    continue
                end
                if(mm <= 0)
                    mm = mm + N/4;
                end
                EcolL = [E(mm,ii);E(N/4+mm,ii)];
                EcolH = [E(N/2+mm,ii);E(3*N/4+mm,ii)];
                RDcolL = RDcolL + L{kk+1}' * EcolL;
                RDcolH = RDcolH + H{kk+1}' * EcolH;
            end
            RD(nn,ii) = RDcolL(1) + RDcolH(1);
            RD(N/2+nn,ii) = RDcolL(2) + RDcolH(2);
        end
    end
    %% Step 1-2 行多小波逆变换
    RC = zeros(N);
    for nn = 1:N/2
        for ii = 1:N
            RCrowL = [0;0];
            RCrowH = [0;0];
            for kk = 0:2
                mm = (nn-kk)/2;
                if(mod(mm, 1) ~= 0)
                    continue
                end
                if(mm <= 0)
                    mm = mm + N/4;
                end
                RDrowL = [RD(ii,mm);RD(ii,N/4+mm)];
                RDrowH = [RD(ii,N/2+mm);RD(ii,3*N/4+mm)];
                RCrowL = RCrowL + L{kk+1}' * RDrowL;
                RCrowH = RCrowH + H{kk+1}' * RDrowH;
            end
            RC(ii,nn) = RCrowL(1) + RCrowH(1);
            RC(ii,nn+N/2) = RCrowL(2) + RCrowH(2);
        end
    end
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
    X(:,:,iChannel) = RA;
end
end

