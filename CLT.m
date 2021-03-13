function Y = CLT(X)
%CLT CL Multi-Wavelet Transform
%   CL多小波变换
Y = zeros(size(X));
N = size(X, 1);

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
Pre = [1/4 1/4;1/(1+sqrt(7)) -1/(1+sqrt(7))];
% 后置滤波器
% Post = [2 (1+sqrt(7))/2;2 -(1+sqrt(7))/2];

% 逐通道处理
for iChannel = 1:size(X, 3)
    A = squeeze(double(X(:,:,iChannel)));
    
    % Step1-1 行前置滤波
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
    % Step1-2 列前置滤波
    C = zeros(N);
    for ii = 1:N
        for nn = 1:N/2
            Bcol = [B(2*nn-1,ii);B(2*nn,ii)];
            Ccol = Pre * Bcol;
            C(nn,ii) = Ccol(1);
            C(N/2+nn,ii) = Ccol(2);
        end
    end
    % Step 2-1 行多小波正变换
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
    % Step 2-2 列多小波正变换
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
    Y(:,:,iChannel) = E;
    
end
end

