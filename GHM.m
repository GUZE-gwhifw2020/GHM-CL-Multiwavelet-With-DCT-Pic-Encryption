function bP = GHM(aP)
%% GHM GHM多小波变换
%   Input:
%   aP      Input Image, N * N * channel
%   Output:
%   bP      Output Image, N * N * channel

% 图片大小N
N = size(aP,1);
nChannel = size(aP,3);

%% 滤波系数
H0=[3/(5*sqrt(2)),4/5;-1/20,-3/(10*sqrt(2))];
H1=[3/(5*sqrt(2)),0;9/20,1/sqrt(2)];
H2=[0,0;9/20,-3/(10*sqrt(2))];
H3=[0,0;-1/20,0];
G0=[-1/20,-3/(10*sqrt(2));1/(10*sqrt(2)),3/10];
G1=[9/20,-1/sqrt(2);-9/(10*sqrt(2)),0];
G2=[9/20,-3/(10*sqrt(2));9/(10*sqrt(2)),-3/10];
G3=[-1/20,0;-1/(10*sqrt(2)),0];

%% construct the W matrix
% 构造系数矩阵组
w = [H0,H1,H2,H3;G0,G1,G2,G3];
clearvars -except w aP N nChannel
% 构造W矩阵
W = kron(eye(N/2), w(:,1:4)) + kron(circshift(eye(N/2),-1), w(:,5:8));

% 列变换提前
e = eye(2*N);
jj = 0:4:2*N-1;
jj = jj + [1;2];
jj = jj(:);
kk = jj + 2;
e = [e(:,jj) e(:,kk)];
W = e' * W;
clearvars e jj kk

%% 图片处理
bP = zeros(2*N, 2*N, nChannel);

for iChannel = 1:nChannel
    % 单通道图片
    a = squeeze(aP(:,:,iChannel));
    
    %% 诡异的列操作1
    X = zeros(N,2*N);
    % X为原图片行重复，变换为2N*N
    X(:,1:2:2*N) = a;
    X(:,2:2:2*N) = a / sqrt(2);
    aa = W * X';
    
    %% 诡异的行操作2
    X = zeros(2*N,2*N);
    % X为新图片aa列重复，变换为2N*2N
    X(:,1:2:2*N) = aa;
    X(:,2:2:2*N) = aa / sqrt(2);
    b = W * X';
    
    %% 诡异的恢复位置
    % Coefficient permutation
    % b = kron(eye(2),e') * b * kron(eye(2),e);
    b = [b(1:2:N,:);b(2:2:N,:); b((1:2:N)+N,:);b((2:2:N)+N,:)];
    b = [b(:,1:2:N) b(:,2:2:N) b(:,(1:2:N)+N) b(:,(2:2:N)+N)];
    
    bP(:,:,iChannel) = b;
    
end
end