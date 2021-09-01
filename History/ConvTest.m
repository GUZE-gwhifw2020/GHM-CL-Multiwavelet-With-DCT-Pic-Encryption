%% Birth Certificate
% ===================================== %
% DATE OF BIRTH:    2021.03.25
% NAME OF FILE:     ConvTest
% FILE OF PATH:     /History
% FUNC:
%   卷积，循环卷积，解卷积。
% ===================================== %
%% 测试数据
A = [1 2 0 1 1 0 0 1 1 1];
C = [0.2 0.5 -0.1 0.2];

%% 普通卷积
N = length(A) + length(C) - 1;
%% 卷积
K1 = conv(A,C);
RA1 = deconv(K1,C);
%% 傅里叶变换
K2 = ifft(fft(A, N) .* fft(C, N));
RA2 = ifft(fft(K2, N) ./ fft(C, N));
%% 循环表达
K3 = zeros(1, N);
ATemp = [zeros([1, length(C) - 1]) A zeros([1, length(C) - 1])]; 
for ii = 1:N
    K3(ii) = sum(ATemp(ii:ii+length(C)-1) .* fliplr(C));
end
%% 循环卷积
N = length(A);
% 傅里叶变化
K4 = ifft(fft(A) .* fft(C, N));
RA3 = ifft(fft(K4, N) ./ fft(C, N));
%% 循环表达
K5 = zeros(1, N);
for ii = 1:N
    K5(ii) = sum(A(mod((ii-length(C)+1:ii)-1, N) + 1) .* fliplr(C));
end
%% 解循环卷积 - 法1
CInv = ifft(1./fft(C, N));
RA5 = zeros(1, N);
for ii = 1:N
    RA5(ii) = sum(K5(mod((ii-N+1:ii)-1, N) + 1) .* fliplr(CInv));
end
%% 解循环卷积 - 法2
Cmat = zeros(length(A));
for ii = 1:N
    Cmat(ii, mod((ii-length(C)+1:ii)-1,N)+1) = fliplr(C);
end
RA6 = Cmat \ K5';

