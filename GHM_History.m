function b = GHM_History(a)
N = size(a,1);
H0 = [3/(5*sqrt(2)),4/5;-1/20,-3/(10*sqrt(2))];
H1 = [3/(5*sqrt(2)),0;9/20,1/sqrt(2)];
H2 = [0,0;9/20,-3/(10*sqrt(2))];
H3 = [0,0;-1/20,0];
G0 = [-1/20,-3/(10*sqrt(2));1/(10*sqrt(2)),3/10];
G1 = [9/20,-1/sqrt(2);-9/(10*sqrt(2)),0];
G2 = [9/20,-3/(10*sqrt(2));9/(10*sqrt(2)),-3/10];
G3 = [-1/20,0;-1/(10*sqrt(2)),0];
% construct the W matrix
w = [H0,H1,H2,H3;G0,G1,G2,G3];
for i = 1:N/2-1
    W(4*(i-1)+1:4*i,4*i-3:4*i+4) = w;
end
W = [W;[[H2,H3;G2,G3],zeros(4,2*N-8),[H0,H1;G0,G1]]];
p = [];X = [];
% oversampling rows (repeated row preprocessing)
X(1:2:2*N,:) = a;X(2:2:2*N,:) = 1/sqrt(2)*a;
% row transformation
z = W*X;
% row vector permutation
ii = 0:4:2*N-1;jj = sort([ii+1,ii+2]);kk = sort([ii+3,ii+4]);
p = [p;z(jj,:);z(kk,:)];
aa = p';X = [];
% oversampling columns (repeated row preprocessing)
X(1:2:2*N,:) = aa;X(2:2:2*N,:) = 1/sqrt(2)*aa;
% column transformation
z = W*X;
% column vector permutation
p = [];
ii = 0:4:2*N-1;jj = sort([ii+1,ii+2]);kk = sort([ii+3,ii+4]);
p = [p;z(jj,:);z(kk,:)];
b = p';
% Coefficient permutation
b1 = b(1:N,1:N);b2 = b(1:N,N+1:2*N);b3 = b(N+1:2*N,1:N);b4 = b(N+1:2*N,N+1:2*N);
T = [b1(1:2:N,:);b1(2:2:N,:)]';b1 = [T(1:2:N,:);T(2:2:N,:)]';
T = [b2(1:2:N,:);b2(2:2:N,:)]';b2 = [T(1:2:N,:);T(2:2:N,:)]';
T = [b3(1:2:N,:);b3(2:2:N,:)]';b3 = [T(1:2:N,:);T(2:2:N,:)]';
T = [b4(1:2:N,:);b4(2:2:N,:)]';b4 = [T(1:2:N,:);T(2:2:N,:)]';
b = [b1,b2;b3,b4];
end