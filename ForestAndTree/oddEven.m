function y = oddEven(x)
%ODDEVEN 奇偶判断至小数
%   x       Odd/Even
%   0       E
%   1       O
%   -1      O
%   0.1     E
%   -0.1    E
%   0.51    O
%   0.49    E
%   0.5     O

y = (x/2+0.25-round(x/2+0.25)) <= 0;
end

