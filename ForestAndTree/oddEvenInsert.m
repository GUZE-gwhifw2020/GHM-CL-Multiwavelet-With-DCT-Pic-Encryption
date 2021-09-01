function y = oddEvenInsert(x, oe)
%ODDEVENINSERT 奇偶嵌入
assert(isequal(size(x), size(oe)), ...
    '输入维度不匹配');
oe = double(oe);
assert(all(ismember(oe(:), [0;1])),...
    'OE输入为BOOL型变量');

c = 0.5 - oe;
y = round(x + c) - oddEven(x + c) + oe;

end

