function hashBit = hashEncode(inputBit)
%HASHENCODE MD5-HASH编码，输出128bit
assert(all(ismember(inputBit(:), [0;1])), ...
    'inputBit不是比特数据');

% bit图片转换为基于ASCII的字符串
charImg = char(sum(reshape(inputBit, 8, []) .* 2.^(7:-1:0)'));
% HASH输出16进制
hashStr = mlreportgen.utils.hash(charImg);
% 进制转换
c = sprintf('%s', dec2bin(hex2dec(char(hashStr)'))');
% 比特转换
hashBit = str2num(c'); %#ok<ST2NM>
end

