function [s,parm] = chaosSequence(N, stype)
%CHAOSSEQUENCE 混沌置乱序列生成器
if(nargin < 2)
    stype = 'cheby';
end

s = zeros(N,1);

if(strcmpi(stype, 'cheby'))
    s(1) = rand(1) * 2 - 1;
    parm = {5};
    for iter = 2:N
        s(iter) = cos(parm{1} * acos(s(iter-1)));
    end
end

end