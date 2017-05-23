function C = corrColumns(A,B)

    a = detrend(A,'constant');
    b = detrend(B,'constant');
    C = sum(a.*b)./sqrt(sum(a.*a).*sum(b.*b));
end