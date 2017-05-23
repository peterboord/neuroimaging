function z = corr2z(corrMap,N)
clamp = 0.999;
corrMap(corrMap > clamp) = clamp;
z = atanh(corrMap)*sqrt(N - 3);
end
