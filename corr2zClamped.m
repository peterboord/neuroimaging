function z = corr2zClamped(N,corrMap)
clamp = 0.999;
corrMap(corrMap > clamp) = clamp;
corrMap(corrMap < -clamp) = -clamp;
z = single(atanh(corrMap)*sqrt(N - 3));
end