function ind = sysresample(w)
wc = cumsum(w);
N = length(w);
u = ([0:N-1]+rand(1))/N;
ind = zeros(1,N);
k = 1;
for j = 1:N
    while(wc(k) < u(j))
        k = k+1;
    end
    ind(j) = k;
end
end