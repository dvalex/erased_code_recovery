% script to draw timings
timing = load('timing.txt' , '-ascii');
N = timing(:,1);
T = timing(:,2); 
figure;
subplot(1,2,1);
semilogy(N,T, '-kx');
ylim([0.1 10000]); xlim([7,20]);
title('a)');
xlabel('n');ylabel('t (cек.)');
grid on; grid minor;
subplot(1,2,2);
steps = zeros(length(N),1);
for i = 1:length(N)
    steps(i) = nchoosek(nchoosek(N(i),3),3);
end
%T1 = T ./ (N-3).^9;
T1 = T ./ steps;
plot(N, T1, '-kx');
ylim([0,max(T1)*1.5]);xlim([7,20]);
title('b)');
xlabel('n');ylabel('t (cек.)');

grid on; grid minor;