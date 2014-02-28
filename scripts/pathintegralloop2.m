%% The path integral loop!
iter = 0;
for c = .00005:.0001:1
iter = iter+1;
F = @(x)(acos(c(iter)/x).*x);
Q(iter) = quad(F,0,1);
end

%%