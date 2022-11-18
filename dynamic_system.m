% first we know that u_1(t) and u_2(t)
% u_1(t) = exp(-t)*cos(t)
% u_2(t) = 
dt = 0.0001;
N = 100000;

% TRAPEZOIDAL RULE
v = zeros(N,1); t = v;
v(1) = 0;
for k = 2:N
    t(k) = (k-1)*dt;
    v_tmp = v(k-1);
    v_tmp_old = v_tmp;
    v_tmp = v(k-1)+dt/2*(cos(v(k-1))+exp(t(k-1))+cos(v_tmp)-exp(t(k)));
    v(k) = v_tmp;
end


% plot everything
plot(t,v,'r')
legend('trapezoidal')
hold off
xlabel('t')
ylabel('solution')
