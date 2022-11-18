dt = 0.0001;
N = 100000;
% integral for u(s)^2 between 1 and 10
u = zeros(N,1); t = u;
s = zeros(N,1);
for k = 1:N
    t(k) = k*dt;
    % indicator functions expression
    if t(k)<1
        u(k) = 1;
    elseif t(k)<3
        u(k) = -2;
    elseif t(k)<10
        u(k) = 1/2;
    end
    s(k) =u(k).^2*dt;
end
total1 = 1/2 * sum(s,"all");

% TRAPEZOIDAL RULE for x(t) & \eta(t)
x = zeros(N,1); t = x;
x(1) = 0;
eta = zeros(N,1); t = eta;
eta_der = zeros(N,1);
eta(1) = 0;
eta_der(1) = 2;
for k = 2:N
    t(k) = (k-1)*dt;
    x_tmp = x(k-1);
    x_tmp_old = x_tmp;
    % indicator functions
    if t<1
        x_tmp = x(k-1)+dt/2*(cos(x(k-1))+1+cos(x_tmp)+1);
        x(k) = x_tmp;
    elseif t<3
        x_tmp = x(k-1)+dt/2*(cos(x(k-1))-2+cos(x_tmp)-2);
        x(k) = x_tmp;
    elseif t<10
        x_tmp = x(k-1)+dt/2*(cos(x(k-1))+1/2+cos(x_tmp)+1/2);
        x(k) = x_tmp;
    end
    ss(k) = x(k).^2*dt;
    eta_tmp = eta(k-1);
    eta_tmp_old = eta_tmp;
    eta_tmp = eta(k-1)+dt/2*(cos(eta(k-1))+exp(-t(k-1))+cos(eta_tmp)+exp(-t(k)));
    eta(k) = eta_tmp;
    eta_der(k) = cos(eta(k))+exp(-t(k));
    result(k) = x(k).*eta_der(k)*dt;
end
total2 = 1/2 * sum(ss,"all")
total3 = sum(result,"all")

total = total1+total2-total3


% plot everything
plot(t,x,'r');
hold on;
plot(t,eta,'g');
plot(t,eta_der);
plot(t,u,"k")
hold off;
legend('x(t)','\eta(t)','\eta der(t)','u(t)');
xlabel('t');
ylabel('solution');
