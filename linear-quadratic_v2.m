clear
N=100000;
JK=20;
L=1;
JM=10;
T=10;
dt=0.1;
M=T/dt;
R=20;
nu=0.05;
lambda=nu*R;
type="gpuArray";

J = zeros(JK, JM, type);
psi = zeros(JK, JM, type);
u = zeros(JK, JM, type);
for mi=1:JM
    m = (mi-1)/(JM-1)*(M-1);
    for ki=1:JK
        x0 = (ki-1)/JK*L; 
        x = x0*ones(N, 1, type);
        for j=m:M-1            
            dxi = sqrt(dt*nu)*randn(N, 1, type);     
            if j == m; dxim = dxi; end
            V = 0.5*x.^2;
            prob = V*dt/lambda;
            idx =  prob < rand(N, 1, type);
            x(~idx) = NaN;            
            x(idx) = x(idx) + x(idx)*dt + dxi(idx);
        end
        psi(ki,mi) = sum(~isnan(x))/N;
        u(ki,mi) = 1/psi(ki,mi)*sum(~isnan(x).*dxim)/N/dt;
    end
    J(:,mi) = -lambda*log(psi(:, mi));    
end
surf(T*(0:JM-1)/JM, L*(0:JK-1)/JK, u);
