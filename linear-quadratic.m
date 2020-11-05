clear
N=100000;
K=100;
T=2;
dt=0.2;
M=T/dt;
L=10;
R=0.01;
nu=1;
lambda=nu*R;
psi=zeros(K+1, M+1, "gpuArray");
vec_m=[0];
stdt = sqrt(dt);
%xx=zeros(M+1,N,"gpuArray");
for m=vec_m
    for k=0:K  
        x0 = 2*L*(k - K/2)/K; 
        x = x0*ones(N,1,"gpuArray");
        %xx(1,:)=x;
        for j=m:M-1            
            dxi = stdt*randn(N,1,"gpuArray");           
            x = x + dxi;
            %xx(j+2,:)=x;
        end
        psi(k+1,m+1)=sum(exp(-0.5*(x.^2)/lambda))/N;
    end
end
%xlim([0 2]);
%ylim([-10 10]);
%plot([0:dt:T], xx);
%return
for m=vec_m
    plot(linspace(-L,L,K+1),-lambda*log(psi(:,m+1)));
    ylim([0 0.5]);
    hold on;
    
    x=-10:10;
    t = m*dt;
    sigma=sqrt(nu*(T-t));
    sigma1=1/sqrt((1/sigma^2+1/(nu*R)));
    if (m < M)
        J=nu*R*log(sigma/sigma1)+1/2*(sigma1^2/sigma^2)*x.^2;
    else 
        J=1/2*x.^2;
    end
    plot((-10:10), J,'r+');
    hold on;
end
hold on;
    
hold off;
