function [shapeError, volError, timeItTook] = advection(PLICMethod,trial,advectMethod,shapeCase)
% function [shapeError, volError] = advection(PLICMethod,trial,advectMethod,shapeCase)
% PLICMethod: 1-2-3-4 = 'Finite Difference','Center of Mass','ELVIRA','LVIRA'
% trial: Resolution, 1 = 64, 2 = 128 ... 5 = 1024
% advectMethod: 1 = Euler, 2 = Lagrange
% shapeCase: 1 = zalezak, 2 = Deformation Field
if nargin==0
    PLICMethod = 4; %PLICMethod: 1-2-3-4 = FDM, COM, ELVIRA, LVIRA
    trial = 2; % Resolution, 1 = 64, 2 = 128 ... 5 = 1024
    advectMethod = 2; %1 = Euler, 2 = Lagrange
    shapeCase = 2; %1 = zalezak, 2 = Deformation Feild
end
PLMT = {'Finite Difference','Center of Mass','ELVIRA','LVIRA'};
M = 2^(trial+5);
N = M;
dim = [0 1 0 1]; dim2 = [0.25 0.75 0.5 1];

dx = (dim(2)-dim(1))/M;
dy = (dim(4)-dim(3))/N;
k = 0; %iteration count
t = 0;
if shapeCase == 1
    DT = 2*pi;
else
    DT = 8;
end
BC3 = @(p,M,N) p([4 3 3:M+2 M+2 M+1],[4 3 3:N+2 N+2 N+1]);
% Initialization
x = linspace(dim(1)-3*dx/2,dim(2)+3*dx/2,M+4); % cell center locations
y = linspace(dim(3)-3*dy/2,dim(4)+3*dx/2,N+4);
xf = linspace(dim(1),dim(2),M+1); % wall locations
yf = linspace(dim(3),dim(4),N+1);

xgrid = [[xf;xf;nan*xf] [yf;yf;yf]*0+[0;1;nan]]; ygrid = [[xf;xf;xf]*0+[0;1;nan] [yf;yf;nan*yf]];

[X, Y] = meshgrid(x,y); X = X'; Y = Y'; % mesh Grid
[XU, YU] = meshgrid(xf,y(2:end-1)); XU = XU'; YU = YU';
[XV, YV] = meshgrid(x(2:end-1),yf); XV = XV'; YV = YV';

psi = phitopsi(M,N,X,Y,dx,dy,shapeCase); %compute psi from phi
psi2  = psi;
if  shapeCase == 1
    ts = (exp(linspace(-pi/2+asin(.025/.15),3*pi/2-asin(.025/.15),M)*1i)*0.15 + 0.5 + 0.75*1i);
    ts = [0.525+.85*1i ts 0.475+.85*1i 0.525+.85*1i];
else
    ts = exp(2i*pi*(0:1/M:1))*0.15 + 0.5 + 0.75*1i;
end
[mx, my, alpha] = PLIC(M,N,dx,psi,PLICMethod); %compute mx my and alpha

%%%completed initialization
tic;
while t < DT
    dt = min(.25*dx*3/pi,DT-t);
    if dt>(DT/2-t)
        if shapeCase~=1
        end
    end
    if shapeCase == 1
        u = (YU-0.5);
        v = -(XV-0.5);
    else
        u =  -sin(pi*XU).^2.*sin(2*pi*YU)*cos(pi*(t+dt/2)/DT);
        v =  sin(2*pi*XV).*sin(pi*YV).^2*cos(pi*(t+dt/2)/DT);
    end

    if advectMethod == 1
        [psi] = octovof         (M,N,u,v,mx,my,alpha,psi,dy,dt);
    else
        [psi] = lagrange_remap_1(M,N,u,v,mx,my,alpha,psi,dy,dt);
    end
    %up = u; vp = v;

    
    psi(alpha>1.5e19 & alpha<2.5e19)=1; psi(isnan(psi))=0;
    psi = BC3(psi,M,N);
    parea = sum(sum(psi(3:M+2,3:N+2)))*dx*dy;
    t = t + dt;
    k = k + 1;
    %ti(k) = t;
    %pa(k) = parea;
    %fprintf("%i\t%.16f\n",k,parea)
    [mx, my, alpha] = PLIC(M,N,dx,psi,PLICMethod);
    mx = BC3(mx,M,N); my = BC3(my,M,N);
    mx([1 2 M+3 M+4],:) = -mx([1 2 M+3 M+4],:);
    my(:,[1 2 N+3 N+4]) = -my(:,[1 2 N+3 N+4]);

    
    if ((mod(k,100)==1) || (t>=DT)) && 1
        %imagesc(x,y,(psi)')
        [X, Y] = interface(M,N,mx,my,alpha,dx);
        plot(X,Y,real(ts),imag(ts),'--')
        title(sprintf('Area = %e,\tt = %.4f',parea,t));
        axis equal; grid on;
        axis(dim);
        %fprintf('%f\n',ta);
        drawnow;
    end
end
timeItTook = toc;
shapeError = dx*dy*sum(abs(psi(:)-psi2(:)));
volError = dx*dy*abs(sum(psi(:))-sum(psi2(:)));
fprintf('%i%i%i%i\t%i\t%.16e\t%.16e\t%f\n',PLICMethod,trial,advectMethod,shapeCase,M,shapeError,volError,timeItTook)
[X, Y] = interface(M,N,mx,my,alpha,dx);
plot(X,Y,'-','linewidth',2); hold on; plot(real(ts),imag(ts),'--'); hold off
legend('Solution','Exact Solution')
title(sprintf('M = %i, Shape Error = %.4e, Volume Error = %.4e',M,shapeError,volError));
xlabel(PLMT{PLICMethod});
axis equal; grid on;
axis(dim2);
saveas(gcf,[pwd '/images/' sprintf('f%i%i%i%i.fig',...
    shapeCase,advectMethod,PLICMethod,trial)]);
saveas(gcf,[pwd '/images/' sprintf('f%i%i%i%i.emf',...
    shapeCase,advectMethod,PLICMethod,trial)]);