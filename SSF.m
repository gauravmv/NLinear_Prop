%Split Step Fourier Transform for propogation of Pulses 


N = 512; %Resolution
beta = -25; %GVD
P0 = 6.76; %Peak Input Power
gamma = 1; %Non-Linearity Factor
Tmax = 40; %Temporal grid size
T0 = 10; %Pulse duration
h = Tmax/N; %Temporal step
dx = 0.01; %Spatial Step
x = 0.63; %Fiber Length in m
niter = round(x/dx);
n = (-256:1:255); %Spatial Resolution
t = n*h; %grid
w = 2*pi*n./Tmax; 
%u = sqrt(P0)*exp(1j*0).*sech(t/T0); %initial pulse shape - secant 
u = sqrt(P0)*exp(1*j*0)*exp(-t.^2/(2*T0*T0)); %Initial pulse shape - gaussian
%u = sqrt(P0)*(1 + 0.01*cos(2.7*10e12*t/T0)); %Initial Pulse shape -
%Constant power with modulation
op_pulse = zeros(niter, 512);

hold on

plot(t, abs(u))
xlabel('Time (ps)'); ylabel('Amplitude (W^{1/2})')
axis tight manual 
set(gca,'nextplot','replacechildren');
axis([-20 20 0 9])
hold on;


v = VideoWriter('P=Dynamics.avi');
open(v);

for i=1:niter
    u = exp(gamma*dx*1j*abs(u).*abs(u)).*u;
    c = fftshift(fft(u));
    c = exp(beta*dx*1j*w.*w/2).*c;
    u = ifft(fftshift(c));
    op_pulse(i,:) = u;
    p1 = plot(t, abs(op_pulse(i,:)));
    title(['Pulse at z =' num2str(10*i) 'm'])
    frame = getframe(gcf);
    writeVideo(v,frame);
    delete(p1);
    
end


plot(t, abs(op_pulse(niter,:)))
frame = getframe(gcf);
writeVideo(v, frame)
close(v);


