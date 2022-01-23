% Solve y'(t)=-2y(t) with y0=3, 4th order Runge Kutta, 0<=t<=2
f = @(t,y)(-2*y);
t0 = 0;
tf = 2;
y0 = 3;
h  = 0.1;   %Step Size
t_size  = t0:h:tf;
t = t0;    
w(1)  = y0;    %Initial condition
tic
for i = 1:(length(t_size)-1)
   k1 = f(t,w(i));
   k2 = f(t+h/4, w(i)+k1*h/4.0);
   k3 = f(t+h/4, w(i)+k1*h/8+k2*h/8);
   k4 = f(t+h/2, w(i)-k2*h/2+k3*h);
   k5 = f(t+3*h/4, w(i)+3*k1*h/16+9*k4*h/16);
   k6 = f(t+h, w(i) - 3*k1*h/7 + 2*k2*h/7 + 12*k3*h/7 - 12*k4*h/7 + 8*k5*h/7);
   
   w(i+1) = w(i)+(1/90)*(7*k1+32*k3+12*k4+32*k5+7*k6)*h;
   t = t0+i*h
   %plot(t,w,'b*');
   %hold on;
end
toc


