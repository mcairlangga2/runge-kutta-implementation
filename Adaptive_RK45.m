%Solve y'(t)=-2y(t) with y0=3, using Runge Kutta 45
f = @(t,y)(-2*y);       %function to be solved
y0 = 3;                 %Initial Condition
h0=0.1;                 %Time step
hmin = 0.01;            %Time step minimum               
hmax = 0.5;             %Time step maximum
t0 = 0;                 %initial time
tf = 2;                 %final time
%t = t0:h0:tf;
nmax = 20;              %maximum iteration
emin = 1e-06;            %error minimum
emax = 5e-06;            %error maksimum
%yexact = 3*exp(-2*t);   %Exact solution (in general we won't know this
%ystar = zeros(1,nmax);  %Preallocate array (good coding practice)
           
k = 1;
t(1) = t0;
h(1) = h0
w(1) = y0;
w_hat(1) = y0;
err = 0;
br = tf-0.00001*abs(tf);
tic
while t(k)<tf
    if t(k)+h(k)>br
       h(k) = tf - t(k);
    end
    if h(k)<hmin 
        h(k) = hmin;
    elseif h(k)>hmax 
        h(k) = hmax;
    end
    %Computing RK4
    k1 = h(k)*f(t(k),w(k));
    k2 = h(k)*f(t(k)+h(k)/4, w(k)+k1/4);
    k3 = h(k)*f(t(k)+(3/8)*h(k), w(k) + (3/32)*k1 + (9/32)*k2);
    k4 = h(k)*f(t(k)+(12/13)*h(k),w(k) + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    k5 = h(k)*f(t(k)+h(k), w(k) +(439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    k6 = h(k)*f(t(k)+(1/2)*h(k), w(k) -(8/27)*k1 + 2*k2 -(3544/2565)*k3 + (1859/4104)*k4 -(11/40)*k5);
    w(k+1)  = w(k) + (25/216)*k1 +(1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5; %RK4
    w_hat(k+1) = w(k) + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6; %RK5
    err = abs(w(k+1)-w_hat(k+1));
    if (err > emax) && (h(k)>hmin)
        h(k) = h(k)/2;
        continue
    else
        k = k+1;
        h(k) = h(k-1);
        if (t(k-1)+h(k))>br
            t(k) = tf;
        else
        t(k) = t(k-1)+h(k);
        end
    end
    if (err < emin)
        h(k) = 2*h(k-1);
    end
end
toc 

