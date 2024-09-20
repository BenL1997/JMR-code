% Input parameters
global t1b t2b t3b wb w1b h1b h2b l1b l2b l3b nb t1h t2h t3h wh w1h h1h h2h l1h l2h l3h nh t1k t2k t3k wk w1k h1k h2k l1k l2k c10 c20 p h3 R0 d F h lb lh;
t1b=1.5;t2b=3; t3b=1;wb=17; w1b=15; h1b=13; h2b=4; l1b=2; l2b=2; nb=11;
t1k=0.5; t2k=3; t3k=1;wk=17; w1k=15; h1k=14; h2k=3; l1k=0.7; l2k=0.9; 
t1h=1.5; t2h=3; t3h=1;wh=17; w1h=15; h1h=13; h2h=3; l1h=2.3; l2h=1.8; nh=9;
c10=0.116;c20=0.025; p=0.04; h=17.5;

% Input coordinates
coordinates=[-71.44  -72.78  -71.55  -68.17  -63.08  -56.72  -49.50  -41.78  -33.88  -26.08  -18.62  -11.69  -5.44  0.00  0	-1.58	-6.03	-12.41	-19.42	-25.40	-29.04	-29.30	-26.09	-20.22	-13.20	-7.17
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-0.12 -0.95	-3.05	-6.78	-12.08	-18.67	-25.70	-32.34	-37.76	-41.43	-43.24
-20.15  -12.40  -4.62  2.51  8.53  13.21  16.40  18.06  18.22  16.95  14.36  10.57  5.73  0.00  -15.5	-23.52	-30.28	-34.78	-36.50	-35.44	-32.26	-28.25	-24.91	-23.68	-25.53	-30.51];

% Divide the input coordinates into three parts: bending, buckling, and helical
x=coordinates(1,:);
y=coordinates(2,:);
z=coordinates(3,:);
bending=[1,14];straight=[14,15];buckling=[14];helical=[15,26];% Define the range of each segment
l_b=bending(2)-bending(1)+1;
l_h=helical(2)-helical(1)+1;
n_k=size(straight, 1);

% Define arrays, where x_b and y_b are the coordinates of the bending segment, x_h, y_h and z_h are the coordinates of the helical segment, and arrays g_b, g_k, and g_h store the parameters for the bending, buckling, and helical segments respectively
x_b=x(bending(1):bending(2));
y_b=z(bending(1):bending(2));
x_h=x(helical(1):helical(2));
y_h=y(helical(1):helical(2));
z_h=z(helical(1):helical(2));
n_1=14;% Define the number of feature points in the bending segment
g_b=zeros(n_1,2);% Store the parameters of the bending segment. In this paper, the iterated parameter we choose is w1b, so only w1b and t3b change, while other parameters remain unchanged. The first column is w1b, and the second column is t3b
g_k=zeros(9,1);
n_2=12;% Define the number of feature points in the helical segment
g_h=zeros(n_2,3);% Store the parameters of the helical segment. Noted that the parameters correspond to those of the straight pneu-net (the first column is w1b, the second column is t3b) and the helical direction angle (the third column)

% Inverse design of the bending segment
% Fit the coordinate points
s=zeros(l_b,1);
s(1)=0;
for i=2:l_b
        s(i)=s(i-1)+(((x_b(i)-x_b(i-1)))^2+((y_b(i)-y_b(i-1)))^2)^0.5;
end
s_1=zeros(100*(n_1-1),1);
l_b1=zeros(n_1-1,1);
ss = linspace(0, max(s), 100*(n_1-1));
xx = spline(s, x_b', ss);% Cubic spline interpolation fitting
yy = spline(s, y_b', ss);
j=1;
for i=2:100*(n_1-1)
        s_1(i)=s_1(i-1)+(((xx(i)-xx(i-1)))^2+(yy(i)-yy(i-1))^2)^0.5;
        qy=mod(i,100);
        if(qy==0)
            if(i~=100)
            l_b1(j)= s_1(i)-s_1(i-100);
            j=j+1;
            else 
            l_b1(j)= s_1(i)-s_1(i-99);
            j=j+1;
            end
        end
end
lb=s_1(end); l3b=lb/nb-2*l2b-l1b;

% Solve for curvature
j=2;
s_e=zeros(n_1,1);
for i=100:100:100*(n_1-2)
    rr=ceshi(xx(i-1),yy(i-1),xx(i),yy(i),xx(i+1),yy(i+1));% Calculate curvature using the three-point circle method
    s_e(j)=rr;
    j=j+1;
end
rr=ceshi(xx(1),yy(1),xx(2),yy(2),xx(3),yy(3));
s_e(1)=rr;
rr=ceshi(xx(100*(n_1-1)-2),yy(100*(n_1-1)-2),xx(100*(n_1-1)-1),yy(100*(n_1-1)-1),xx(100*(n_1-1)),yy(100*(n_1-1)));
s_e(n_1)=rr;

% Solve for cross-sectional parameters
for i=1:n_1
b0 = 0.5;
b = fsolve(@(b)fvn(b),b0); 
r=lb/b-0.5*t2b;
if(r<s_e(i))
    while(r<s_e(i))
    w1b=w1b-0.02;
    t3b=t3b+0.01;
    b0 = 1;
    b = fsolve(@(b)fvn(b),b0); 
    r=lb/b-0.5*t2b;
    end
else
    while(r>s_e(i))
    w1b=w1b+0.02;
    t3b=t3b-0.01;
    b0 = 1;
    b = fsolve(@(b)fvn(b),b0); 
    r=lb/b-0.5*t2b;
    end
end
    g_b(i,1)=w1b;g_b(i,2)=t3b;
end

% Calculate the inflection point where the deformation direction changes
dt1=zeros(100*n_1,1);
j=1;
for i=4:100*(n_1-1)
A = [xx(i-2)-xx(i-3), yy(i-2)-yy(i-3)];
B = [xx(i-1)-xx(i-2), yy(i-1)-yy(i-2)];
C = [xx(i)-xx(i-1), yy(i)-yy(i-1)];
cross_product= (A(1) * B(2) - A(2) * B(1))*(B(1) * C(2) - B(2) * C(1));
if(cross_product<0)
    dt1(j)=i;
    j=j+1;
end
end

% Determine whether there is sidewall contact; if pc is less than the working pressure, sidewall contact may occur. This can be resolved by increasing l2b or decreasing l3b in the corresponding chamber. Since l2b and l3b have minimal impact on bending, the curvature at this location can be considered unchanged
h3=h1b-h2b;
pc=zeros(n_1,1);
for i=1:n_1
w1b=g_b(i,1);t3b=g_b(i,2);
q= max(w1b,h3)/ min(w1b,h3);
d=(p* (0.5*max(w1b,h3))^4/2/c10/l2b)^(1/3)* (0.183*q*q-0.9555*q+1.5905);
    c0 = [1,3,9,0.584,11.17,0.03];
    c = fsolve(@(c)ffn(c),c0); 
pc(i)=c(6);
end


% Inverse design of the buckling segment
% Calculate the buckling angle
u = [xx((n_1-1)*100)-xx((n_1-1)*100-1), yy((n_1-1)*100)-yy((n_1-1)*100-1)]; 
v = [x(straight(2))-x(straight(1)), z(straight(2))-z(straight(1))];  
cosTheta = dot(u, v) / (norm(u) * norm(v));  
theta = acos(cosTheta); 

% Solve for geometric parameters
a=[0.259	0.259	0.259	0.343	0.343	0.343	0.343	0.408	0.408	0.408	0.408	0.408	0.466	0.466	0.466	0.466	0.466	0.522	0.522	0.522	0.522	0.522	0.522	0.573	0.573	0.573	0.573	0.573	0.573	0.624	0.624	0.624	0.624	0.624	0.624	0.624	0.675	0.675	0.675	0.675	0.675	0.675	0.675	0.72	0.72	0.72	0.72	0.72	0.72	0.72	0.72	0.766	0.766	0.766	0.766	0.766	0.766	0.766	0.766	0.816	0.816	0.816	0.816	0.816	0.816	0.816	0.816	0.816	0.863	0.863	0.863	0.863	0.863	0.863	0.863	0.863	0.863	0.91	0.91	0.91	0.91	0.91	0.91	0.91	0.91	0.91	0.91	0.953	0.953	0.953	0.953	0.953	0.953	0.953	0.953	0.953	0.953
0.1	0.2	0.259	0.1	0.2	0.3	0.343	0.1	0.2	0.3	0.4	0.408	0.1	0.2	0.3	0.4	0.466	0.1	0.2	0.3	0.4	0.5	0.522	0.1	0.2	0.3	0.4	0.5	0.573	0.1	0.2	0.3	0.4	0.5	0.6	0.624	0.1	0.2	0.3	0.4	0.5	0.6	0.675	0.1	0.2	0.3	0.4	0.5	0.6	0.7	0.72	0.1	0.2	0.3	0.4	0.5	0.6	0.7	0.766	0.1	0.2	0.3	0.4	0.5	0.6	0.7	0.8	0.816	0.1	0.2	0.3	0.4	0.5	0.6	0.7	0.8	0.863	0.1	0.2	0.3	0.4	0.5	0.6	0.7	0.8	0.9	0.91	0.1	0.2	0.3	0.4	0.5	0.6	0.7	0.8	0.9	0.953
1.184	0.659	0	1.338	1.02	0.545	0	1.438	1.181	0.841	0.302	0	1.511	1.287	0.993	0.586	0	1.538	1.344	1.074	0.721	0.343	0	1.563	1.386	1.201	0.882	0.601	0	1.599	1.428	1.257	0.985	0.738	0.326	0	1.622	1.456	1.3	1.106	0.819	0.526	0	1.65	1.503	1.372	1.131	0.941	0.691	0.338	0	1.695	1.518	1.399	1.203	1.011	0.794	0.448	0	1.726	1.578	1.423	1.217	1.071	0.848	0.592	0.288	0	1.732	1.61	1.456	1.305	1.109	0.941	0.731	0.472	0	1.736	1.627	1.483	1.356	1.168	1.017	0.832	0.572	0.186	0	1.74	1.639	1.512	1.409	1.211	1.076	0.911	0.679	0.433	0];
F = scatteredInterpolant(a(1,:)', a(2,:)', a(3,:)', 'natural');
for i=1:n_k
    h3=h1k-h2k;
q= max(w1k,h3)/ min(w1k,h3);
d=(p* (0.5*max(w1k,h3))^4/2/c10/l2k)^(1/3)* (0.183*q*q-0.9555*q+1.5905);
R0=(0.25*h3)*((0.5*h3)/d+d/(0.5*h3));
x0 = [0.5,1,4,8,3,3,1.3];
k = fsolve(@(x)fvzn(x),x0); 
buk_angle=k(1);
while(buk_angle>theta)
    w1k=w1k-0.1;
    q= max(w1k,h3)/ min(w1k,h3);
    d=(p* (0.5*max(w1k,h3))^4/2/c10/l2k)^(1/3)* (0.183*q*q-0.9555*q+1.5905);
    R0=(0.25*h3)*((0.5*h3)/d+d/(0.5*h3));
    x0 = [0.5,1,4,8,3,3,1.3];
    k = fsolve(@(x)fvzn(x),x0); 
    buk_angle=k(1);
    if(buk_angle<theta)
    break;
    end
    h1k=h1k-0.1;
    h3=h1k-h2k;
    q= max(w1k,h3)/ min(w1k,h3);
    d=(p* (0.5*max(w1k,h3))^4/2/c10/l2k)^(1/3)* (0.183*q*q-0.9555*q+1.5905);
    R0=(0.25*h3)*((0.5*h3)/d+d/(0.5*h3));
    x0 = [0.5,1,4,8,3,3,1.3];
    k = fsolve(@(x)fvzn(x),x0); 
    buk_angle=k(1);
    if(buk_angle<theta)
    break;
    end
    h2k=h2k+0.1;
    h3=h1k-h2k;
    q= max(w1k,h3)/ min(w1k,h3);
    d=(p* (0.5*max(w1k,h3))^4/2/c10/l2k)^(1/3)* (0.183*q*q-0.9555*q+1.5905);
    R0=(0.25*h3)*((0.5*h3)/d+d/(0.5*h3));
    x0 = [0.5,1,4,8,3,3,1.3];
    k = fsolve(@(x)fvzn(x),x0); 
    buk_angle=k(1);
if(buk_angle<theta)
    break;
end
    l1k=l1k+0.05;
    x0 = [0.5,1,4,8,3,3,1.3];
    k = fsolve(@(x)fvzn(x),x0); 
    buk_angle=k(1);
if(buk_angle<theta)
    break;
end
    l2k=l2k+0.05;
    d=(p* (0.5*max(w1k,h3))^4/2/c10/l2k)^(1/3)* (0.183*q*q-0.9555*q+1.5905);
    R0=(0.25*h3)*((0.5*h3)/d+d/(0.5*h3));
    x0 = [0.5,1,4,8,3,3,1.3];
    k = fsolve(@(x)fvzn(x),x0); 
    buk_angle=k(1);
if(buk_angle<theta)
    break;
end
end
end
t1k=h-t2k-h1k;t3k=(wk-w1k)*0.5;
g_k(1)=t1k; g_k(2)=t2k; g_k(3)=t3k; g_k(4)=wk; g_k(5)=w1k; g_k(6)=h1k; g_k(7)=h2k; g_k(8)=l1k; g_k(9)=l2k;


%Inverse design of the helical segment
% Fit the coordinate points
s=zeros(length(x_h),1);
for i=2:length(x_h)
        s(i)=s(i-1)+(((x_h(i)-x_h(i-1)))^2+((y_h(i)-y_h(i-1)))^2+(z_h(i)-z_h(i-1))^2)^0.5;
end
lh=max(s);
l3h=lh/nh-2*l2h-l1h;
ss = linspace(0, max(s), 10*((n_2-1)+1));
xx = spline(s, x_h', ss);
yy = spline(s, y_h', ss);
zz = spline(s, z_h', ss);

% Adjust the curve's pose
X = xx-xx(1);
Y = yy-yy(1);
Z = zz-zz(1);
X1 = x_h-x_h(1);
Y1 = y_h-y_h(1);
Z1 = z_h-z_h(1);
points=[X;Y;Z];
points1=[X1;Y1;Z1];
if (Y(2)>0)
    c=-acos(X(2)/(X(2)^2+Y(2)^2)^0.5);
else
    c=acos(X(2)/(X(2)^2+Y(2)^2)^0.5);
end
    b=pi-acos(Z(2)/(X(2)^2+Y(2)^2+Z(2)^2)^0.5);

Ry = [cos(b), 0, sin(b);
          0, 1, 0;
          -sin(b), 0, cos(b)];
Rz = [cos(c), -sin(c), 0;
          sin(c), cos(c), 0;
          0, 0, 1];
new = Ry*Rz * points;
newh = Ry*Rz * points1;
c=-0.5*pi;
Rz = [cos(c), -sin(c), 0;
          sin(c), cos(c), 0;
          0, 0, 1];
new = Rz * new;
newh = Rz * newh;
temp = new(1, :);  
new(1, :) = -new(2, :);  
new(2, :) = temp;  
temp = newh(1, :);  
newh(1, :) = -newh(2, :);  
newh(2, :) = temp;  

% Calculate the direction of the helical axis
S=zeros(3,1);
mm=99999999999;
for i=0:1:359 
    Rz = [cos(i/180*pi), -sin(i/180*pi), 0;
          sin(i/180*pi), cos(i/180*pi), 0;
          0, 0, 1];
for j=0:-1:-89
Rx = [1, 0, 0;
          0, cos(j/180*pi), -sin(j/180*pi);
          0, sin(j/180*pi), cos(j/180*pi)];
new1 = Rx *Rz *new;
x_h1=new1(1,:);z_h1=new1(3,:);
s2=zeros(length(x_h1),1);
l_h1=zeros(n_2-1,1);
m=1;
for k=2:length(x_h1)
        s2(k)=s2(k-1)+(((x_h1(k)-x_h1(k-1)))^2+(z_h1(k)-z_h1(k-1))^2)^0.5;
        qy=mod(k,10);
        if(qy==0)
            if(k~=10)
            l_h1(m)= s2(k)-s2(k-10);
            m=m+1;
            else 
            l_h1(m)= s2(k)-s2(k-9);
            m=m+1;
            end
        end
end
aver = mean(l_h1);
S1=0;
for ii=1:n_2-1
    S1=S1+(l_h1(ii)-aver)^2;
end
if(S1<mm)
    mm=S1;
S(1)=S1;S(2)=i;S(3)=j;
end
end
end
Rz = [cos(S(2)/180*pi), -sin(S(2)/180*pi), 0;
          sin(S(2)/180*pi), cos(S(2)/180*pi), 0;
          0, 0, 1];
Rx = [1, 0, 0;
          0, cos(S(3)/180*pi), -sin(S(3)/180*pi);
          0, sin(S(3)/180*pi), cos(S(3)/180*pi)];
new1 = Rx *Rz *newh;
new2 =Rz *newh;
x_h1=new1(1,:);z_h1=new1(3,:);
s=zeros(l_h,1);
s(1)=0;
for i=2:l_h
        s(i)=s(i-1)+(((x_h1(i)-x_h1(i-1)))^2+((z_h1(i)-z_h1(i-1)))^2)^0.5;
end
s_2=zeros(100*(n_2-1),1);
l_h1=zeros(n_2-1,1);
ss = linspace(0, max(s), 100*(n_2-1)); 
xx = spline(s, x_h1', ss);
yy = spline(s, z_h1', ss);
j=1;
for i=2:100*(n_2-1)
        s_2(i)=s_2(i-1)+(((xx(i)-xx(i-1)))^2+(yy(i)-yy(i-1))^2)^0.5;
        qy=mod(i,100);
        if(qy==0)
            if(i~=100)
            l_h1(j)= s_2(i)-s_2(i-100);
            j=j+1;
            else 
            l_h1(j)= s_2(i)-s_2(i-99);
            j=j+1;
            end
        end
end

% Calculate the curvature
j=2;
s_e=zeros(n_2,1);
for i=100:100:100*(n_2-2)
    rr=ceshi(xx(i-1),yy(i-1),xx(i),yy(i),xx(i+1),yy(i+1));
    s_e(j)=rr;
    j=j+1;
end
rr=ceshi(xx(1),yy(1),xx(2),yy(2),xx(3),yy(3));
s_e(1)=rr;
rr=ceshi(xx(100*(n_2-1)-2),yy(100*(n_2-1)-2),xx(100*(n_2-1)-1),yy(100*(n_2-1)-1),xx(100*(n_2-1)),yy(100*(n_2-1)));
s_e(n_2)=rr;

% Solve for cross-sectional parameters
t1b=t1h;t2b=t2h; t3b=t3h/cos(-S(3)/180*pi);wb=wh/cos(-S(3)/180*pi); w1b=w1h/cos(-S(3)/180*pi); h1b=h1h; h2b=h2h; l1b=l1h*cos(-S(3)/180*pi); l2b=l2h*cos(-S(3)/180*pi); nb=nh;lb=lh*cos(-S(3)/180*pi);l3b=l3h*cos(-S(3)/180*pi);
for i=1:n_2
b0 = 0.2;
b = fsolve(@(b)fvn(b),b0); 
r=lb/b-0.5*t2b;
if(r<s_e(i))
    while(r<s_e(i))
    w1b=w1b-0.02;
    t3b=t3b+0.01;
    b0 = 0.2;
    b = fsolve(@(b)fvn(b),b0); 
    r=lb/b-0.5*t2b;
    end
else
    while(r>s_e(i)) 
    w1b=w1b+0.02;
    t3b=t3b-0.01;
    b0 = 0.2;
    b = fsolve(@(b)fvn(b),b0); 
    r=lb/b-0.5*t2b;
    end
end
    g_h(i,1)=w1b;g_h(i,2)=t3b;
end
g_h(:,3)=-S(3);

% Calculate the inflection point where the deformation direction changes
dt2=zeros(100*(n_2-1),1);
j=1;
for i=4:100*(n_2-1)
A = [xx(i-2)-xx(i-3), yy(i-2)-yy(i-3)];
B = [xx(i-1)-xx(i-2), yy(i-1)-yy(i-2)];
C = [xx(i)-xx(i-1), yy(i)-yy(i-1)];
cross_product= (A(1) * B(2) - A(2) * B(1))*(B(1) * C(2) - B(2) * C(1));
if(cross_product<0)
    dt2(j)=i;
    j=j+1;
end
end


% Functions
% Curvature solving function
function [rr,xo,yo]=ceshi(x1,y1,x2,y2,x3,y3)
A=x1*(y2-y3)-y1*(x2-x3)+x2*y3-x3*y2;
B=(x1*x1+y1*y1)*(y3-y2)+(x2*x2+y2*y2)*(y1-y3)+(x3*x3+y3*y3)*(y2-y1);
C=(x1*x1+y1*y1)*(x2-x3)+(x2*x2+y2*y2)*(x3-x1)+(x3*x3+y3*y3)*(x1-x2);
D=(x1*x1+y1*y1)*(x3*y2-x2*y3)+(x2*x2+y2*y2)*(x1*y3-x3*y1)+(x3*x3+y3*y3)*(x2*y1-x1*y2);
yo=-C/(2*A);
xo=-B/(2*A);
rr=sqrt((B*B+C*C-4*A*D)/(4*A*A));
end

% Bending angle solving function
function b = fvn(b) 
global lb l3b l1b h1b h2b t3b t2b l2b c10 nb w1b wb c20 h p t1b;
b(1)=2*((((0.5*t2b))*lb*wb/nb*((((lb-(0.25*t2b)*b(1))/lb)^2+1/((lb-(0.25*t2b)*b(1))/lb)^2-2)*c10+c20*(((lb-(0.25*t2b)*b(1))/lb)^2+1/((lb-(0.25*t2b)*b(1))/lb)^2-2)^2))+(((0.5*t2b))*lb*wb/nb*((((lb+(0.25*t2b)*b(1))/lb)^2+1/((lb+(0.25*t2b)*b(1))/lb)^2-2)*c10+c20*(((lb+(0.25*t2b)*b(1))/lb)^2+1/((lb+(0.25*t2b)*b(1))/lb)^2-2)^2))+(((((lb+0.5*(h2b+t2b)*b(1))/(lb+0.5*(h2b+0.5*t2b)*b(1)))^2+1/((lb+0.5*(h2b+t2b)*b(1))/(lb+0.5*(h2b+0.5*t2b)*b(1)))^2-2)*c10+c20*(((lb+0.5*(h2b+t2b)*b(1))/(lb+0.5*(h2b+0.5*t2b)*b(1)))^2+1/((lb+0.5*(h2b+t2b)*b(1))/(lb+0.5*(h2b+0.5*t2b)*b(1)))^2-2)^2)*((h2b+t2b)-t2b)*(l1b+2*l2b)*w1b)+((l3b+2*l2b)*((h-(h2b+t2b))-t1b)*2*(2.2+(t3b-2.5)*0.4)*(((((((lb+0.5*(h2b+0.5*t2b)*b(1))/nb-l1b)/(lb+0.5*(h2b+0.5*t2b)*b(1))*(lb+(h2b+0.5*t2b)*b(1)))+l3b+2*l2b)/2/(l3b+2*l2b))^2+1/(((((lb+0.5*(h2b+0.5*t2b)*b(1))/nb-l1b)/(lb+0.5*(h2b+0.5*t2b)*b(1))*(lb+(h2b+0.5*t2b)*b(1)))+l3b+2*l2b)/2/(l3b+2*l2b))^2-2)*c10+c20*((((((lb+0.5*(h2b+0.5*t2b)*b(1))/nb-l1b)/(lb+0.5*(h2b+0.5*t2b)*b(1))*(lb+(h2b+0.5*t2b)*b(1)))+l3b+2*l2b)/2/(l3b+2*l2b))^2+1/(((((lb+0.5*(h2b+0.5*t2b)*b(1))/nb-l1b)/(lb+0.5*(h2b+0.5*t2b)*b(1))*(lb+(h2b+0.5*t2b)*b(1)))+l3b+2*l2b)/2/(l3b+2*l2b))^2-2)^2))+ (((((lb+0.5*(h2b+t2b)*b(1))/lb)^2+1/((lb+0.5*(h2b+t2b)*b(1))/lb)^2-2)*c10+c20*(((lb+0.5*(h2b+t2b)*b(1))/lb)^2+1/((lb+0.5*(h2b+t2b)*b(1))/lb)^2-2)^2)*lb*((h2b+t2b)-t2b)*2*(2.2+(t3b-2.5)*0.4)/nb))/(((h2b+t2b)-t2b)*w1b*0.5*((h2b+t2b)-t2b)*(b(1)/nb-(l1b+2*l2b)*b(1)/(lb+0.5*(h2b+0.5*t2b)*b(1)))+((h2b+t2b)-t2b)*w1b*((((lb/b(1))+((0.5*t2b)))*tan(b(1)/nb-(l1b+2*l2b)*b(1)/(lb+0.5*(h2b+0.5*t2b)*b(1)))-l3b)*cos(b(1)/nb-(l1b+2*l2b)*b(1)/(lb+0.5*(h2b+0.5*t2b)*b(1))))+(h1b+t2b-(h2b+t2b))^2*w1b*asin(((((lb+0.5*(h2b+0.5*t2b)*b(1))/nb-(l1b+2*l2b))/(lb+0.5*(h2b+0.5*t2b)*b(1))*(lb+(h2b+0.5*t2b)*b(1)))-l3b)/2/(h1b+t2b-(h2b+t2b))))-p;
end

% Critical lateral wall contact pressure solving function
function Q = ffn(c)
global lb l3b l1b h1b h2b t3b l2b t1b t2b c10 nb w1b wb c20 h d ;
Q(1) = c(6)- 2*((((0.5*t2b))*lb*wb/nb*((((lb-(0.5*(0.5*t2b))*c(1))/lb)^2+1/((lb-(0.5*(0.5*t2b))*c(1))/lb)^2-2)*c10+c20*(((lb-(0.5*(0.5*t2b))*c(1))/lb)^2+1/((lb-(0.5*(0.5*t2b))*c(1))/lb)^2-2)^2))+(((0.5*t2b))*lb*wb/nb*((((lb+(0.5*(0.5*t2b))*c(1))/lb)^2+1/((lb+(0.5*(0.5*t2b))*c(1))/lb)^2-2)*c10+c20*(((lb+(0.5*(0.5*t2b))*c(1))/lb)^2+1/((lb+(0.5*(0.5*t2b))*c(1))/lb)^2-2)^2))+(((((lb+0.5*(h2b+t2b)*c(1))/(lb+0.5*((h2b+t2b)-((0.5*t2b)))*c(1)))^2+1/((lb+0.5*(h2b+t2b)*c(1))/(lb+0.5*((h2b+t2b)-((0.5*t2b)))*c(1)))^2-2)*c10+c20*(((lb+0.5*(h2b+t2b)*c(1))/(lb+0.5*((h2b+t2b)-((0.5*t2b)))*c(1)))^2+1/((lb+0.5*(h2b+t2b)*c(1))/(lb+0.5*((h2b+t2b)-((0.5*t2b)))*c(1)))^2-2)^2)*((h2b+t2b)-(2*(0.5*t2b)))*(l1b+2*l2b)*w1b)+((l3b+2*l2b)*((h-(h2b+t2b))-t1b)*2*(2.2+(t3b-2.5)*0.4)*(((((((lb+0.5*((h2b+t2b)-((0.5*t2b)))*c(1))/nb-l1b)/(lb+0.5*((h2b+t2b)-((0.5*t2b)))*c(1))*(lb+((h2b+t2b)-((0.5*t2b)))*c(1)))+l3b+2*l2b)/2/(l3b+2*l2b))^2+1/(((((lb+0.5*((h2b+t2b)-((0.5*t2b)))*c(1))/nb-l1b)/(lb+0.5*((h2b+t2b)-((0.5*t2b)))*c(1))*(lb+((h2b+t2b)-((0.5*t2b)))*c(1)))+l3b+2*l2b)/2/(l3b+2*l2b))^2-2)*c10+c20*((((((lb+0.5*((h2b+t2b)-((0.5*t2b)))*c(1))/nb-l1b)/(lb+0.5*((h2b+t2b)-((0.5*t2b)))*c(1))*(lb+((h2b+t2b)-((0.5*t2b)))*c(1)))+l3b+2*l2b)/2/(l3b+2*l2b))^2+1/(((((lb+0.5*((h2b+t2b)-((0.5*t2b)))*c(1))/nb-l1b)/(lb+0.5*((h2b+t2b)-((0.5*t2b)))*c(1))*(lb+((h2b+t2b)-((0.5*t2b)))*c(1)))+l3b+2*l2b)/2/(l3b+2*l2b))^2-2)^2))+ (((((lb+0.5*(h2b+t2b)*c(1))/lb)^2+1/((lb+0.5*(h2b+t2b)*c(1))/lb)^2-2)*c10+c20*(((lb+0.5*(h2b+t2b)*c(1))/lb)^2+1/((lb+0.5*(h2b+t2b)*c(1))/lb)^2-2)^2)*lb*((h2b+t2b)-(2*(0.5*t2b)))*2*(2.2+(t3b-2.5)*0.4)/nb))/(((h2b+t2b)-(2*(0.5*t2b)))*w1b*0.5*((h2b+t2b)-(2*(0.5*t2b)))*(c(1)/nb-(l1b+2*l2b)*c(1)/(lb+0.5*((h2b+t2b)-((0.5*t2b)))*c(1)))+((h2b+t2b)-(2*(0.5*t2b)))*w1b*((((lb/c(1))+((0.5*t2b)))*tan(c(1)/nb-(l1b+2*l2b)*c(1)/(lb+0.5*((h2b+t2b)-((0.5*t2b)))*c(1)))-l3b)*cos(c(1)/nb-(l1b+2*l2b)*c(1)/(lb+0.5*((h2b+t2b)-((0.5*t2b)))*c(1))))+(h1b+(2*(0.5*t2b))-(h2b+t2b))^2*w1b*asin(((((lb+0.5*((h2b+t2b)-((0.5*t2b)))*c(1))/nb-(l1b+2*l2b))/(lb+0.5*((h2b+t2b)-((0.5*t2b)))*c(1))*(lb+((h2b+t2b)-((0.5*t2b)))*c(1)))-l3b)/2/(h1b+(2*(0.5*t2b))-(h2b+t2b))));
Q(2) = 2*(lb/c(1)+(h2b+t2b)-(0.5*t2b))*sin((l1b)/2/(lb/c(1)+0.5*((h2b+t2b)-(0.5*t2b))))-(c(2)-18.248148*c(6));
Q(3) = 2*(((h2b+t2b)-(0.5*t2b))+lb/c(1))*sin(0.5*(c(1)/nb-2*asin(0.5*c(2)/(((h2b+t2b)-(0.5*t2b))+lb/c(1)))))-c(3);
Q(4) = acos(0.5*(c(3)-l3b-2*l2b)/((h-(h2b+t2b))-t1b))-c(1)/2/nb-acos(((h-(h2b+t2b))-t1b)/2/c(5))-c(4);
Q(5) = c(5)/(0.5*(h-(h2b+t2b))-0.5*t1b)-0.5*((0.5*(h-(h2b+t2b))-0.5*t1b)/d+d/(0.5*(h-(h2b+t2b))-0.5*t1b));
Q(6) = 2*c(5)-2*c(5)*cos(c(4))-c(2);
end

% Buckling angle calculation function
function Q = fvzn(x)
global l2k w1k h3 h2k l1k t2k c10 p d R0 F;
l=0.25*min(w1k,h3)*F(d/0.5/min(w1k,h3),x(5)/0.5/min(w1k,h3));
Q(1) = 0.5*pi-0.5*x(1)-acos((h3)/2/R0)-x(2);
Q(2) = R0*sin(x(2))-x(3);
Q(3) = x(3)+ (h2k+0.5*t2k)-(l1k/x(1)+ (h2k+0.5*t2k))*(1-cos(0.5*x(1)))-x(4);
Q(4) = (R0^2-x(3)^2)^0.5+(l1k/x(1)+ (h2k+0.5*t2k))*sin(0.5*x(1))-(R0^2-0.25*(h3)^2)^0.5-x(5);
Q(5) = pi*l^2*(max(w1k,h3)/ min(w1k,h3))*p*x(4)-x(6);
Q(6)= (l1k+2*l2k+0.5*(h2k+0.5*t2k)*x(1))/(l1k+2*l2k)-x(7);
Q(7) = x(6)-2*c10*(x(7)-1/(x(7)^3))*l1k*(h2k+0.5*t2k)^2*l*0.5;
end