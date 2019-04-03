A=xlsread('C:\lagrange.xlsx');
B=xlsread('C:\qiebixuefu.xlsx');
C=xlsread('C:\lab1.xlsx');
a1=A(:,1);
b1=A(:,2);
a2=B(:,1);
b2=B(:,2);
a3=C(:,1);
b3=C(:,2);
x=-5:0.1:5;
y=1./(1+x.*x);
plot(x,y,'k');
hold on;
plot(a1,b1);
plot(a2,b2);
plot(a3,b3);
legend('origin','lagrange','chebyshev','spline');



