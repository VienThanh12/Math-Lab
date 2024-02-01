clear
close all
x=[1, 4, 3, 1, 0]
y=[1, 2, 4, 5, 3]
n=length(x)

%tapaus=case in finnish
%tapaus=1 %1=natural, 2=clamped, 3=periodic

%tapaus=3 clamped spline
%k1=1 %Y'(x1)
%kn=2 %Y'(xn)

tapaus = 3
% then requirement y(1)=y(n)


%coefficients from a function spline_curves_coefficients.m
%which has to be on the current folder
%see spline_curves_coefficients.pdf
n=length(x)

%tapaus=3 %1=natural, 2=clamped, 3=periodic
%system of equations for  zk=Y''(xk)

A=zeros(n,n); %coefficient matrix
B=zeros(n,1); %right hand side
h=x(2:end)-x(1:end-1); %lenghts of the intervals x_k...x_(k+1)

for k=2:n-1
  alfak=h(k-1)/6;
  betak=(h(k-1)+h(k))/3;
  gammak=h(k)/6;
  deltak=(y(k+1)-y(k))/h(k)-(y(k)-y(k-1))/h(k-1);
  A(k,:)=[zeros(1,k-2),alfak,betak,gammak,zeros(1,n-(k-2)-3)];
  B(k)=deltak;
end


if tapaus==1
%1) natural spline

A(1,:)=[1,zeros(1,n-1)];
A(n,:)=[zeros(1,n-1),1];
B(1)=0;
B(n)=0;

elseif tapaus==2
%2) clamped
%k1=0 %Y'(x1)
%kn=-2 %Y'(xn)

alfa1=h(1)/3
beta1=h(1)/6
delta1=(y(2)-y(1))/h(1)-k1
alfan=h(n-1)/6
betan=h(n-1)/3
deltan=kn-(y(n)-y(n-1))/h(n-1)

A(1,:)=[alfa1,beta1,zeros(1,n-2)];
A(n,:)=[zeros(1,n-2),alfan,betan];
B(1)=delta1;
B(n)=deltan;

else
%3) periodic

alfa1=h(1)/3
beta1=h(1)/6
mu1=h(n-1)/6
la1=h(n-1)/3
delta1=(y(2)-y(1))/h(1)-(y(n)-y(n-1))/h(n-1)


A(1,:)=[alfa1,beta1,zeros(1,n-4),mu1,la1];
A(n,:)=[1,zeros(1,n-2),-1];
B(1)=delta1;
B(n)=0;
end

%solution i.e z1,z2,...,zn

z=A^-1*B

%coefficients for curves Yk=ak(x-xk)^3+bk(x-xk)^2+ck(x-xk)+dk, k=1,2,...,n-1,

for k=1:n-1
a(k)=(z(k+1)-z(k))/(6*h(k));
b(k)=z(k)/2;
c(k)=-h(k)/3*z(k)-h(k)/6*z(k+1)+(y(k+1)-y(k))/h(k);
d(k)=y(k);
end

abcd=[a;b;c;d]'

%Y_k=a_k*x^3+b_k*x^2+c_k*x+d_k
%coefficients a_k,b_k,c_k,d_k
%as rows in matrix abcd
a=abcd(:,1);
b=abcd(:,2);
c=abcd(:,3);
d=abcd(:,4);
%% graph
figure(1)
plot(x,y,'r.','markersize',18)
hold on
for k=1:n-1
  dx=(x(k+1)-x(k))/100;
  xx=x(k):dx:x(k+1);
  yy=a(k)*(xx-x(k)).^3+b(k)*(xx-x(k)).^2+c(k)*(xx-x(k))+d(k);
  plot(xx,yy,'linewidth',1.5)
end
hold off
grid
xlabel('x','fontsize',14)
ylabel('y','fontsize',14)
if tapaus==1
  title('natural spline','fontsize',10)
elseif tapaus==2
  title('clamped spline','fontsize',10)
else
  title('periodic spline','fontsize',10)
end
