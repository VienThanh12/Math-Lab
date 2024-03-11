function [y,dy,d2y,tau]=second_degree_homogeneous(a,b,c,y0,y1)
%y=y(t) is the solution of the second degree homogeneous equation
%a*y'(t)+b*y(t)+c*y(t)=0,y(0)=y0,y'(0)=y1
%dy=y'(t)
%d2y=y''(t)
%tau is time constant
D=b^2-4*a*c

if D<0
    alpha=-b/(2*a);
    tau=-1/alpha;
    w=sqrt(-D)/(2*a);
    C1=y0;
    C2=(y1-alpha*y0)/w;
    C=sqrt(C1^2+C2^2);
    th=atan2(C1,C2);
    y=@(t) C*exp(alpha*t).*sin(w*t+th);
    dy=@(t) C*exp(alpha*t)*alpha.*sin(w*t+th)+C*exp(alpha*t).*cos(w*t+th)*w;
    d2y=@(t) C*exp(alpha*t)*alpha^2.*sin(w*t+th) +2*C*exp(alpha*t)*alpha.*cos(w*t+th)*w-C*exp(alpha*t).*sin(w*t+th)*w^2;





elseif D==0
    la=-b/(2*a);
    tau=-1/la;
    C1=y0;
    C2=y1+la*y0;
    y=@(t) C1*exp(la*t)+C2*t.*exp(la*t);
    dy=@(t) C1*exp(la*t)*la+C2*exp(la*t)+C2*t.*exp(la*t)*la;
    d2y=@(t) C1*exp(la*t)*la^2+2*C2*la*exp(la*t) + C2*t.*exp(la*t)*la^2;


else
    la1=(-b+sqrt(D))/(2*a);
    la2=(-b-sqrt(D))/(2*a);
    tau=-1/la1;
    C1=(y1-la2*y0)/(la1-la2);
    C2=(y1-la1*y0)/(la2-la1);
    y=@(t) C1*exp(la1*t)+C2*exp(la2*t);
    dy=@(t) C1*exp(la1*t)*la1+C2*exp(la2*t)*la2;
    d2y=@(t) C1*exp(la1*t)*la1^2+C2*exp(la2*t)*la2^2;
end
