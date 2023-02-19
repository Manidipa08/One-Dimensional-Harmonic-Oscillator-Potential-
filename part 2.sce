clc
clear
n=2001//input("enter number of datapoints : ")
h=6.63D-34
b=50
rmin=-b*1D-6//input("enter r_min : ")
rmax=b*1D-6//input("enter r_max : ")
counter = 200
e=1.6D-19 
m=1.67*10^(-27)
//vt=xt_2
//v0=-100
E0 = 9*1D9
r=linspace(rmin,rmax,n)//range
new = r(2:n-1)*1D-9
rnew = new'
dr=((rmax-rmin)/(n-1))*1D-9
C=-((h/(2*%pi))^2)/(2*m*(dr^2)*e)
//rtp calculation
w=5.34*10^(21)
for i=1:4
    rtp(i)=sqrt((2*i-1)*(h/(2*%pi))/(m*w))
end
V=zeros(1,n)
for i=1:n
    V(i)=(1/2)*m*w*w*r(i)*1D-9.*r(i)*1D-9/(e)
end
A=eye(n-2,n-2)
v=diag(V(2:n-1))
//disp(A)
//----------------------By inbuilt command-------------------------
//tic()
D=(-2*C)*ones(n-2,1)
A1=diag(D)
//disp(A1)
C1=C*ones(n-3,1)
A2=diag(C1,1)
//disp(A2)
A3=diag(C1,-1)
//disp(A3)
//Hermitian matrix (tridiagonal)
H=A1+A2+A3+v
//disp("By inbuilt : ",H)
[a1 a2]=spec(H)
Z=spec(H)
//disp("Corresponding eigen vectors : ",a1)
//disp("Eigen values : ",a2)
//counter = 0
//for i=1:n-2
    //if a2(i,i)<0
        //counter=counter+1
    //end
//end
//disp(counter)

//--------------------------Normalizations & Expectation values-------------------------
vnew = V(2:n-1)
for i = 1:counter
    eigenvec = a1(:,i).*a1(:,i)
    Nrm(i)=inttrap(rnew,eigenvec)
    u_r(:,i)=(1/sqrt(Nrm(i))).*a1(:,i)
    neigen(:,i) = u_r(:,i).*u_r(:,i)
end
plot(rnew*1D9,neigen(:,130),'b')
title("<<|u_r(130)|^2 vs r(in nm) Plot >>",'color','brown','font',2,'fontsize',4)
xlabel("r(in nm)--------------->",'color','brown','font',2,'fontsize',4)
ylabel("|u_r(130)|^2----------->",'color','brown','font',2,'fontsize',4)
