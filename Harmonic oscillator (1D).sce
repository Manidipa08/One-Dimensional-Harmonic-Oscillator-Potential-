clc
clear
n=4001//input("enter number of datapoints : ")
h=6.63D-34
b=12
rmin=-b*1D-6//input("enter r_min : ")
rmax=b*1D-6//input("enter r_max : ")
counter = 4
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
w=5.34D+21
for i=1:4
    rtp(i)=sqrt((2*i-1)*(h/(2*%pi))/(m*w))
end
//disp("Classical turning points : ",rtp)
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
Z0=(((1:(n-2))-0.5)*(h/(2*%pi))*w)/e
//disp("E_Th, E_comp & their ratio for first four values of Energy- ")
//disp(" E_th         E_comp Energy            eigen value Ratio ")
//disp(string(Z0(1:5)')+"        "+string(Z(1:5))+"         "+string(Z(1:5)./Z0'(1:5)))
//--------------------------Normalizations & Expectation values-------------------------
vnew = V(2:n-1)
for i = 1:100
    eigenvec = a1(:,i).*a1(:,i)
    Nrm(i)=inttrap(rnew,eigenvec)
    u_r(:,i)=(1/sqrt(Nrm(i))).*a1(:,i)
    neigen(:,i) = u_r(:,i).*u_r(:,i)
    neigen1(:,i) = rnew.*neigen(:,i)
    neigen2(:,i) = rnew.*neigen1(:,i)
    neigen3(:,i) = vnew'.*neigen(:,i)
end
for i=1:counter
    ex_r(i)=inttrap(rnew,neigen1(:,i))
    ex_r2(i)= inttrap(rnew,neigen2(:,i))
    ex_V(i) = inttrap(rnew,neigen3(:,i))
end
for i = 1:n-3
    for j = 1:counter
        r2(i) = (rnew(i)+rnew(i+1))/2
        mid_u(i,j) = (u_r(i,j)+u_r(i+1,j))/2
        diff_u(i,j) = (u_r(i+1,j)-u_r(i,j))/dr
    end
end

for i = 1:n-4
    for j = 1:counter
        r3(i) = (r2(i) + r2(i+1))/2
        mid2_u(i,j) = (u_r(i+2,j) + 2*u_r(i+1,j)+ u_r(i,j))/4
        diff2_u(i,j) = (u_r(i+2,j) - 2*u_r(i+1,j) + u_r(i,j))/(dr*dr)
    end
end
//----------------Uncertainty check--------------------------------------
Un = (4*%pi)/h

for i=1:counter
    y2(:,i) = mid_u(:,i).*diff_u(:,i)
    y3(:,i) = mid2_u(:,i).*diff2_u(:,i)
    ex_p(i) = -1*%i*(h/(2*%pi))*inttrap(r2,y2(:,i))
    ex_p2(i) = -1*(h/(2*%pi))**2*inttrap(r3,y3(:,i))
    sig_r(i) = sqrt(ex_r2(i) - (ex_r(i)*ex_r(i)))
    sig_p(i) = sqrt(ex_p2(i) - (ex_p(i)*ex_p(i)))
    ex_K(i) = ex_p2(i)/(2*m*e)
    un(i) = Un*(sig_r(i).* sig_p(i))
end
//-----------------total energy------------------------
E = ex_V + ex_K
disp("Expectation value <r> =",ex_r)
disp("Expectation value <r2> =",ex_r2)
disp("Expectation value <p> =",ex_p)
disp("Expectation value <p2> =",ex_p2)
disp("Expectation value <V> =",ex_V)
disp("Expectation value <KE> =",ex_K)
disp("Standard deviation of r =",sig_r)
disp("Standard deviation of p =",sig_p)
disp("Uncertanity Product (hbar/2) = ",un)
disp("E_Th, E_comp(<KE>+<V>) for first four values of Energy- ")
disp(" E_th         E_comp (<KE>+<V>) Energy")
disp(string(Z0(1:4)')+"        "+string(E(1:4)))
show_window(1)
for i=1:counter
    subplot(2,2,i)
    plot(rnew*1D9,u_r(:,i),'m','linewidth',2)
    title("<<u_r vs r(in nm) Plot for bound state "+string(i)+">>",'color','brown','font',2,'fontsize',4)
    xlabel("r(in nm)--------------->",'color','brown','font',2,'fontsize',4)
    ylabel("u_r----------->",'color','brown','font',2,'fontsize',4)
    xgrid()
end
show_window(2)
for i=1:counter
    subplot(2,2,i)
    plot(rnew*1D9,neigen(:,i),'k','linewidth',2)
    title("<<|u_r|^2 vs r(in nm) Plot for bound state "+string(i)+">>",'color','brown','font',2,'fontsize',4)
    xlabel("r(in nm)--------------->",'color','brown','font',2,'fontsize',4)
    ylabel("|u_r|^2----------->",'color','brown','font',2,'fontsize',4)
    xgrid()
end
//Comparison of energy eigen value ratio
/*
Z = Z(1:counter)
Z0=(1/2).*6.582*10^(-16).*w
ratio=Z0/Z(1)
disp("Ratio : ",ratio)
*/
//Hermite Polynomial
show_window(3)
si=sqrt(m*w/(h/(2*%pi)))*(r*1D-9)
he(:,1)=si./si
he(:,2)=2.*si
he(:,3)=-4*si.*si+2
he(:,4)=-8*si.*si.*si+12*si

for i=1:n
    hermite_prac(i)=exp(-si(i)*si(i)/2)
end
for i=1:counter
    subplot(2,2,i)
    plot(rnew,u_r(:,i),'-*c')
    Prac_si(i,:)=(((m*w)/(3.1416*(h/(2*%pi))))^(1/4))*(1/sqrt((2^(i-1))*factorial(i-1)))*(he(:,i).*hermite_prac)
    subplot(2,2,i)
    plot((r*1D-9),Prac_si(i,:))
    at=gca()
    at.tight_limits=["on"]
    at.x_ticks=tlist(["ticks","locations","labels"],[-1D-14,0,1D-14],["-1D-14","0","1D-14"])
    title('<- u_r('+string(i)+') vs r plot ->','color','blue','Fontsize','4')
 xlabel('Position(r) --->>>>','color','brown','Fontsize','4')
 ylabel('u_r --->>>>','color','brown','Fontsize','4')
 legend('Theoritical','Computational')
end
