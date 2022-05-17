
disp('**************************')
disp(' With this program you will transformate WGS84---> Hatt (Old GR Datum)')
disp('Put the geodetic coordinates of the center of the sheet')
disp('**************************')
Latitude_degrees=input('Give the Degrees of the Latitude = ')
Latitude_minutes=input('Give the Minutes of the Latitude = ')
disp('**************************')
disp('NOTE : if you must input negative values,please type (-) before degrees and
minutes e.g (-1 deg -45 min)')
Longtitude_degrees=input('Give the Degrees of the Longtitude = ')
Longtitude_minutes=input('Give the Minutes of the Longtitude = ')
C_Rad_Lat=(Latitude_degrees)*(pi/180)+(Latitude_minutes/60)*(pi/180);
C_Rad_Lon=(Longtitude_degrees)*(pi/180)+(Longtitude_minutes/60)*(pi/180)+
0.413928175889201;
load CAR.m
V=CAR;
Q=size(V);
n=Q(1);
% Approximate transformation %
for i=1:n
V(i,2)=V(i,2)-456;
V(i,3)=V(i,3)-373;
V(i,4)=V(i,4)-497;
end
e1=0.006674372;
e2=0.006719219;
a=6377397.155;
for i=1:n
lon(i)=atan(V(i,3)/V(i,2));
end
Rad_lon=lon';
Deg_lon=Rad_lon*(180/pi);
for i=1:n
P(i)=sqrt((V(i,2)^2)+(V(i,3)^2));
end
P=P';
for i=1:n
fo(i)=atan((V(i,4)*(1+e2))/P(i));
end
fo=fo';
for i=1:n
Wo(i)=1/(sqrt(1-e1*(sin(fo(i))^2)));
end
Wo=Wo';
for i=1:n
No(i)=a*Wo(i);
21/01/06 22:34 C:\MATLAB7\gps\Hatt.m 2 of 7
end
No=No';
for i=1:n
ho(i)=(V(i,4)/sin(fo(i)))-(1-e1)*No(i);
end
ho=ho';
for i=1:n
Po(i)=(No(i)+ho(i))*cos(fo(i));
end
Po=Po';
for i=1:n
hh(i)=ho(i)+(P(i)-Po(i))*cos(fo(i));
end
hh=hh';
for i=1:n
f(i)=atan((tan(fo(i)))/(1+e2*(hh(i)/(No(i)+hh(i)))));
end
Rad_f=f';
Deg_f=Rad_f*(180/pi);
for i=1:n
W(i)=1/(sqrt(1-e1*(sin(f(i))^2)));
end
W=W';
for i=1:n
N(i)=a*W(i);
end
N=N';
for i=1:n
h(i)=(V(i,4)/sin(f(i)))-(1-e1)*N(i);
end
h=h';
for i=1:n
f_l_h_GEOCENTRIC_Radians(i,1)=[V(i,1)];
f_l_h_GEOCENTRIC_Radians(i,2)=[Rad_f(i)];
f_l_h_GEOCENTRIC_Radians(i,3)=[Rad_lon(i)];
f_l_h_GEOCENTRIC_Radians(i,4)=[h(i)];
end
for i=1:n
phi_lam_h_GEOCENTRIC_Degrees(i,1)=[V(i,1)];
phi_lam_h_GEOCENTRIC_Degrees(i,2)=[Deg_f(i)];
phi_lam_h_GEOCENTRIC_Degrees(i,3)=[Deg_lon(i)];
phi_lam_h_GEOCENTRIC_Degrees(i,4)=[h(i)];
end
WC=1/(sqrt(1-e1*(sin(C_Rad_Lat)^2)));
NC=a*WC;
MC=a*(1-e1)*(WC^3);
A1=NC*cos(C_Rad_Lat);
A2=(-MC*sin(C_Rad_Lat));
A3=((-MC*cos(C_Rad_Lat))/6)*(2+9*e2*((sin(C_Rad_Lat))^2));
A4=(-NC*cos(C_Rad_Lat)*((sin(C_Rad_Lat))^2))/6;
A5=(-NC*sin(C_Rad_Lat)*(1-2*(sin(C_Rad_Lat)^2)))/6;
B1=MC;
B2=NC*cos(C_Rad_Lat)*sin(C_Rad_Lat)/2;
B3=(3*e2*(MC^2)*cos(C_Rad_Lat)*sin(C_Rad_Lat))/(2*NC);
B4=MC*(1-4*(sin(C_Rad_Lat)^2)+e2*(cos(C_Rad_Lat)^4))/6;
B5=e2*MC*(1-2*(sin(C_Rad_Lat)^2))/2;
B6=(NC*sin(C_Rad_Lat)*cos(C_Rad_Lat)*(1-2*(sin(C_Rad_Lat)^2)));
B7=(-NC*sin(C_Rad_Lat)*cos(C_Rad_Lat))/3;
for i=1:n
df(i)=Rad_f(i)-C_Rad_Lat;
dl(i)=Rad_lon(i)-C_Rad_Lon;
end
df=df';
dl=dl';
for i=1:n
x(i)=A1*dl(i)+(A2*dl(i)*df(i))+(A3*dl(i)*df(i)^2)+(A4*dl(i)^3)+A5*(dl(i)^3)*df(i);
y(i)=B1*df(i)+B2*(dl(i)^2)+B3*(df(i)^2)+B4*(df(i)*dl(i)^2)+B5*(df(i)^3)+B6*(dl(i)
^4)+B6*(df(i)^2)*(dl(i)^2);
end
x=x';
y=y';
disp(' Do you want to find the transformations parameters (similarity transformation)')
disp(' 1:YES 2:NO')
ANSWER=input('ANSWER :')
% LOADING COORDINATES OF LOCAL DATUM %
if ANSWER==1
load HATT_COR_2.m
load itrf_2.m
G=HATT_COR_2;
IT=itrf_2;
si=size(IT);
li=si(1);
for i=1:li
AP(i,1)=IT(i,1);
AP(i,2)=IT(i,2)-456;
AP(i,3)=IT(i,3)-373;
AP(i,4)=IT(i,4)-497;
end
% finding longtitude (Greece) %
for i=1:li
lont(i)=atan(AP(i,3)/AP(i,2));
end
%radians%
Rad_lont=lont';
%degrees%
Deg_lont=Rad_lont*(180/pi);
for i=1:li
Pt(i)=sqrt((AP(i,2)^2)+(AP(i,3)^2));
end
Pt=Pt';
for i=1:li
fot(i)=atan((AP(i,4)*(1+e2))/Pt(i));
end
fot=fot';
for i=1:li
Wot(i)=1/(sqrt(1-e1*(sin(fot(i))^2)));
end
Wot=Wo';
for i=1:li
Not(i)=a*Wot(i);
end
Not=Not';
for i=1:li
hot(i)=(AP(i,4)/sin(fot(i)))-(1-e1)*Not(i);
end
hot=hot';
for i=1:li
Pot(i)=(Not(i)+hot(i))*cos(fot(i));
end
Pot=Pot';
for i=1:li
hht(i)=hot(i)+(Pt(i)-Pot(i))*cos(fot(i));
end
hht=hht';
for i=1:li
ft(i)=atan((tan(fot(i)))/(1+e2*(hht(i)/(Not(i)+hht(i)))));
end
Rad_ft=ft';
Deg_ft=Rad_ft*(180/pi);
for i=1:li
Wt(i)=1/(sqrt(1-e1*(sin(ft(i))^2)));
end
Wt=Wt';
for i=1:li
Nt(i)=a*Wt(i);
21/01/06 22:34 C:\MATLAB7\gps\Hatt.m 5 of 7
end
Nt=Nt';
for i=1:li
ht(i)=(AP(i,4)/sin(ft(i)))-(1-e1)*Nt(i);
end
ht=ht';
for i=1:li
f_l_h_GEOCENTRIC_Radianst(i,1)=[AP(i,1)];
f_l_h_GEOCENTRIC_Radianst(i,2)=[Rad_ft(i)];
f_l_h_GEOCENTRIC_Radianst(i,3)=[Rad_lont(i)];
f_l_h_GEOCENTRIC_Radianst(i,4)=[ht(i)];
end
for i=1:li
phi_lam_h_GEOCENTRIC_Degreest(i,1)=[AP(i,1)];
phi_lam_h_GEOCENTRIC_Degreest(i,2)=[Deg_ft(i)];
phi_lam_h_GEOCENTRIC_Degreest(i,3)=[Deg_lont(i)];
phi_lam_h_GEOCENTRIC_Degreest(i,4)=[ht(i)];
end
for i=1:li
dft(i)=Rad_ft(i)-C_Rad_Lat;
dlt(i)=Rad_lont(i)-C_Rad_Lon;
end
dft=dft';
dlt=dlt';
for i=1:li
xt(i)=A1*dlt(i)+(A2*dlt(i)*dft(i))+(A3*dlt(i)*dft(i)^2)+(A4*dlt(i)^3)+A5*(dlt(i)^3)
*dft(i);
yt(i)=B1*dft(i)+B2*(dlt(i)^2)+B3*(dft(i)^2)+B4*(dft(i)*dlt(i)^2)+B5*(dft(i)^3)+B6*
(dlt(i)^4)+B6*(dft(i)^2)*(dlt(i)^2);
end
xt=xt';
yt=yt';
% Known %
xm=mean(xt);
ym=mean(yt);
for i=1:li
u(i)=xt(i)-xm;
v(i)=yt(i)-ym;
end
u=u';
v=v';
for i=1:li;
u2(i)=u(i)^2;
v2(i)=v(i)^2;
end
u2=u2';
v2=v2';
for i=1:li
R1(i)=u(i)*G(i,2);
R2(i)=v(i)*G(i,3);
R3(i)=v(i)*G(i,2);
R4(i)=u(i)*G(i,3);
R5(i)=u2(i)+v2(i);
end
f1=sum(R1);
f2=sum(R2);
f3=sum(R3);
f4=sum(R4);
f5=sum(R5);
c=(f1+f2)/f5
d=(f3-f4)/f5
sx=0;sy=0;
for i=1:li
sx=sx+G(i,2);
sy=sy+G(i,3);
end
sx=sx/li;
sy=sy/li;
tx=sx-c*xm-d*ym
ty=sy+d*xm-c*ym
for i=1:li
resE(i)=G(i,2)-(c*xt(i)+d*yt(i)+tx);
resN(i)=G(i,3)-(-d*xt(i)+c*yt(i)+ty);
end
resE=resE'
resN=resN'
end
disp('******************************')
disp('Do you want to continue?')
disp('1:YES 2:NO')
cont=input('cont = ')
if cont == 1
for i=1:n
X(i)=c*x(i)+d*y(i)+tx;
Y(i)=-d*x(i)+c*y(i)+ty;
end
X=X';
Y=Y';
TRANSCOR=[X Y];
else
disp('END')
end