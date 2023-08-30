// RK4 for solving a system of ODE//
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
float f(float xn,float yn,float zn) // function f(x)=y'
{
  float f;
  f=2*(yn/xn)-xn*zn;
  return f;
}
float g(float xn,float yn,float zn) // function g(x)=z'
{
  float g;
  g=(zn+yn)/xn;
  return g;
}
float lx(float xn) // analytical y
{
  float l;
  l=pow(xn,2)*cos(xn);
  return l;
}
float mx(float xn) // analytical z
{
  float m;
  m=xn*sin(xn);
  return m;
}
float* yRK4(float xo,float yo,float zo,float L,int N) //runge-kutta method of forth order
{
  float h,eyL2,ezL2,xn,yn,zn,zn1,yn1,esumy,esumz,k1,k2,k3,k4,l1,l2,l3,l4;
  float Yval[100],Zval[100];
  static float valy[100];
  Yval[0]=yo;
  Zval[0]=zo;
  esumy=esumz=0;
  h=(L-xo)/N;
  xn=xo;
  yn=yo;
  zn=zo;
  for(int i=1;i<=N;i++)
   {
     k1=h*f(xn,yn,zn);
     l1=h*g(xn,yn,zn);
     k2=h*f(xn+h/2,yn+k1/2,zn+l1/2);
     l2=h*g(xn+h/2,yn+k1/2,zn+l1/2);
     k3=h*f(xn+h/2,yn+k2/2,zn+l2/2);
     l3=h*g(xn+h/2,yn+k2/2,zn+l2/2);
     k4=h*f(xn+h,yn+k3,zn+l3);
     l4=h*g(xn+h,yn+k3,zn+l3);
     yn1=yn+(k1+2*k2+2*k3+k4)/6;
     zn1=zn+(l1+2*l2+2*l3+l4)/6;
     valy[i]=Yval[i]=yn1;
     Zval[i]=zn1;
     esumy=esumy+pow((yn1-lx(xn+h)),2);
     esumz=esumz+pow((zn1-mx(xn+h)),2);
     xn=xn+h;
     yn=yn1;
     zn=zn1;
   }
  eyL2=sqrt(esumy)/N;
  ezL2=sqrt(esumz)/N;
  printf("%f \t %f \t %f \t %f\n ",Yval[N],Zval[N],eyL2,ezL2);
  return valy;
}
float* zRK4(float xo,float yo,float zo,float L,int N) //runge-kutta method of forth order
{
  float h,eyL2,ezL2,xn,yn,zn,zn1,yn1,esumy,esumz,k1,k2,k3,k4,l1,l2,l3,l4;
  float Yval[100],Zval[100];
  static float valz[100];
  Yval[0]=yo;
  Zval[0]=zo;
  esumy=esumz=0;
  h=(L-xo)/N;
  xn=xo;
  yn=yo;
  zn=zo;
  for(int i=1;i<=N;i++)
   {
     k1=h*f(xn,yn,zn);
     l1=h*g(xn,yn,zn);
     k2=h*f(xn+h/2,yn+k1/2,zn+l1/2);
     l2=h*g(xn+h/2,yn+k1/2,zn+l1/2);
     k3=h*f(xn+h/2,yn+k2/2,zn+l2/2);
     l3=h*g(xn+h/2,yn+k2/2,zn+l2/2);
     k4=h*f(xn+h,yn+k3,zn+l3);
     l4=h*g(xn+h,yn+k3,zn+l3);
     yn1=yn+(k1+2*k2+2*k3+k4)/6;
     zn1=zn+(l1+2*l2+2*l3+l4)/6;
     Yval[i]=yn1;
     valz[i]= Zval[i]=zn1;
     esumy=esumy+pow((yn1-lx(xn+h)),2);
     esumz=esumz+pow((zn1-mx(xn+h)),2);
     xn=xn+h;
     yn=yn1;
     zn=zn1;
   }
  eyL2=sqrt(esumy)/N;
  ezL2=sqrt(esumz)/N;
//printf("%f \t %f \t %f \t %f\n ",Yval[N-1],Zval[N-1],eyL2,ezL2);
  return valz;
}
int main(){
// declaring variables
float  *valy;
float *valz;
float  *vay;
float *vaz;
float  *valuey;
float *valuez;
int k;
float Yan[6],Zan[6];
float pi =3.1416;
float j=pi/2;
FILE *fptr;
fptr=fopen("result_2a.txt","w"); // to write the results in a file
//to find y,z using RK4 for N=5
valy=yRK4(pi/2,0,pi/2,pi,5);
valz=zRK4(pi/2,0,pi/2,pi,5);
fprintf(fptr,"Table for N=5 \n x \tAnalytical,y(x)\t yn(RK4)\t Analytical,z(x)\tzn(RK4)\n");
for(int i=0;i<6;i++)
 {
   Yan[i]=lx(j);
   Zan[i]=mx(j);
   fprintf(fptr,"%.3f \t %f\t%f\t %f \t %f \n",j,Yan[i],valy[i],Zan[i],valz[i]);
   j=j+pi/10;
 }
//to find y,z using RK4 for N=10
vay=yRK4(pi/2,0,pi/2,pi,10);
vaz=zRK4(pi/2,0,pi/2,pi,10);
j=pi/2;
k=0;
fprintf(fptr,"Table for N=10 \n x \tAnalytical,y(x)\t yn(RK4)\t Analytical,z(x)\tzn(RK4)\n");
for(int i=0;i<11;i=i+2)
 {
   fprintf(fptr,"%.3f \t %f\t%f\t %f \t %f \t\n",j,Yan[k],vay[i],Zan[k],vaz[i]);
   j=j+pi/10;
   k=k+1;
 }
//to find y,z using RK4 for N=20
valuey=yRK4(pi/2,0,pi/2,pi,20);
valuez=zRK4(pi/2,0,pi/2,pi,20);
j=pi/2;
k=0;
fprintf(fptr,"Table for N=20 \n x \tAnalytical,y(x)\t yn(RK4)\t Analytical,z(x)\t zn(RK4)\n");
for(int i=0;i<21;i=i+4)
 {
   fprintf(fptr,"%.3f\t %f\t%f\t %f \t %f\n",j,Yan[k],valuey[i],Zan[k],valuez[i]);
   j=j+pi/10;
   k=k+1;
 }
fclose(fptr);
printf("Results are written in a file named result_2a");
return 0;
}
