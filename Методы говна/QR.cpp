#include <stdio.h>
#include <math.h>
#include <time.h>
#define _N 10 // число уравнений
#define _M 10 // число неизвестных

float **a;
float *b;
float *x;
int n,m;
// QR-алгоритм
void qr () {
int l,k;
float c[_N][_M];
float s[_N][_M];
float akk,akl,alk,all,bk,bl;

// Прямой ход
for (k=0; k<n-1; k++) {
for (l=k+1; l<n; l++) {
c[k][l] = a[k][k] / (sqrt( a[k][k]*a[k][k] + a[l][k]*a[l][k] ));
s[k][l] = a[l][k] / (sqrt( a[k][k]*a[k][k] + a[l][k]*a[l][k] ));

// Умножение матрицы a[][] на T[k][l]
akk=a[k][k]; alk=a[l][k]; akl=a[k][l]; all=a[l][l];
a[k][k] = akk*c[k][l] + alk*s[k][l];
a[k][l] = akl*c[k][l] + all*s[k][l];
a[l][k] = -akk*s[k][l] + alk*c[k][l];
a[l][l] = -akl*s[k][l] + all*c[k][l];

// Вектор свободных членов умножается на T[k][l]
bk = b[k]; bl = b[l];
b[k] = bk*c[k][l] + bl*s[k][l];
b[l] = -bk*s[k][l] + bl*c[k][l];
}
}
// Теперь матрица a[][] — верхняя диагональная.

// Обратный ход
float h;
x[n-1]=b[n-1]/a[n-1][n-1];
for (l=(n-1);l>=1;l--) {
h=b[l-1];
for (k=(l+1);k<=n;k++) h=h-x[k-1]*a[l-1][k-1];
x[l-1]=h/a[l-1][l-1];
}

}
// —--------------------------------------------------—
// —--------------------------------------------------—
float matrix(int N,int M) {
n = N; m = M;
a = new float* [n];
for (int i=0; i<n; i++) {
a[i] = new float [m];
for (int j=0; j<m; j++) a[i][j]=0;
b = new float [n];
x = new float [n];
{ b[i]=0; x[i]=0; }}

}

// —--------------------------------------------------—
// Вывод СЛАУ на экран
// —--------------------------------------------------—
void print_slau() {
for (int i=0; i<n; i++) {
for (int j=0; j<m; j++) printf("%.2f ",a[i][j]);
printf("| %.3f ",b[i]);
printf("\n");
}
}

// —--------------------------------------------------—
// Печать решения и невязки решения
// —--------------------------------------------------—
void print_x() {
float max=0,h;
for (int i=0; i<n; i++) {
h=0;
for (int j=0; j<n; j++) h=h+x[j]*a[i][j];
if (max<fabs(b[i]-h)) max=fabs(b[i]-h);
printf("x[%i]=%.3f ",i,x[i]);
if (i==5) printf("\n");
}
printf("\nMax. Nevyaka: %f\n",max);
}
// —--------------------------------------------------—
// Установить значение матрицы a[][]
// —--------------------------------------------------—
void seta(int i, int j, float value) {
a[i][j] = value;
}
// —--------------------------------------------------—
// Установить значение вектора b[]
// —--------------------------------------------------—
void setb(int i, float value) {
b[i] = value;
}

// —--------------------------------------------------—

int main() {
double t;

matrix (_N,_M);
for (int i=0; i<_N; i++)
{for (int j=0; j<_M; j++)
seta(i,j,i+j);
setb(i,1);}
printf("Solusion (QR-razlojeniy)\n");
printf("=Ishodnaya sistema========================================\n");
print_slau();
printf("\n=Solusion=========================================\n");
t=clock();
qr();
t=clock()-t;
printf("Vremayraboti funk %lqr sekundi\n",t/CLOCKS_PER_SEC);
print_x();

}