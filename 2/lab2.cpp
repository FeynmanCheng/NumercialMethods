#include<iostream>
#include<iomanip>
#include<cmath>
#define N 6
using namespace std;
double xy[N][2]={0,1,1,2.7,2,5.8,3,6.6,4,7.5,5,9.9};
void makematrix(int n);
double sumright(double a[N][2],double i);
double sumleft(double a[N][2],double i);
void mtrxmlt(double n[N][N],double m[N]);
double getA(double arcs[N][N],int n);
void getAStart(double arcs[N][N],int n,double ans[N][N]);
void mtrxmlt(double n[N][N],double m[N],double ans[]);
double offset(double a[],int n);
main(){
	int i=1;
	while(i<N+1){
		makematrix(i);//调用该函数 求不同阶数的拟合函数 
		i++;
	}
}
double offset(double a[],int n){
	double num[6]={0},dvt=0;
	int i=0,j=0;
	while(j<6){//求每一点的拟合值 
		for(i=0;i<n;i++){
		num[j]+=a[i]*pow((double)j,i);
		}
		//cout<<j<<" "<<setw(8)<<num[j]<<endl;;
		j++;
	}
	i=0;
	while(i<6){//对每一点产生的误差平方求和 
		dvt+=pow((num[i]-xy[i][1]),2);
		i++;
	}
	return dvt;
}
void makematrix(int n){
	/*构造n阶矩阵方程*/ 
	double sqare[N][N]={0},vector[N]={0};
	sqare[0][0]=N;
	int a,b;
	for(a=0;a<n;a++){
			for(b=a;b<n;b++){
					sqare[a][b]=sumleft(xy,a+b);
					sqare[b][a]=sqare[a][b];
			}
	vector[a]=sumright(xy,a);
	}
	/*解方程，方程两边左乘方阵的逆*/
	double inverse[N][N]={0},cmp[N][N]={0},ans[N]={0};//ans存放求解出的列向量 
	getAStart(sqare,n,cmp);
	double det=getA(sqare,n); 
	for(a=0;a<n;a++){
		for(b=0;b<n;b++){
			inverse[a][b]=cmp[a][b]/det;
		}
	} 
	mtrxmlt(inverse,vector,ans);
	/*求误差*/ 
	double result=offset(ans,n);
	cout<<"The polynominal is:"; 
	cout<<ans[0]<<"+";
	for(a=1;a<n;a++)
		a==n-1?cout<<ans[a]<<"x^"<<a:cout<<ans[a]<<"x^"<<a<<"+";
	cout<<endl;
	cout<<n-1<<" grades' offset:"<<result<<endl;	
}
double sumleft(double a[N][2],double i){//求和函数，结果是方程左边的方阵的元素 
	int j;
	double val=0;
	for(j=0;j<N;j++)
		val+=pow(xy[j][0],i);
	return val;
}
double sumright(double a[N][2],double i){//求和函数，结果是方程右边的列向量的元素 
	int j;
	double val=0;
	for(j=0;j<N;j++)
		val+=pow(xy[j][0],i)*xy[j][1];
	return val;
}
void mtrxmlt(double n[N][N],double m[N],double ans[])//矩阵乘法 
{
	int i=0,j=0;
	for (i=0;i<N;i++)
		for (j=0;j<N;j++)
			ans[i]+=n[i][j]*m[j];
} 
double getA(double arcs[N][N],int n)//矩阵:按第一行展开递归计算|A|
{
	if(n==1)
		return arcs[0][0];
	double ans = 0;
	double temp[N][N];
	int i,j,k;
	for(i=0;i<n;i++){
		for(j=0;j<n-1;j++)
			for(k=0;k<n-1;k++)
				temp[j][k] = arcs[j+1][(k>=i)?k+1:k];	
		double t = getA(temp,n-1);
		if(i%2==0)
			ans += arcs[0][i]*t;
		else
			ans -= arcs[0][i]*t;
	}
		
	return ans;
}
void getAStart(double arcs[N][N],int n,double ans[N][N])//矩阵：计算每一行每一列的每个元素所对应的余子式，组成A*
{
	if(n==1){
		ans[0][0] = 1;
		return;
	}
		
	int i,j,k,t;
	double temp[N][N];
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			for(k=0;k<n-1;k++)
				for(t=0;t<n-1;t++)
					temp[k][t] = arcs[k>=i?k+1:k][t>=j?t+1:t];
			ans[j][i]  =  getA(temp,n-1);
			if((i+j)%2 == 1)
				ans[j][i] = - ans[j][i];
		}
	}
		
}

