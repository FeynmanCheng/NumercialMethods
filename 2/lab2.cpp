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
		makematrix(i);//���øú��� ��ͬ��������Ϻ��� 
		i++;
	}
}
double offset(double a[],int n){
	double num[6]={0},dvt=0;
	int i=0,j=0;
	while(j<6){//��ÿһ������ֵ 
		for(i=0;i<n;i++){
		num[j]+=a[i]*pow((double)j,i);
		}
		//cout<<j<<" "<<setw(8)<<num[j]<<endl;;
		j++;
	}
	i=0;
	while(i<6){//��ÿһ����������ƽ����� 
		dvt+=pow((num[i]-xy[i][1]),2);
		i++;
	}
	return dvt;
}
void makematrix(int n){
	/*����n�׾��󷽳�*/ 
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
	/*�ⷽ�̣�����������˷������*/
	double inverse[N][N]={0},cmp[N][N]={0},ans[N]={0};//ans��������������� 
	getAStart(sqare,n,cmp);
	double det=getA(sqare,n); 
	for(a=0;a<n;a++){
		for(b=0;b<n;b++){
			inverse[a][b]=cmp[a][b]/det;
		}
	} 
	mtrxmlt(inverse,vector,ans);
	/*�����*/ 
	double result=offset(ans,n);
	cout<<"The polynominal is:"; 
	cout<<ans[0]<<"+";
	for(a=1;a<n;a++)
		a==n-1?cout<<ans[a]<<"x^"<<a:cout<<ans[a]<<"x^"<<a<<"+";
	cout<<endl;
	cout<<n-1<<" grades' offset:"<<result<<endl;	
}
double sumleft(double a[N][2],double i){//��ͺ���������Ƿ�����ߵķ����Ԫ�� 
	int j;
	double val=0;
	for(j=0;j<N;j++)
		val+=pow(xy[j][0],i);
	return val;
}
double sumright(double a[N][2],double i){//��ͺ���������Ƿ����ұߵ���������Ԫ�� 
	int j;
	double val=0;
	for(j=0;j<N;j++)
		val+=pow(xy[j][0],i)*xy[j][1];
	return val;
}
void mtrxmlt(double n[N][N],double m[N],double ans[])//����˷� 
{
	int i=0,j=0;
	for (i=0;i<N;i++)
		for (j=0;j<N;j++)
			ans[i]+=n[i][j]*m[j];
} 
double getA(double arcs[N][N],int n)//����:����һ��չ���ݹ����|A|
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
void getAStart(double arcs[N][N],int n,double ans[N][N])//���󣺼���ÿһ��ÿһ�е�ÿ��Ԫ������Ӧ������ʽ�����A*
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

