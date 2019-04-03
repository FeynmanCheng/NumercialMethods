#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define _USE_MATH_DEFINES
#define N 8
#define eps 1e-6
double a1[N]={-5.0,-5+10.0/7,-5+20.0/7,-5+30.0/7,5-30.0/7,5-20.0/7,5-10.0/7,5.0};
double q[N]={0};//用于保存切比雪夫插值点 
double s[101]={0};//用于计算插值结果的保存 
double g[N]={0};//存放gk的值 
double ans[N]={0};//矩阵方程右边左乘左边方阵的逆，结果保存至此，即为M矩阵的值 
double runge(double x);
double ss[N-1][4]={0};//存放三次样条差值多项式系数 
double h[N-1]={0};//存放hk 
int DinV(double A[N][N]);
inline void swap(double &a,double &b){double c=a;a=b;b=c;};
void makes(double a[N-1][4],double an[])//生成插值函数的系数 
{
	int  i,j;
	for (i=0;i<N;i++){
		a[i][0]=(an[i+1]-an[i])/6/h[i];
		a[i][1]=an[i]/2-3*a[i][0]*a1[i];
		a[i][2]=((runge(a1[i+1])-runge(a1[i]))/h[i]-(2*an[i]+an[i+1])/6*h[i])\
					+3*a[i][0]*a1[i]*a1[i] -an[i]*a1[i];
		a[i][3]=runge(a1[i])-a1[i]*((runge(a1[i+1])-runge(a1[i]))/h[i]-(2*an[i]+an[i+1])/6*h[i])\
				+an[i]/2*a1[i]*a1[i]-a[i][0]*pow(a1[i],3);
	}
}
void mtrxmlt(double n[N][N],double m[N])//矩阵乘法 
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

void spline(double mtrx[N][N]){
	double lumda[N-1]={0},miu[N-1]={0};
	int i=0,j=0;
	while(i<N-1){
		h[i]=a1[i+1]-a1[i];
		i++;
	}
	for (i=1;i<N-1;i++){
		lumda[i]=h[i]/(h[i]+h[i-1]);
		miu[i]=h[i-1]/(h[i]+h[i-1]);
		g[i]=6*(runge(a1[i+1])-2*runge(a1[i])+runge(a1[i-1]))/20*7/10*7;
	}
	g[0]=6/h[0]*((runge(a1[1])-runge(a1[0]))/h[0]-10/pow(26,2));
	g[7]=6*(-(runge(a1[7])-runge(a1[6]))/h[6]-10/pow(26,2))/10*7;
	mtrx[0][1]=1;
	mtrx[N-1][N-2]=1;
	for(i=0;i<N;i++)
		mtrx[i][i]=2;
	for(i=1;i<N-1;i++){
		mtrx[i][i-1]=miu[i];
		mtrx[i][i+1]=lumda[i];
	}

}
double qiebixuefu(double s,double e,double i,double n){
	return (s+e)/2+(e-s)/2*cos((2*i+1)/(2*n+2)*M_PI);
}
double runge(double x){
	return 1/(1+pow(x,2));
}
double lagrange(double a,double t1[]){
	int i=0,j=0;
	double count=0,multi=1;
	while(i<8){
		j=0;
		multi=1;
		while(j<N){
			if(i!=j)
				multi*=(a-t1[j])/(t1[i]-t1[j]); 
			j++;
		}
		multi*=runge(t1[i])	;
	count+=multi;
	i++;
	}
	return count;
}
int main(){
	FILE *fp;
	int i=0;
	/*拉格朗日插值*/ 
	while(i<101){
		s[i]=lagrange(i*0.1-5,a1);
		i++;
	}
	fp=fopen("C:\\lagrange.txt","w+");
	for (i=0;i<101;i++)
		fprintf(fp,"%lf\t%lf\n",i*0.1-5,s[i]);
	fclose(fp);	
	for(i=0;i<8;i++)//求切比雪夫点 
		q[i]=qiebixuefu(-5,5,i,7);
	for (i=0;i<101;i++)//切比雪夫拉格朗日插值 
		s[i]=lagrange(i*0.1-5,q);
	fp=fopen("C:\\qiebixuefu.txt","w+");
	for (i=0;i<101;i++)
		fprintf(fp,"%lf\t%lf\n",i*0.1-5,s[i]);
	fclose(fp);
	/* 三次样条插值*/ 
	double arcs[N][N]={0};//要求逆的矩阵 
	double astar[N][N]={0};//保存伴随矩阵 
	double m[N]={0}; 
	int j;
	spline(arcs) ;
	double a = getA(arcs,N);//开始求逆矩阵，将结果保存至inverse	
	if(abs(a-0)<eps)
		printf("can not transform!\n");
	else
		getAStart(arcs,N,astar);
	double inverse[N][N];
	for(i=0;i<N;i++)
		for(j=0;j<N;j++)
			 inverse[i][j]=astar[i][j]/a;
	mtrxmlt(inverse,g); 
	makes(ss,ans);
	double ii=0;
	for(i=0;i<101;i++)
	{
		s[i]=0;
		ii=-5+i*0.1;
		for(j=0;j<N-1;j++){
			if(a1[j]<=ii&&ii<=a1[j+1]){
				s[i]=ss[j][0]*pow(ii,3)+ss[j][1]*pow(ii,2)+ss[j][2]*ii+ss[j][3];
				break;
			}
		}
			
	}	
	fp=fopen("C:\\lab1.txt","w+");
	for(i=0;i<101;i++)
		fprintf(fp,"%lf\t%lf\n",-5+0.1*i,s[i]);
	fclose(fp);	
}
