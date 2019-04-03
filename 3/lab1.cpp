#include<stdio.h>
#include<math.h>
#include<time.h>
#define eps 1e-6
#define PI 3.14159265
#define N 3
double A[3][3] = {1,1,0.5,1,1,0.25,0.5,0.25,2};
bool judge(double a[N][N]);//判断方阵是否非对角线元素为0 
double max(double a[]);//返回向量中最大元素的值 
double makeP(double a[3][3],double b[3][3]);//构造Jacobi法中的P 
double detA(double arcs[N][N],int n);//求方阵的行列式值 
void mtrmlt(double a[N][N],double b[],double c[]);//方阵乘向量 ，返回至c 
double vecmlt(double *a,double *b,int n);//向量内积 
void division(double a[],double n,double result[]);//幂法中对向量进行标准化 ，返回至result 
void company(double arcs[N][N],int n,double ans[N][N]);//伴随矩阵 ,返回至ans 
void sqrmlt(double a[N][N],double b[N][N],double c[N][N]);//将乘积放入C矩阵
void clear(double a[N][N]);//将＜eps的项置零 
void makeQ(double a[N][N],double q[N][N]);//构造Q矩阵,返回至q 
void transe(double a[N][N],double b[N][N]);//转置 
int main(){
	/*幂法*/ 
	double v[3] = {-0.5315,-0.4615,-0.7103};
	double ans[3] = {0};
	double P[3][3] = {0},Pni[3][3] = {0},detP = 0,inter[N][N] = {0},end[N][N] = {0};
	double max1 = 1,max2 = 0;
	int i,j,k,m,n;
	while(1){
		mtrmlt(A,v,ans);
		max1 = max(ans);
		if(fabs(max1-max2) < eps)
			break;
		division(ans,max1,v);
		max2 = max1;
	}
	printf("lamda = %lf\n",max1);	
	division(ans,max1,v);
	printf("vector is %lf %lf %lf\n",v[0],v[1],v[2]);
	printf("\n");
	/*Jacobi*/
	double multiply = 1;
	while(1)
	{
		makeP(A,P);
		company(P,N,Pni);
		detP = detA(P,N);
		for(i = 0;i < N;i ++)
		{
			for(j = 0;j < N;j++)
			{
				Pni[i][j] =  Pni[i][j] / detP;
			}
		}
		sqrmlt(Pni,A,inter);
		sqrmlt(inter,P,A);
		clear(A);
		if(judge(A))
		{
			break;
		}
	}	
	printf("Jacobi lamda is\n");
	for(i = 0;i < N;i++)
	{
		printf("%lf\t",A[i][i]);
	}
	printf("%\n");
	/*QR*/
	double A[3][3] = {1,1,0.5,1,1,0.25,0.5,0.25,2};
	double Q[N][N],R[N][N],Qt[N][N],AA[N][N];
	makeQ(A,Q); 
	transe(Q,Qt);
	sqrmlt(Qt,A,R);
	sqrmlt(R,Q,AA);
	while(1)
	{
		makeQ(AA,Q); 
		transe(Q,Qt);
		sqrmlt(Qt,AA,R);
		sqrmlt(R,Q,AA);
		clear(AA);
		if(judge(AA))
		{
			break;
		}
	}	 
	printf("\n");
	printf("QR lamda is:\n");
	for(i = 0;i < N;i++)
	{
		printf("%lf\t",AA[i][i]);
	}
}
void transe(double a[N][N],double b[N][N])
{
	int i,j;
	for(i = 0;i < N - 1;i++)
	{
		for(j = i + 1;j < N;j++)
		{
			b[i][j] = a[j][i];
			b[j][i] = a[i][j];
		} 
	}
	for(i = 0;i < N;i++)
	{
		b[i][i] = a[i][i];
	}
}
double vecmlt(double *a,double *b,int n)
{
	int i = 0;
	double result = 0;
	for(;i < n;i++)
	{
		result += *(a + i) * *(b + i);
	}
	return result;
}
void makeQ(double a[N][N],double q[N][N])
{
	double a1[N],a2[N],a3[N];
	int i;
	for(i = 0;i < N;i++)
	{
		a1[i] = a[i][0];
		a2[i] = a[i][1];
		a3[i] = a[i][2];
	}
	double para2 = -vecmlt(a1,a2,N)/vecmlt(a1,a1,N);
	for(i = 0;i < N;i++)
	{
		a2[i] = a[i][1] + para2 * a1[i];
	}
	double para3 = -vecmlt(a1,a3,N)/vecmlt(a1,a1,N);
	para2 = -vecmlt(a2,a3,N)/vecmlt(a2,a2,N);
	for(i = 0;i < N;i++)
	{
		a3[i] = a[i][2] + para3 * a1[i] + para2 * a2[i];
	}
	for(i = 0;i < N;i++)
	{
		q[i][0] = a1[i]/sqrt(vecmlt(a1,a1,N));
		q[i][1] = a2[i]/sqrt(vecmlt(a2,a2,N));
		q[i][2] = a3[i]/sqrt(vecmlt(a3,a3,N));
	}
}
bool judge(double a[N][N])
{
	int i = 0,j = 0;
	for(i = 0;i < N - 1;i++)
	{
		for(j = i + 1;j < N;j++)
		{
			if(fabs(a[i][j]) > eps)
			{
				return 0;
			}	
		}

	}
	return 1;
}
void clear(double a[N][N])
{
	int i,j;
	for(i = 0;i < N - 1;i++)
	{
		for(j = i + 1;j < N;j++)
		{
			if(fabs(a[i][j]) <= eps)
			{
				a[i][j] = 0;
				a[j][i] = 0;
			}
		
		}
	}
}
void mtrmlt(double a[N][N],double b[],double c[])
{
	int i = 0,j = 0;
	for(j = 0;j<3;j++)
		c[j] = 0;
	for(j = 0;j<3;j++)
		for(i = 0;i<3;i++)
			c[j] += a[j][i] * b[i];	
}
double max(double a[])
{
	int i = 1;
	double max = a[0];
	while(i < N){
		if (max < fabs(a[i]))
			max = a[i];
		i++;
	}
	return max;
}
void division(double a[],double n,double result[])
{
	int i = 0;
	while(i < 3){
		result[i] = a[i] / n;
		i++;
	}
}
double makeP(double a[N][N],double b[N][N])
{
	int i = 0,j = 0;
	int l = 0,c = 0;
	double max = 0;
	for(;i < N;i++)
		for(j = 0;j < N;j++)
			b[i][j] = 0;
	for(i = 0;i < N;i++){
		for(j = i + 1;j < N;j++){
			if(fabs(a[i][j]) > max){
				l = i;
				c = j;
				max = fabs(a[i][j]);
			}			
		}
	}
	double theta = 0;
	if(a[l][l] == a[c][c]){
		theta = PI / 4;
	}
	else{
		theta = atan(2 * a[l][c] / (a[l][l] - a[c][c])) / 2;
	}
	b[l][l] = cos(theta);
	b[c][c] = cos(theta);
	b[l][c] = -sin(theta);
	b[c][l] = sin(theta);
	for(i = 0;i < N;i++)
	{
		if(b[i][i] == 0)
		{
			b[i][i] = 1;
		}
	}
}
double detA(double arcs[N][N],int n)
{
	if(n==1)
	{
		return arcs[0][0];
	}
	double ans = 0;
	double temp[N][N];
	int i,j,k;
	for(i=0;i<n;i++)
	{
		for(j=0;j<n-1;j++)
		{
			for(k=0;k<n-1;k++)
			{
				temp[j][k] = arcs[j+1][(k>=i)?k+1:k];
				
			}
		}
		double t = detA(temp,n-1);
		if(i%2==0)
		{
			ans += arcs[0][i]*t;
		}
		else
		{
			ans -=  arcs[0][i]*t;
		}
	}
	return ans;
}
void company(double arcs[N][N],int n,double ans[N][N])
{
	if(n==1)
	{
		ans[0][0] = 1;
		return;
	}
	int i,j,k,t;
	double temp[N][N];
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			for(k=0;k<n-1;k++)
			{
				for(t=0;t<n-1;t++)
				{
					temp[k][t] = arcs[k>=i?k+1:k][t>=j?t+1:t];
				}
			}
			ans[j][i]  =  detA(temp,n-1);
			if((i+j)%2 == 1)
			{
				ans[j][i] = - ans[j][i];
			}
		}
	}
}
void sqrmlt(double a[N][N],double b[N][N],double c[N][N])
{
	int i,j,k;
	for(i = 0;i < N;i++)
	{
		for(j = 0;j < N;j++)
		{
			c[i][j] = 0;
			for(k = 0;k < N;k++)
			{
				c[i][j] += a[i][k] * b[k][j];
			}
		}
	}
}

