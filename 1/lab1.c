#include<stdio.h>
#include<math.h>
int main(){
	double val=0,ins=1;
	int i;
	for(i=0;i<10;i++){
		printf("%.8f\n",2.718281828-val-ins);
		val=ins+val;
		ins=ins/(i+1);
	}
}
