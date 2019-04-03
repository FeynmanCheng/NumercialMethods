#include<stdio.h>
#include<time.h>
#include<math.h>
double straight(double x){
    double a;
    a=0.00125*pow(x,5)+0.0625*pow(x,4)+0.425*pow(x,3)+1.2215*pow(x,2)+1.912*x+2.196;
    return a;
}
double qin(double x){
    double a;
    a=((((0.0125*x+0.0625)*x+0.425)*x+1.2215)*x+1.912)*x+2.196;
    return a;
}
int main(){
    int i=0;
    clock_t start,finish;
    double m,n;
    start=clock();
    for(i;i<1000000;i++)
        m=straight(2.0);
    finish=clock();
    printf("Straight time: %f",(double)(finish-start)/CLOCKS_PER_SEC);
    printf("\n");
    start=clock();
    for(i=0;i<1000000;i++)
        n=qin(2.0);
    finish=clock();
    printf("Qin time: %f",(double)(finish-start)/CLOCKS_PER_SEC);
    printf("\n");
}

