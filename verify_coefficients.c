#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
/* Regression check: compare saved iterative coefficient CSVs, never solve LU. */
int main(int argc,char **argv){
 if(argc!=3)return 2;
 const char *suffix[]={"_velocity.csv","_pressure.csv"};
 for(int field=0;field<2;field++){
  char a[2048],b[2048],la[2048],lb[2048];
  snprintf(a,sizeof(a),"%s%s",argv[1],suffix[field]);snprintf(b,sizeof(b),"%s%s",argv[2],suffix[field]);
  FILE *fa=fopen(a,"r"),*fb=fopen(b,"r");if(!fa||!fb)return 3;
  if(!fgets(la,sizeof(la),fa)||!fgets(lb,sizeof(lb),fb))return 4;
  double diff=0,ref=0;int rows=0;
  while(fgets(la,sizeof(la),fa)){
   if(!fgets(lb,sizeof(lb),fb))return 5;
   char *pa=la,*pb=lb,*ea,*eb;
   for(int j=0;j<(field?4:6);j++){
    double x=strtod(pa,&ea),y=strtod(pb,&eb);
    if(ea==pa||eb==pb||!isfinite(x)||!isfinite(y))return 6;
    if(j<3 && fabs(x-y)>1e-12)return 7;
    if(j>=3){diff+=(x-y)*(x-y);ref+=y*y;}
    pa=ea+1;pb=eb+1;
   }rows++;
  }
  if(fgets(lb,sizeof(lb),fb))return 8;
  fclose(fa);fclose(fb);
  double relative=sqrt(diff)/fmax(sqrt(ref),1e-30);
  printf("%s rows=%d relative coefficient difference=%.4e\n",suffix[field],rows,relative);
  if(relative>1e-10)return 9;
 }return 0;
}
