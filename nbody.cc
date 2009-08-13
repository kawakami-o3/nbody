#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
using namespace std;

#define DTime (0.01)
#define TimeStep (1000)
//#define NUM (1024)
//#define NUM (2048)
//#define NUM (4096)
#define NUM (3)
//#define NUM (512)


int swap(int numObj, double r[NUM][3], double v[NUM][3], double m[NUM])
{
  // swap により numObj にあるデータは末尾と交換される
  double tmp;
  for (int i=0 ; i<3 ; i++) {
    tmp = r[NUM-1][i];
    r[NUM-1][i] = r[numObj][i];
    r[numObj][i] = tmp;

    tmp = v[NUM-1][i];
    v[NUM-1][i] = v[numObj][i];
    v[numObj][i] = tmp;
  }
  tmp=m[NUM-1];
  m[NUM-1] = m[numObj];
  m[numObj] = tmp;
  return 0;
}



int printer(int stage, double r[NUM][3], double v[NUM][3], double m[NUM]) {
  ofstream fout;
  int number=NUM;
  //int timestep=TimeStep;
  char chr[256];
  char format[256]="../dat/st%08i.bin";

  //sprintf(chr,"./dat/st%08i.bin",stage);
  sprintf(chr,format,stage);

  fout.open(chr,ios::out|ios::binary|ios::trunc);
  fout.write((char*) &number, sizeof(int));
  //fout.write((char*) &timestep, sizeof(int));
  fout.write((char*) &stage, sizeof(int));


  for (int i=0 ; i<NUM ; i++) {
    fout.write((char*) r[i], sizeof(double)*3);
    fout.write((char*) v[i], sizeof(double)*3);
    fout.write((char*) &m[i], sizeof(double));
  }

  fout.close();
/*
  cout << "log" << endl; 
  for (int i=0 ; i<NUM ; i++) {
    for (int j=0 ; j<3 ; j++)
      cout << r[i][j] << "  ";
    cout << endl;
  }
*/

  /*
  for (int i=0 ; i<NUM ; i++) {
//   printf("%12.5lf %12.5lf %12.5lf   ",r[i][0],r[i][1],r[i][2]);
   printf("%17.10e %17.10e ",r[i][0],r[i][1]);
  }
  printf("\n");
  */
  return 0;
}



int rad(double r[NUM][3], double length[NUM-1])
{
  // 末尾の要素との距離を求め、r^2 を設定する。
  for (int j=0 ; j<NUM-1 ; j++) {
    double ret=0, tmp;
    for (int i=0 ; i<3 ; i++) {
      ret += (r[j][i]-r[NUM-1][i])*(r[j][i]-r[NUM-1][i]);
    }
    length[j] = ret;
  }
  return 0;
}

int force(double f[3], double r[NUM][3], double v[NUM][3], double m[NUM]) {
  double length[NUM-1];
  double soft=1e-1;
  double eps=1e-1;

  for (int i=0 ; i<3 ; i++)
    f[i] = 0.0;
  rad(r,length);
  for (int i=0 ; i<(NUM-1) ; i++) {
    for (int j=0 ; j<3 ; j++) {
      // plummer softenning
      f[j] += (r[i][j]-r[NUM-1][j])/(length[i]+soft)/sqrt(length[i]+soft);
      
      /*
      //spline softenning
      double u2=length[i]/eps/eps;
      double gspline=0;
      if ( u2 >= 4)
        gspline = 1.0/length[i]/sqrt(length[i]);
      else if ( u2 >= 1 )
        gspline = 1.0/length[i]/sqrt(length[i]) * ( -1.0/15 + 8.0/3.0*u2*sqrt(u2)
            - 3*u2*u2 + 6.0/5.0*u2*u2*sqrt(u2) - 1.0/6.0*u2*u2*u2 );
      else
        gspline = 1.0/eps/sqrt(eps) * ( 4.0/3.0 - 6.0/5.0*u2 + 0.5*u2*sqrt(u2));
      f[j] += (r[i][j]-r[NUM-1][j])*gspline;
      */
    }
    /*
    f[0] +=  0.0;
    f[1] += -0.1;
    f[2] +=  0.0;
    */
  }
 
  /*
  for (int j=0 ; j<3 ; j++) {
    f[j] -= v[NUM-1][j]*0.1;
  }
  */

#ifdef DEBUG
  puts("force"); 
  printf("%12.7lf %12.7lf %12.7lf\n",f[0],f[1],f[2]);
  for (int i=0 ; i<NUM-1 ; i++)
    printf("%12.7lf ",length[0]);
  puts("");
#endif
  return 0;
}


int leapfrog(double r[NUM][3], double v[NUM][3], double m[NUM]) {
  int i,j;
  double f[3]; //, rtmp[NUM][3], vtmp[NUM][3];
  const double dt=DTime;

  for (i=0 ; i<NUM ; i++) {
    swap(i,r,v,m);
    for (j=0 ; j<3 ; j++)
      r[NUM-1][j] = r[NUM-1][j] + 0.5*dt*v[NUM-1][j];
    swap(i,r,v,m);
  }
  for (i=0 ; i<NUM ; i++) {
    swap(i,r,v,m);
    force(f,r,v,m);
    for (j=0 ; j<3 ; j++)
      v[NUM-1][j] = v[NUM-1][j] + dt*f[j];
    swap(i,r,v,m);
  }
  for (i=0 ; i<NUM ; i++) {
    swap(i,r,v,m);
    for (j=0 ; j<3 ; j++)
      r[NUM-1][j] = r[NUM-1][j] + 0.5*dt*v[NUM-1][j];
    swap(i,r,v,m);
  }
  return 0;
}

int boundary(double r[NUM][3], double v[NUM][3])
{
  const double length=20;
  const double alpha=0.7;
  double tmp;
  for (int i=0 ; i<NUM ; i++) {
    /*
    if (r[i][0] > 20) {
      r[i][0] = 2*length-r[i][0];
      v[i][0] = -v[i][0]*alpha;
    } else if (r[i][0] < 0) {
      r[i][0] = -r[i][0];
      v[i][0] = -v[i][0]*alpha;
    }
    if (r[i][1] < 0) {
      r[i][1] = -r[i][1];
      v[i][1] = -v[i][1]*alpha;
    }

    if (r[i][2] > 20) {
      r[i][2] = 2*length-r[i][2];
      v[i][2] = -v[i][2]*alpha;
    } else if (r[i][2] < 0) {
      r[i][2] = -r[i][2];
      v[i][2] = -v[i][2]*alpha;
    }
    */


    for (int j=0 ; j<3 ; j++) {
      /*
      if (r[i][j] > 20) {
        r[i][j] = 2*length-r[i][j];
        v[i][j] = -v[i][j]*alpha;
      } else if (r[i][j] < 0) {
        r[i][j] = -r[i][j];
        v[i][j] = -v[i][j]*alpha;
      }
      */
      tmp = fmod(r[i][j],length);
      if (tmp > 0)
        r[i][j]=tmp;
      else
        r[i][j]=length+tmp;
    }
  }
  return 0;
}

int main(void)
{
  double r[NUM][3],v[NUM][3],m[NUM];


  r[0][0] =  1.0; r[0][1] =  0.0; r[0][2] = 0.0;
  r[1][0] = -1.5; r[1][1] =  0.0; r[1][2] = 0.0;
  r[2][0] =  2.0; r[2][1] =  2.0; r[2][2] = 0.0;

  v[0][0] =  0.5; v[0][1] =  0.1; v[0][2] = 0.0;
  v[1][0] =  0.2; v[1][1] = -0.5; v[1][2] = 0.0;
  v[2][0] = 00.7; v[2][1] = -0.3; v[2][2] = 0.0;
  
  m[0] = 1.0;
  m[1] = 1.0;
  m[2] = 1.0;


//#include "init.h"
//#include "init-twodisk.h"
//#include "init-singledisk.h"
//  boundary(r,v);
//  for (int j=0 ; j<NUM ; j++) {
//    printf("%i,%i: %e %e, %e %e, %e\n",0,j,r[j][0],r[j][1],v[j][0],v[j][1],m[j]);
//  }




  printer(0,r,v,m);
  for (int i=1 ; i<TimeStep ; i++) {
    printf("%i\n",i);
    leapfrog(r,v,m);
//    boundary(r,v);
    printer(i,r,v,m);
//  for (int j=0 ; j<NUM ; j++) {
//    printf("%i,%i: %e %e, %e %e, %e\n",i,j,r[j][0],r[j][1],v[j][0],v[j][1],m[j]);
//  }
  }
  return 0;
}
