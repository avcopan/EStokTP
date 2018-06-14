double sqr(double x);
void VectorNormalisation(double *v);
double Norm(double *v);
double DotProduct(double *v, double *w);
void CrossProduct(double *v, double *w, double *cross);

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

int dihedral_(void) {

  double PI=3.14159;
  double angle;
  double dangle;
  int i;
  FILE *fp;
    
  if((fp=fopen("dihed.dat","r"))==NULL) {
    printf("********FAILED OPENING FILE %s*************\n","dihed.dat");
    exit(1);
  }

  double xc[4];
  double yc[4];
  double zc[4];

  for(i=0;i<4;i++) {
    fscanf(fp,"%lf %lf %lf",&xc[i],&yc[i],&zc[i]);
  } 

  fclose(fp);
  
  double a_b[3];
  double b_c[3];
  double c_d[3];

  a_b[0] = xc[0] - xc[1];
  a_b[1] = yc[0] - yc[1];
  a_b[2] = zc[0] - zc[1];

  b_c[0] = xc[1] - xc[2];
  b_c[1] = yc[1] - yc[2];
  b_c[2] = zc[1] - zc[2];

  c_d[0] = xc[2] - xc[3];
  c_d[1] = yc[2] - yc[3];
  c_d[2] = zc[2] - zc[3];

  double axb;
  double anorm;
  double bnorm;

  axb= DotProduct(a_b,b_c);
  anorm=Norm(a_b);
  bnorm=Norm(b_c);

  angle=180.0 / PI* (PI-acos(axb/anorm/bnorm));

  
  double n1[3];
  double n2[3];
  double m[3];
  double x, y;

  VectorNormalisation(b_c);
  VectorNormalisation(a_b);
  VectorNormalisation(c_d);

  CrossProduct(a_b, b_c, n1);
  CrossProduct(b_c, c_d, n2);
  CrossProduct(b_c,n1, m);

  x = DotProduct(n1, n2);
  y = DotProduct(m, n2);
  
  dangle = 180.0 / PI * atan2(y, x);

  //  printf(" angle is %.4le \n",angle);
  //  printf(" dangle is %.4le \n",dangle);

  if((fp=fopen("dihed.res","w"))==NULL) {
    printf("********FAILED OPENING FILE %s*************\n","dihed.res");
    exit(1);
  }

  fprintf(fp,"%lf \n",angle);
  fprintf(fp,"%lf \n",dangle);
  fclose(fp);

  return 0;
}

double sqr(double x){ return x*x; }

void VectorNormalisation(double *v) { double lenght = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2])); v[0] /= lenght; v[1] /= lenght; v[2] /= lenght; }

double Norm(double *v) { double lenght = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2])); return lenght; }

double DotProduct(double *v, double *w) { return (v[0] * w[0] + v[1] * w[1] + v[2] * w[2]); }

void CrossProduct(double *v, double *w, double *cross) { //
  cross[0] = w[1] * v[2] - w[2] * v[1];
  cross[1] = w[2] * v[0] - w[0] * v[2];
  cross[2] = w[0] * v[1] - w[1] * v[0];  
}


/*

Given the coordinates of the four points, obtain the vectors b1, b2, and b3 by vector subtraction.

Compute n1 and n2, the normal vectors to the planes containing b1 and b2, and b2 and b3 respectively. The angle we seek is the same as the angle between n1 and n2.

The three vectors n1 ,n2, and m1= form an orthonormal frame. Compute the coordinates of n2 in this frame: x= and y=. (You don't need to compute  as it should always be zero.)

The dihedral angle, with the correct sign, is atan2(y,x).

(The reason I recommend the two-argument atan2 function to the traditional cos in this case is both because it naturally produces an angle over a range of , and because cos is poorly conditioned when the angle is close to 0 or )

 */
