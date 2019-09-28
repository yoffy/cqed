#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

int Jacobi(int n, int ct, double eps, double **A, double **A1, double **A2,
           double **X1, double **X2) {
  double max, s, t, v, sn, cs1;
  int i1, i2, k = 0, ind = 1, p = 0, q = 0;
  // 初期設定
  for (i1 = 0; i1 < n; i1++) {
    for (i2 = 0; i2 < n; i2++) {
      A1[i1][i2] = A[i1][i2];
      X1[i1][i2] = 0.0;
    }
    X1[i1][i1] = 1.0;
  }
  // 計算
  while (ind > 0 && k < ct) {
    // 最大要素の探索
    max = 0.0;
    for (i1 = 0; i1 < n; i1++) {
      for (i2 = 0; i2 < n; i2++) {
        if (i2 != i1) {
          if (fabs(A1[i1][i2]) > max) {
            max = fabs(A1[i1][i2]);
            p = i1;
            q = i2;
          }
        }
      }
    }
    // 収束判定
    // 収束した
    if (max < eps)
      ind = 0;
    // 収束しない
    else {
      // 準備
      s = -A1[p][q];
      t = 0.5 * (A1[p][p] - A1[q][q]);
      v = fabs(t) / sqrt(s * s + t * t);
      sn = sqrt(0.5 * (1.0 - v));
      if (s * t < 0.0)
        sn = -sn;
      cs1 = sqrt(1.0 - sn * sn);
      // Akの計算
      for (i1 = 0; i1 < n; i1++) {
        if (i1 == p) {
          for (i2 = 0; i2 < n; i2++) {
            if (i2 == p)
              A2[p][p] = A1[p][p] * cs1 * cs1 + A1[q][q] * sn * sn -
                         2.0 * A1[p][q] * sn * cs1;
            else if (i2 == q)
              A2[p][q] = 0.0;
            else
              A2[p][i2] = A1[p][i2] * cs1 - A1[q][i2] * sn;
          }
        } else if (i1 == q) {
          for (i2 = 0; i2 < n; i2++) {
            if (i2 == q)
              A2[q][q] = A1[p][p] * sn * sn + A1[q][q] * cs1 * cs1 +
                         2.0 * A1[p][q] * sn * cs1;
            else if (i2 == p)
              A2[q][p] = 0.0;
            else
              A2[q][i2] = A1[q][i2] * cs1 + A1[p][i2] * sn;
          }
        } else {
          for (i2 = 0; i2 < n; i2++) {
            if (i2 == p)
              A2[i1][p] = A1[i1][p] * cs1 - A1[i1][q] * sn;
            else if (i2 == q)
              A2[i1][q] = A1[i1][q] * cs1 + A1[i1][p] * sn;
            else
              A2[i1][i2] = A1[i1][i2];
          }
        }
      }
      // Xkの計算
      for (i1 = 0; i1 < n; i1++) {
        for (i2 = 0; i2 < n; i2++) {
          if (i2 == p)
            X2[i1][p] = X1[i1][p] * cs1 - X1[i1][q] * sn;
          else if (i2 == q)
            X2[i1][q] = X1[i1][q] * cs1 + X1[i1][p] * sn;
          else
            X2[i1][i2] = X1[i1][i2];
        }
      }
      // 次のステップへ
      k++;
      for (i1 = 0; i1 < n; i1++) {
        for (i2 = 0; i2 < n; i2++) {
          A1[i1][i2] = A2[i1][i2];
          X1[i1][i2] = X2[i1][i2];
        }
      }
    }
  }

  return ind;
}

void LuvtoRGB(double L, double u, double v, double *R, double *G, double *B) {
  double X, Y, Z, ud, vd, u0, v0, TEMP, L1;
  double RR, GG, BB;
  double eps = 216.0 / 24389.0;
  double k = 24389.0 / 27.0;
  double Xr = 0.964221; // reference white D50
  double Yr = 1.0;
  double Zr = 0.825211;

  u0 = 4.0 * Xr / (Xr + 15.0 * Yr + 3.0 * Zr);
  v0 = 9.0 * Yr / (Xr + 15.0 * Yr + 3.0 * Zr);
  L1 = ((double)(L)) / 1.0;
  if ((double)(L1) > k * eps) {
    TEMP = (((double)(L1) + 16.0) / 116.0);
    Y = TEMP * TEMP * TEMP;
  } else {
    Y = ((double)(L1)) / k;
  }
  if ((L == 0) && (u == 0) && (v == 0)) {
    X = 0;
    Y = 0;
    Z = 0;
  } else {
    ud = (u / (13.0 * L1) + u0);
    vd = (v / (13.0 * L1) + v0);
    X = (ud / vd) * Y * 9.0 / 4.0;
    Z = (Y / vd - ((ud / vd) * Y / 4.0 + 15.0 * Y / 9.0)) * 3.0;
  }

  // a=((52.0*L)/(u+13.0*L*u0)-1.0)/3.0;
  // b=-5.0*Y;
  // c=-1.0/3.0;
  // d=Y*( (39.0*((double)(L))/(v+13.0*L*v0))-5.0 );
  // X=(d-b)/(a-c);
  // Z=X*a+b;
  // fY = pow((L + 16.0) / 116.0, 3.0);
  // if (fY < 0.008856)
  // fY = L / 903.3;
  // Y = fY;

  // if (fY > 0.008856)
  // fY = pow(fY, 1.0/3.0);
  // else
  // fY = 7.787 * fY + 16.0/116.0;

  // fX = a / 500.0 + fY;
  // if (fX > 0.206893)
  // X = pow(fX, 3.0);
  // else
  // X = (fX - 16.0/116.0) / 7.787;

  // fZ = fY - b /200.0;
  // if (fZ > 0.206893)
  // Z = pow(fZ, 3.0);
  // else
  // Z = (fZ - 16.0/116.0) / 7.787;

  X *= 255.0; //(0.950456 * 255);
  Y *= 255.0;
  Z *= 255.0; //(1.088754 * 255);

  RR = (3.2404813432005 * X - 1.5371515162713 * Y - 0.49853632616889 * Z);
  GG = (-0.96925494999657 * X + 1.8759900014899 * Y + 0.041555926558293 * Z);
  BB = (0.055646639135177 * X - 0.20404133836651 * Y + 1.0573110696453 * Z);

  *R = RR; // < 0 ? 0 : RR > 255 ? 255 : RR;
  *G = GG; // < 0 ? 0 : GG > 255 ? 255 : GG;
  *B = BB; // < 0 ? 0 : BB > 255 ? 255 : BB;
}

void RGBtoLuv(double R, double G, double B, double *L, double *u, double *v) {

  double X, Y, Z; //, fX, fY, fZ;
  /* double r,g,b;
  r = R/255.f; //R 0..1
  g = G/255.f; //G 0..1
  b = B/255.f; //B 0..1
  if (r <= 0.04045)
  r = r/12;
  else
  r = (double) pow((r+0.055)/1.055,2.4);

  if (g <= 0.04045)
  g = g/12;
  else
  g = (double) pow((g+0.055)/1.055,2.4);

  if (b <= 0.04045)
  b = b/12;
  else
  b = (double) pow((b+0.055)/1.055,2.4);
  */
  X = 0.412453 * R + 0.357580 * G + 0.180423 * B;
  Y = 0.212671 * R + 0.715160 * G + 0.072169 * B;
  Z = 0.019334 * R + 0.119193 * G + 0.950227 * B;

  X /= 255.0; //(255 * 0.950456);
  Y /= 255.0; // 255;
  Z /= 255.0; //(255 * 1.088754);
  // X = 0.412453*R + 0.357580*G + 0.180423*B;
  // Y = 0.212671*R + 0.715160*G + 0.072169*B;
  // Z = 0.019334*R + 0.119193*G + 0.950227*B;

  // X = 0.436052025f*r + 0.385081593f*g + 0.143087414f *b;
  // Y = 0.222491598f*r + 0.71688606f *g + 0.060621486f *b;
  // Z = 0.013929122f*r + 0.097097002f*g + 0.71418547f *b;

  double L1, u1, v1, u2, v2, ur2, vr2, yr, eps, k;
  eps = 216.0 / 24389.0;
  k = 24389.0 / 27.0;

  double Xr = 0.964221; // reference white D50
  double Yr = 1.0;      // 1.0;
  double Zr = 0.825211;

  // u2 = 4.0*X / (X + 15.0*Y + 3.0*Z);
  // v2 = 9.0*Y / (X + 15.0*Y + 3.0*Z);
  u2 = 4.0 * X / (X + 15.0 * Y + 3.0 * Z);
  v2 = 9.0 * Y / (X + 15.0 * Y + 3.0 * Z);

  ur2 = 4.0 * Xr / (Xr + 15.0 * Yr + 3.0 * Zr);
  vr2 = 9.0 * Yr / (Xr + 15.0 * Yr + 3.0 * Zr);

  yr = Y / Yr;

  if (yr > eps) {
    L1 = (116.0 * pow(yr, 1.0 / 3.0) - 16.0);
    // L1 = (double) (116.0 * pow(yr, 0.3333333333333) - 16.0);
    // L1 = (double) (116.0 * yr - 16.0);
    // debug start
    // fprintf(stderr,"YATTER\n");

    // debug end
  } else {
    L1 = k * yr;
  }

  // u1 = 13.0*(L1)*(u2 -ur2);
  // v1 = 13.0*(L1)*(v2 -vr2);
  u1 = 13.0 * (L1) * (u2 - ur2); // -ur2);
  v1 = 13.0 * (L1) * (v2 - vr2); // -vr2);
  // u1 = L1;
  // v1 = L1;
  if (X == 0.0 && Y == 0.0 && Z == 0.0) {
    *L = 0.0;
    *u = 0.0;
    *v = 0.0;
  } else {
    // *L = (int) (2.55*(L1) + 0.5);
    *L = /*(int) ((*/ L1 /*) + 0.5)*/;
    *u = /*(int) (*/ u1 /* + 0.5)*/;
    *v = /*(int) (*/ v1 /* + 0.5)*/;
  }

  // *L = R;
  // *u = G;
  // *v = B;

  // X /= (255 * 0.950456);
  // Y /= 255;
  // Z /= (255 * 1.088754);

  // if (Y > 0.008856)
  // {
  // fY = pow(Y, 1.0/3.0);
  // *L = (int)(116.0*fY - 16.0 + 0.5);
  // }
  // else
  // {
  // fY = 7.787*Y + 16.0/116.0;
  // *L = (int)(903.3*Y + 0.5);
  // }

  // if (X > 0.008856)
  // fX = pow(X, 1.0/3.0);
  // else
  // fX = 7.787*X + 16.0/116.0;

  // if (Z > 0.008856)
  // /fZ = pow(Z, 1.0/3.0);
  // else
  // fZ = 7.787*Z + 16.0/116.0;

  // *a = (int)(500.0*(fX - fY) + 0.5);
  // *b = (int)(200.0*(fY - fZ) + 0.5);
}

void LABtoRGB(double L, double a, double b, double *R, double *G, double *B) {
  // Convert between RGB and CIE-Lab color spaces
  // Uses ITU-R recommendation BT.709 with D65 as reference white.
  // algorithm contributed by "Mark A. Ruzon" <ruzon@CS.Stanford.EDU>

  double X, Y, Z, fX, fY, fZ;
  double RR, GG, BB;
  fY = (L + 16.0) / 116.0;
  fX = fY + a / 500.0;
  fZ = fY - b / 200.0;
  if (fY > pow(0.008856, 1.0 / 3.0)) {
    Y = pow(fY, 3.0);
  } else {
    Y = 27.0 / 29.0 / 29.0 / 29.0 * (116.0 * fY - 16.0);
  }
  if (fX > pow(0.008856, 1.0 / 3.0)) {
    X = pow(fX, 3.0);
  } else {
    X = 27.0 / 29.0 / 29.0 / 29.0 * (116.0 * fX - 16.0);
  }
  if (fZ > pow(0.008856, 1.0 / 3.0)) {
    Z = pow(fZ, 3.0);
  } else {
    Z = 27.0 / 29.0 / 29.0 / 29.0 * (116.0 * fZ - 16.0);
  }
  X *= (0.950456 * 255.0);
  Y *= 255;
  Z *= (1.088754 * 255.0);

  RR = (3.2408109640905 * X - 1.5373099569082 * Y - 0.49858604825854 * Z);
  GG = (-0.96924116992192 * X + 1.8759664973723 * Y + 0.041553922456843 * Z);
  BB = (0.05563752111349 * X - 0.20400735547233 * Y + 1.0571298579506 * Z);

  *R = RR; // < 0 ? 0 : RR > 255 ? 255 : RR;
  *G = GG; // < 0 ? 0 : GG > 255 ? 255 : GG;
  *B = BB; // < 0 ? 0 : BB > 255 ? 255 : BB;
}

void RGBtoLAB(double R, double G, double B, double *L, double *a, double *b) {
  // Convert between RGB and CIE-Lab color spaces
  // Uses ITU-R recommendation BT.709 with D65 as reference white.
  // algorithm contributed by "Mark A. Ruzon" <ruzon@CS.Stanford.EDU>

  double X, Y, Z, tempX, tempY, tempZ;
  /*
  X = 0.412453*R + 0.357580*G + 0.180423*B;
  Y = 0.212671*R + 0.715160*G + 0.072169*B;
  Z = 0.019334*R + 0.119193*G + 0.950227*B;
  */
  X = 0.412411 * R + 0.357585 * G + 0.180454 * B;
  Y = 0.212649 * R + 0.715169 * G + 0.072182 * B;
  Z = 0.019332 * R + 0.119195 * G + 0.950390 * B;
  X /= (255.0 * 0.950456);
  Y /= 255.0;
  Z /= (255.0 * 1.088754);

  if (Y > 0.008856) {
    tempY = pow(Y, 1.0 / 3.0);
  } else {
    tempY = (29.0 * 29.0 * 29.0 / 27.0 * Y + 16.0) / 116.0;
  }

  if (X > 0.008856) {
    tempX = pow(X, 1.0 / 3.0);
  } else {
    tempX = (29.0 * 29.0 * 29.0 / 27.0 * X + 16.0) / 116.0;
  }

  if (Z > 0.008856) {
    tempZ = pow(Z, 1.0 / 3.0);
  } else
    tempZ = (29.0 * 29.0 * 29.0 / 27.0 * Z + 16.0) / 116.0;

  *a = ((500.0 * (tempX - tempY)));
  *b = ((200.0 * (tempY - tempZ)));
  *L = ((116.0 * tempY - 16.0));
}

template <typename T>
std::vector<T> Histogram(const T *X, int NUM, T MIN, int NODEHANI) {
  std::vector<int> HIST(NODEHANI);
  int TMP;
  for (int i = 0; i < NUM; i++) {
    TMP = *(X + i) - MIN;
    HIST[TMP]++;
  }
  return HIST;
}

int ohtsu(int NUM, const int *X) {
  std::pair<const int *, const int *> MINMAX =
      std::minmax_element(&X[0], &X[NUM]);
  int MAX = *MINMAX.second;
  int MIN = *MINMAX.first;
  int NODEHANI = MAX - MIN + 1;
  std::vector<int> HIST = Histogram(X, NUM, MIN, NODEHANI);

  float MAX2 = 0.0;
  float AVE1, AVE2;
  int TOTAL1, TOTAL2;
  float HANTEI;
  int THRESH;
  for (int i = 1; i < NODEHANI; i++) {
    AVE1 = 0.0;
    TOTAL1 = 0;
    for (int j = 0; j < i; j++) {
      AVE1 += HIST[j] * (j + MIN);
      TOTAL1 += HIST[j];
    }
    AVE1 /= (float)(TOTAL1);
    AVE2 = 0.0;
    TOTAL2 = 0;
    for (int j = i; j < NODEHANI; j++) {
      AVE2 += HIST[j] * (j + MIN);
      TOTAL2 += HIST[j];
    }
    AVE2 /= (float)(TOTAL2);

    HANTEI = (float)(TOTAL1) * (float)(TOTAL2) * (AVE1 - AVE2) * (AVE1 - AVE2);

    if (MAX2 < HANTEI) {
      MAX2 = HANTEI;
      THRESH = i;
    }
  }
  THRESH += MIN;

  return (THRESH);
} // main kansuu end

int ohtsu2(int NUM, const double *X, const double *Y, const double *Z,
           int omh) {
  int i, j;
  // NUM=7;
  // omh=3;

  // double XX[7] = {-10.5,-9.0,-8.0,-7.0,-6.0,-5.5,-3.2};
  // double XX[7] = {-10.0,-7.5,-6.8,-6.5,-6.3,-5.5,-3.2};
  //    double XX[7] = {-9.0,-7.5,-6.8,-6.5,-6.3,-5.5,-3.2};
  // double XX[7] = {-10.0,-7.5,-6.8,-6.5,-6.3,-5.5,5.3};
  // double XX[7] = {-10.0,-7.5,-6.8,-6.5,-6.3,-5.5,5.0};
  // double XX[7] = {10.5,13.0,14.8,16.5,18.3,19.5,25.0};
  // double XX[7] = {10.0,13.0,14.8,16.5,18.3,19.5,25.0};
  // double XX[7] = {0.0,2.0,14.8,16.5,18.3,19.5,25.0};
  /// double XX[7] = {-10.5,-9.0,-8.0,-7.0,-6.0,-5.5,0.0};

  //    for(i=0;i<7;i++){
  //        X[i] = XX[i];
  //    }

  // -10,-10,-10,-10,-10,
  // -3,-3,-3,-3,-3,
  // -1,-1,0,1,1,
  // 3,3,3,3,3,
  // 10,10,10,10,10,
  // 10,10,10,10,10};
  // printf("OHTSU2 START!!\n");
  std::pair<const double *, const double *> MINMAX =
      std::minmax_element(&X[0], &X[NUM]);
  double MAX = *MINMAX.second;
  double MIN = *MINMAX.first;
  int CMAX, CMIN;

  if (MAX > 0.0) {
    CMAX = (int)MAX;
  } else if (MAX == 0.0) {
    CMAX = 0.0;
  } else {
    CMAX = (int)MAX - 1;
  }
  if (MIN > 0.0) {
    CMIN = (int)MIN;
  } else if (MIN == 0.0) {
    CMIN = 0.0;
  } else {
    if (MIN - (double)((int)(MIN)) != 0.0) {
      CMIN = (int)MIN - 1;
    } else {
      CMIN = (int)MIN;
    }
  }

  int NODEHANI = CMAX - CMIN + 1;
  // int *HIST;
  // HIST = (int *)malloc(sizeof(int)*NODEHANI);
  // ND
  // 0- -11
  // 1- -10
  //   -6.7 - 5.3
  //    -7 - 6  14
  // 11- 0

  // 22- 11
  // for(i=0;i<NODEHANI;i++){
  //*(HIST+i) = 0;
  //}
  // int TMP;
  // for(i=0;i<NUM;i++){
  // TMP = *(X+i)-MIN;
  //(*(HIST+TMP))++;
  //}
  // debug start
  // for(i=0;i<NODEHANI;i++){
  // fprintf(stderr,"NUM ATAI HIST %d %d %d\n",i,i+MIN,*(HIST+i));
  //}
  // while(1);
  // debug end

  double MAX2 = 0.0;
  double MIN2 = 99999.99e64;
  double AVE1, AVE2, BUN1, BUN2, AVEY1, AVEY2, BUNY1, BUNY2, AVEZ1, AVEZ2,
      BUNZ1, BUNZ2;
  int TOTAL1, TOTAL2;
  double HANTEI, HANTEI2;
  int THRESH, THRESH2;
  // printf("CMIN=%d,CMAX=%d,NODEHANI=%d\n",CMIN,CMAX,NODEHANI);
  for (i = 1; i < NODEHANI; i++) {
    // printf("i=%d\n",i);
    AVE1 = 0.0;
    TOTAL1 = 0;
    BUN1 = 0.0;
    AVE2 = 0.0;
    TOTAL2 = 0;
    BUN2 = 0.0;
    AVEY1 = 0.0;
    BUNY1 = 0.0;
    AVEY2 = 0.0;
    BUNY2 = 0.0;
    AVEZ1 = 0.0;
    BUNZ1 = 0.0;
    AVEZ2 = 0.0;
    BUNZ2 = 0.0;

    // for(j=0;j<i;j++){
    for (j = 0; j < NUM; j++) {
      if (X[j] - (double)(CMIN) < (double)i) {
        AVE1 += X[j];
        AVEY1 += Y[j];
        AVEZ1 += Z[j];
        TOTAL1++;
      } else {
        AVE2 += X[j];
        AVEY2 += Y[j];
        AVEZ2 += Z[j];
        TOTAL2++;
      }
    }
    //    printf("TOTAL1 = %d,TOTAL2=%d\n",TOTAL1,TOTAL2);
    AVE1 /= (double)(TOTAL1);
    AVE2 /= (double)(TOTAL2);
    AVEY1 /= (double)(TOTAL1);
    AVEY2 /= (double)(TOTAL2);
    AVEZ1 /= (double)(TOTAL1);
    AVEZ2 /= (double)(TOTAL2);
    //    printf("TOTAL1 = %d,TOTAL2 = %d\n",TOTAL1,TOTAL2);
    // for(j=i;j<NODEHANI;j++){
    // AVE2 += (*(HIST+j))*(j+MIN);
    // TOTAL2 += (*(HIST+j));
    //}
    // AVE2 /= (float)(TOTAL2);
    for (j = 0; j < NUM; j++) {
      if (X[j] - (double)(CMIN) < (double)i) {
        BUN1 += (X[j] - AVE1) * (X[j] - AVE1);
        BUNY1 += (Y[j] - AVEY1) * (Y[j] - AVEY1);
        BUNZ1 += (Z[j] - AVEZ1) * (Z[j] - AVEZ1);
      } else {
        BUN2 += (X[j] - AVE2) * (X[j] - AVE2);
        BUNY2 += (Y[j] - AVEY2) * (Y[j] - AVEY2);
        BUNZ2 += (Z[j] - AVEZ2) * (Z[j] - AVEZ2);
      }
    }
    BUN1 /= TOTAL1;
    BUN2 /= TOTAL2;
    BUNY1 /= TOTAL1;
    BUNY2 /= TOTAL2;
    BUNZ1 /= TOTAL1;
    BUNZ2 /= TOTAL2;
    if (omh == 3) {
      HANTEI =
          (double)(TOTAL1) * (double)(TOTAL2) * (AVE1 - AVE2) * (AVE1 - AVE2);
    }
    if (omh == 4) {
      HANTEI2 = (BUN1 + BUNY1 + BUNZ1) * (double)TOTAL1 +
                (BUN2 + BUNY2 + BUNZ2) * (double)TOTAL2;
    }
    //    printf("HANTEI=%f\n",HANTEI);
    if (omh == 3) {
      if (MAX2 < HANTEI) {
        MAX2 = HANTEI;
        THRESH = i;
      }
    }
    if (omh == 4) {
      if (MIN2 > HANTEI2) {
        MIN2 = HANTEI2;
        THRESH2 = i;
      }
    }
  }
  THRESH += CMIN;
  THRESH2 += CMIN;
  // debug start
  fprintf(stderr, "CMIN=%d,THRESH=%d,THRESH2=%d,CMAX=%d\n", CMIN, THRESH,
          THRESH2, CMAX);
  //        while(1);
  // debug end
  //    printf("ohtsu2 end!!\n");
  if (omh == 3) {
    return (THRESH);
  } else if (omh == 4) {
    return (THRESH2);
  }
  return (0);
} // main kansuu end

int media(int NUM, const int *X) {
  int i, j;
  // int SUM1;
  // int X[NUM] = {-10,-10,-10,-10,-10,
  // -10,-10,-10,-10,-10,
  // -3,-3,-3,-3,-3,
  // -1,-1,0,1,1,
  // 3,3,3,3,3,
  // 10,10,10,10,10,
  // 10,10,10,10,10};
  std::pair<const int *, const int *> MINMAX =
      std::minmax_element(&X[0], &X[NUM]);
  int MAX = *MINMAX.second;
  int MIN = *MINMAX.first;
  int NODEHANI = MAX - MIN + 1;
  std::vector<int> HIST = Histogram(X, NUM, MIN, NODEHANI);
  int GOUKEI, MED;
  GOUKEI = 0;
  for (i = 0; i < NODEHANI; i++) {
    GOUKEI += HIST[i];
    if (NUM / 2 < GOUKEI) {
      MED = i - 1;
      break;
    }
  }
  if (MED == -1) {
    MED = MIN;
  } else {
    MED += MIN;
  }

  float MAX2 = 0.0;
  float AVE1, AVE2;
  int TOTAL1, TOTAL2;
  float HANTEI;
  int THRESH;
  goto klk;
  for (i = 1; i < NODEHANI; i++) {
    AVE1 = 0.0;
    TOTAL1 = 0;
    for (j = 0; j < i; j++) {
      AVE1 += HIST[j] * (j + MIN);
      TOTAL1 += HIST[j];
    }
    AVE1 /= (float)(TOTAL1);
    AVE2 = 0.0;
    TOTAL2 = 0;
    for (j = i; j < NODEHANI; j++) {
      AVE2 += HIST[j] * (j + MIN);
      TOTAL2 += HIST[j];
    }
    AVE2 /= (float)(TOTAL2);

    HANTEI = (float)(TOTAL1) * (float)(TOTAL2) * (AVE1 - AVE2) * (AVE1 - AVE2);

    if (MAX2 < HANTEI) {
      MAX2 = HANTEI;
      THRESH = i;
    }
  }
  THRESH += MIN;
// debug start
// fprintf(stderr,"THRE=%d",THRESH);
// debug end
klk:
  THRESH = MED;
  printf("MIN=%d,MED=%d,MAX=%d\n", MIN, THRESH, MAX); // while(1);
  return (THRESH);
  // return (0);
} // main kansuu end

void modoshi(double V[3], double R, double G, double B, double *X, double *Y,
             double *Z) {
  double COSTHIETA, SINTHIETA, SQV1V2, ONEMINCOSTHIETA, NX, NY, NZ;

  SINTHIETA = -sqrt((V[1] * V[1] + V[2] * V[2]) /
                    (V[0] * V[0] + V[1] * V[1] + V[2] * V[2]));
  COSTHIETA = V[0] / sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
  SQV1V2 = sqrt(V[1] * V[1] + V[2] * V[2]);
  NX = 0.0;
  NY = V[2] / SQV1V2;
  NZ = -V[1] / SQV1V2;
  ONEMINCOSTHIETA = 1.0 - COSTHIETA;
  *X = (NX * NX * ONEMINCOSTHIETA + COSTHIETA) * R +
       (NX * NY * ONEMINCOSTHIETA - NZ * SINTHIETA) * G +
       (NZ * NX * ONEMINCOSTHIETA + NY * SINTHIETA) * B;
  *Y = (NX * NY * ONEMINCOSTHIETA + NZ * SINTHIETA) * R +
       (NY * NY * ONEMINCOSTHIETA + COSTHIETA) * G +
       (NY * NZ * ONEMINCOSTHIETA - NX * SINTHIETA) * B;
  *Z = (NZ * NX * ONEMINCOSTHIETA - NY * SINTHIETA) * R +
       (NY * NZ * ONEMINCOSTHIETA + NX * SINTHIETA) * G +
       (NZ * NZ * ONEMINCOSTHIETA + COSTHIETA) * B;
} // fprintf(stderr,"X1=%f,Y1=%f,Z1=%fTHIETA=%fFAI=%f,%f\n",X1,Y1,Z1,THIETA/3.14,FAI/3.14,sqrt(V[0]*V[0]+V[1]*V[1])/V[2]);

void kaiten(double V[3], double X, double Y, double Z, double *R, double *G,
            double *B) {
  double SINTHIETA, COSTHIETA;
  double NX, NY, NZ;
  double R1, G1, B1;
  double SQV1V2, ONEMINCOSTHIETA;
  SINTHIETA = sqrt((V[1] * V[1] + V[2] * V[2]) /
                   (V[0] * V[0] + V[1] * V[1] + V[2] * V[2]));
  COSTHIETA = V[0] / sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
  SQV1V2 = sqrt(V[1] * V[1] + V[2] * V[2]);
  NX = 0.0;
  NY = V[2] / SQV1V2;
  NZ = -V[1] / SQV1V2;
  ONEMINCOSTHIETA = 1.0 - COSTHIETA;
  R1 = (NX * NX * ONEMINCOSTHIETA + COSTHIETA) * X +
       (NX * NY * ONEMINCOSTHIETA - NZ * SINTHIETA) * Y +
       (NZ * NX * ONEMINCOSTHIETA + NY * SINTHIETA) * Z;
  G1 = (NX * NY * ONEMINCOSTHIETA + NZ * SINTHIETA) * X +
       (NY * NY * ONEMINCOSTHIETA + COSTHIETA) * Y +
       (NY * NZ * ONEMINCOSTHIETA - NX * SINTHIETA) * Z;
  B1 = (NZ * NX * ONEMINCOSTHIETA - NY * SINTHIETA) * X +
       (NY * NZ * ONEMINCOSTHIETA + NX * SINTHIETA) * Y +
       (NZ * NZ * ONEMINCOSTHIETA + COSTHIETA) * Z;

  // *R=R1+0.5;
  // *G=G1+0.5;
  //*B=B1+0.5;
  *R = R1;
  *G = G1;
  *B = B1;
}
