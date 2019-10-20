#include "mediancut_util.h"

#include <algorithm>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

// 二乗
template <typename T> static inline T Square(T v) { return v * v; }

// PIXTYPES[i] != PIXTYPE である画素を使って A に共分散を設定する
void JacobiSetupNotEquals(int SIZE, const int16_t *PIXTYPES, int16_t PIXTYPE,
                          const double *YIN, const double *UIN,
                          const double *VIN,
                          const std::tuple<double, double, double> &U_RGB,
                          double A[3][3]) {
  // Aに分散共分散を入れる
  int NUM_IGNORES = 0;
  double TMP_RR = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(RIN+i))>>RDIV)-U_R);
  double TMP_GG = 0.0; //(((*(GIN+i))>>GDIV)-U_G)*(((*(GIN+i))>>GDIV)-U_G);
  double TMP_BB = 0.0; //(((*(BIN+i))>>BDIV)-U_B)*(((*(BIN+i))>>BDIV)-U_B);
  double TMP_RG = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(GIN+i))>>GDIV)-U_G);
  double TMP_RB = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(BIN+i))>>BDIV)-U_B);
  double TMP_GB = 0.0; //(((*(GIN+i))>>GDIV)-U_G)*(((*(BIN+i))>>BDIV)-U_B);
  double U_R = std::get<0>(U_RGB);
  double U_G = std::get<1>(U_RGB);
  double U_B = std::get<2>(U_RGB);
  for (int i = 0; i < SIZE; i++) {
    if (PIXTYPES[i] == PIXTYPE) {
      NUM_IGNORES++;
      continue;
    }
    TMP_RR += (YIN[i] - U_R) * (YIN[i] - U_R);
    TMP_GG += (UIN[i] - U_G) * (UIN[i] - U_G);
    TMP_BB += (VIN[i] - U_B) * (VIN[i] - U_B);
    TMP_RG += (YIN[i] - U_R) * (UIN[i] - U_G);
    TMP_RB += (YIN[i] - U_R) * (VIN[i] - U_B);
    TMP_GB += (UIN[i] - U_G) * (VIN[i] - U_B);
  }
  SIZE -= NUM_IGNORES;
  TMP_RR /= SIZE;
  TMP_GG /= SIZE;
  TMP_BB /= SIZE;
  TMP_RG /= SIZE;
  TMP_RB /= SIZE;
  TMP_GB /= SIZE;

  A[0][0] = TMP_RR;
  A[0][1] = TMP_RG;
  A[1][0] = TMP_RG;
  A[1][1] = TMP_GG;
  A[2][2] = TMP_BB;
  A[0][2] = TMP_RB;
  A[2][0] = TMP_RB;
  A[1][2] = TMP_GB;
  A[2][1] = TMP_GB;
}

// PIXTYPES[i] == PIXTYPES である画素を使って A に共分散を設定する
void JacobiSetupEquals(int SIZE, const int16_t *PIXTYPES, int16_t PIXTYPE,
                       const double *YIN, const double *UIN, const double *VIN,
                       const std::tuple<double, double, double> &U_RGB,
                       double A[3][3]) {
  // Aに分散共分散を入れる
  int NUM_IGNORES = 0;
  double TMP_RR = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(RIN+i))>>RDIV)-U_R);
  double TMP_GG = 0.0; //(((*(GIN+i))>>GDIV)-U_G)*(((*(GIN+i))>>GDIV)-U_G);
  double TMP_BB = 0.0; //(((*(BIN+i))>>BDIV)-U_B)*(((*(BIN+i))>>BDIV)-U_B);
  double TMP_RG = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(GIN+i))>>GDIV)-U_G);
  double TMP_RB = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(BIN+i))>>BDIV)-U_B);
  double TMP_GB = 0.0; //(((*(GIN+i))>>GDIV)-U_G)*(((*(BIN+i))>>BDIV)-U_B);
  double U_R = std::get<0>(U_RGB);
  double U_G = std::get<1>(U_RGB);
  double U_B = std::get<2>(U_RGB);
  for (int i = 0; i < SIZE; i++) {
    if (PIXTYPES[i] != PIXTYPE) {
      NUM_IGNORES++;
      continue;
    }
    TMP_RR += (YIN[i] - U_R) * (YIN[i] - U_R);
    TMP_GG += (UIN[i] - U_G) * (UIN[i] - U_G);
    TMP_BB += (VIN[i] - U_B) * (VIN[i] - U_B);
    TMP_RG += (YIN[i] - U_R) * (UIN[i] - U_G);
    TMP_RB += (YIN[i] - U_R) * (VIN[i] - U_B);
    TMP_GB += (UIN[i] - U_G) * (VIN[i] - U_B);
  }
  SIZE -= NUM_IGNORES;
  TMP_RR /= SIZE;
  TMP_GG /= SIZE;
  TMP_BB /= SIZE;
  TMP_RG /= SIZE;
  TMP_RB /= SIZE;
  TMP_GB /= SIZE;

  A[0][0] = TMP_RR;
  A[0][1] = TMP_RG;
  A[1][0] = TMP_RG;
  A[1][1] = TMP_GG;
  A[2][2] = TMP_BB;
  A[0][2] = TMP_RB;
  A[2][0] = TMP_RB;
  A[1][2] = TMP_GB;
  A[2][1] = TMP_GB;
}

int Jacobi(int ct, double eps, double A[3][3], double A1[3][3], double A2[3][3],
           double X1[3][3], double X2[3][3]) {
  int n = 3;
  double max, s, t, v, sn, cs1;
  int k = 0, ind = 1, p = 0, q = 0;
  // 初期設定
  for (int i1 = 0; i1 < n; i1++) {
    for (int i2 = 0; i2 < n; i2++) {
      A1[i1][i2] = A[i1][i2];
      X1[i1][i2] = 0.0;
    }
    X1[i1][i1] = 1.0;
  }
  // 計算
  while (ind > 0 && k < ct) {
    // 最大要素の探索
    max = 0.0;
    for (int i1 = 0; i1 < n; i1++) {
      for (int i2 = 0; i2 < n; i2++) {
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
      for (int i1 = 0; i1 < n; i1++) {
        if (i1 == p) {
          for (int i2 = 0; i2 < n; i2++) {
            if (i2 == p)
              A2[p][p] = A1[p][p] * cs1 * cs1 + A1[q][q] * sn * sn -
                         2.0 * A1[p][q] * sn * cs1;
            else if (i2 == q)
              A2[p][q] = 0.0;
            else
              A2[p][i2] = A1[p][i2] * cs1 - A1[q][i2] * sn;
          }
        } else if (i1 == q) {
          for (int i2 = 0; i2 < n; i2++) {
            if (i2 == q)
              A2[q][q] = A1[p][p] * sn * sn + A1[q][q] * cs1 * cs1 +
                         2.0 * A1[p][q] * sn * cs1;
            else if (i2 == p)
              A2[q][p] = 0.0;
            else
              A2[q][i2] = A1[q][i2] * cs1 + A1[p][i2] * sn;
          }
        } else {
          for (int i2 = 0; i2 < n; i2++) {
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
      for (int i1 = 0; i1 < n; i1++) {
        for (int i2 = 0; i2 < n; i2++) {
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
      for (int i1 = 0; i1 < n; i1++) {
        for (int i2 = 0; i2 < n; i2++) {
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
  int THRESH = 0;
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

    HANTEI = (float)(TOTAL1) * (float)(TOTAL2)*Square(AVE1 - AVE2);

    if (MAX2 < HANTEI) {
      MAX2 = HANTEI;
      THRESH = i;
    }
  }
  THRESH += MIN;

  return (THRESH);
} // main kansuu end

struct Threshold {
  double value;  // 大小比較する値
  int index = 0; // threshold

  Threshold() : value(0.0) {}
  Threshold(double v) : value(v) {}
};
bool operator<(const Threshold &lhs, const Threshold &rhs) {
  return lhs.value < rhs.value;
}

struct XYZ_t {
  double X, Y, Z;

  XYZ_t &operator+=(const XYZ_t &rhs) {
    this->X += rhs.X;
    this->Y += rhs.Y;
    this->Z += rhs.Z;
    return *this;
  }
  // ベクトル版 {X1+X2, Y1+Y2, Z1+Z2}
  friend XYZ_t operator+(XYZ_t lhs, const XYZ_t &rhs) {
    lhs += rhs;
    return lhs;
  }
};
// ソート用 (X しか見ないことに注意)
bool operator<(const XYZ_t &lhs, const XYZ_t &rhs) { return lhs.X < rhs.X; }
// ベクトル版 {X1-X2, Y1-Y2, Z1-Z2}
XYZ_t operator-(const XYZ_t &lhs, const XYZ_t &rhs) {
  return XYZ_t{lhs.X - rhs.X, lhs.Y - rhs.Y, lhs.Z - rhs.Z};
}
XYZ_t operator*(const XYZ_t &lhs, const XYZ_t &rhs) {
  return XYZ_t{lhs.X * rhs.X, lhs.Y * rhs.Y, lhs.Z * rhs.Z};
}
// スカラ版 {X*a, Y*a, Z*a}
XYZ_t operator/(const XYZ_t &lhs, double rhs) {
  return XYZ_t{lhs.X / rhs, lhs.Y / rhs, lhs.Z / rhs};
}

struct Statistics {
  XYZ_t SUM = {0.0, 0.0, 0.0};   // 和
  XYZ_t SQSUM = {0.0, 0.0, 0.0}; // 二乗和
  int TOTAL = 0;

  // 加算
  void Add(const XYZ_t &v) {
    this->SUM += v;
    this->SQSUM += v * v;
    this->TOTAL++;
  }
  // 平均
  XYZ_t Mean() const { return this->SUM / this->TOTAL; }
  // 分散
  // (X[j] - mean(X[:]))**2 = mean(X[:]**2) - mean(X[:])**2;
  XYZ_t Variance() const {
    XYZ_t MEAN = this->Mean();
    XYZ_t MEANSQ = MEAN * MEAN;
    XYZ_t SQMEAN = this->SQSUM / this->TOTAL;
    return SQMEAN - MEANSQ;
  }
};

int ohtsu2(int NUM, const double *X, const double *Y, const double *Z,
           int omh) {
  assert(omh == 3 || omh == 4);

  std::pair<const double *, const double *> MINMAX =
      std::minmax_element(&X[0], &X[NUM]);
  double MAX = *MINMAX.second;
  double MIN = *MINMAX.first;
  int CMAX = (int)floor(MAX); // マイナス時にも底を取る
  int CMIN;

  if (MIN > -1.0) {
    // -1.0 より大きい時は 0 に近い値に丸める
    CMIN = (int)MIN;
  } else {
    // -1.0 以下の場合に底を取る
    CMIN = floor(MIN);
  }

  int NODEHANI = CMAX - CMIN + 1;

  // X が一定以下のみ集計するため、事前にソートしておくことで効率よくする
  std::vector<XYZ_t> XYZ(NUM);
  for (int i = 0; i < NUM; i++) {
    XYZ[i] = XYZ_t{X[i], Y[i], Z[i]};
  }
  std::sort(XYZ.begin(), XYZ.end());

  std::vector<Statistics> STAT1(NODEHANI);
  {
    Statistics STAT;
    int j = 0;
    for (int i = 1; i < NODEHANI; i++) {
      for (; j < NUM; j++) {
        // X < i + CMIN なもののみ集計する
        const XYZ_t &XYZ_j = XYZ[j];
        if (!(XYZ_j.X < i + CMIN)) {
          break; // X < i + CMIN を満たしていないので終了
        }
        STAT.Add(XYZ_j);
      }
      STAT1[i] = STAT;
    }
  }

  std::vector<Statistics> STAT2(NODEHANI);
  {
    Statistics STAT;
    int j = NUM - 1;
    for (int i = NODEHANI - 1; 0 < i; i--) {
      for (; 0 <= j; j--) {
        // X >= i + CMIN なもののみ集計する
        const XYZ_t &XYZ_j = XYZ[j];
        if (!(XYZ_j.X >= i + CMIN)) {
          break; // X >= i + CMIN を満たしていないので終了
        }
        STAT.Add(XYZ_j);
      }
      STAT2[i] = STAT;
    }
  }

  Threshold THRESH_MAX(0.0);     // 最大値
  Threshold THRESH_MIN(DBL_MAX); // 最小値

  for (int i = 1; i < NODEHANI; i++) {
    if (omh == 3) {
      Threshold HANTEI;
      XYZ_t AVE1 = STAT1[i].Mean();
      XYZ_t AVE2 = STAT2[i].Mean();
      HANTEI.value = (double)(STAT1[i].TOTAL) * (double)(STAT2[i].TOTAL) *
                     Square(AVE1.X - AVE2.X);
      HANTEI.index = i;
      if (THRESH_MAX < HANTEI) {
        THRESH_MAX = HANTEI;
      }
    }
    if (omh == 4) {
      Threshold HANTEI;
      XYZ_t BUN1 = STAT1[i].Variance();
      XYZ_t BUN2 = STAT2[i].Variance();
      HANTEI.value = (BUN1.X + BUN1.Y + BUN1.Z) * (double)STAT1[i].TOTAL +
                     (BUN2.X + BUN2.Y + BUN2.Z) * (double)STAT2[i].TOTAL;
      HANTEI.index = i;
      if (HANTEI < THRESH_MIN) {
        THRESH_MIN = HANTEI;
      }
    }
  }

  THRESH_MAX.index += CMIN;
  THRESH_MIN.index += CMIN;
  // debug start
  fprintf(stderr, "CMIN=%d,THRESH=%d,THRESH2=%d,CMAX=%d\n", CMIN,
          THRESH_MAX.index, THRESH_MIN.index, CMAX);
  //        while(1);
  // debug end

  if (omh == 3) {
    return (THRESH_MAX.index);
  } else if (omh == 4) {
    return (THRESH_MIN.index);
  }
  return (0);
} // main kansuu end

int media(int NUM, const int *X) {
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
  int MED = 0;
  int GOUKEI = 0;
  for (int i = 0; i < NODEHANI; i++) {
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

    HANTEI = (float)(TOTAL1) * (float)(TOTAL2)*Square(AVE1 - AVE2);

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

static void kaiten_one(double V[3], double X, double Y, double Z, double *R,
                       double *G, double *B) {
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

void kaiten(int IMAGE_SIZE, double V[3], const double *X, const double *Y,
            const double *Z, double *R, double *G, double *B) {
  for (int i = 0; i < IMAGE_SIZE; i++) {
    kaiten_one(V, X[i], Y[i], Z[i], &R[i], &G[i], &B[i]);
  }
}

void kaitenEquals(int IMAGE_SIZE, double V[3], const double *X, const double *Y,
                  const double *Z, const int16_t *INDEXED_COLOR,
                  int16_t PALET_NO, double *R, double *G, double *B) {
  for (int i = 0; i < IMAGE_SIZE; i++) {
    if (INDEXED_COLOR[i] == PALET_NO) {
      kaiten_one(V, X[i], Y[i], Z[i], &R[i], &G[i], &B[i]);
    }
  }
}

// TODO: 愚直な引き算に直す
void SobelFilterHorizontal(int hsize, int vsize, const int *IN, int *OUT) {
  int SOBEL2[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
  for (int i = 0; i < vsize; i++) {   // suityoku
    for (int j = 0; j < hsize; j++) { // suihei
      for (int k = 0; k < 3; k++) {   // suityoku
        for (int l = 0; l < 3; l++) { // suihei
          int m = i + k - 1;          //-1~1suiityoku
          int n = j + l - 1;          //-1~1suihei
          if (n < 0) {
            n = -n;
          }
          if (n > (hsize - 1)) {
            n = 2 * (hsize - 1) - n;
          }
          if (m < 0) {
            m = -m;
          }
          if (m > (vsize - 1)) {
            m = 2 * (vsize - 1) - m;
          }
          OUT[i * hsize + j] += SOBEL2[3 * k + l] * IN[m * hsize + n];
        }
      }
    }
  }
}

// TODO: 愚直な引き算に直す
void SobelFilterVertical(int hsize, int vsize, const int *IN, int *OUT) {
  int SOBEL[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
  for (int i = 0; i < vsize; i++) {   // suityoku
    for (int j = 0; j < hsize; j++) { // suihei
      for (int k = 0; k < 3; k++) {   // suityoku
        for (int l = 0; l < 3; l++) { // suihei
          int m = i + k - 1;          //-1~1suiityoku
          int n = j + l - 1;          //-1~1suihei
          if (n < 0) {
            n = -n;
          }
          if (n > (hsize - 1)) {
            n = 2 * (hsize - 1) - n;
          }
          if (m < 0) {
            m = -m;
          }
          if (m > (vsize - 1)) {
            m = 2 * (vsize - 1) - m;
          }
          OUT[i * hsize + j] += SOBEL[3 * k + l] * IN[m * hsize + n];
        }
      }
    }
  }
}

// エッジ検出
//
// 垂直に微分した後、絶対値が et 以上になるものをエッジとする。
void PrewitAbsoluteFilterVertical(int hsize, int vsize, int et, const int *IN,
                                  int *OUT, int *EDGE) {
  for (int i = 0; i < vsize; i++) {
    for (int j = 0; j < hsize; j++) {
      int k = i - 1;
      int l = i + 1;
      if (k < 0) {
        k = -k;
      }
      if (k > (vsize - 1)) {
        k = 2 * (vsize - 1) - k;
      }
      if (l < 0) {
        l = -l;
      }
      if (l > (vsize - 1)) {
        l = 2 * (vsize - 1) - l;
      }
      int DIFF = abs(IN[k * hsize + j] - IN[l * hsize + j]);
      bool IS_EDGE = (DIFF >= et);
      OUT[i * hsize + j] = DIFF;
      EDGE[i * hsize + j] = IS_EDGE;
    }
  }
}

// エッジ検出
//
// 水平に微分した後、絶対値が et 以上になるものをエッジとする。
void PrewitAbsoluteFilterHorizontal(int hsize, int vsize, int et, const int *IN,
                                    int *OUT, int *EDGE) {
  for (int i = 0; i < vsize; i++) {
    for (int j = 0; j < hsize; j++) {
      int k = j - 1;
      int l = j + 1;
      if (k < 0) {
        k = -k;
      }
      if (l > (hsize - 1)) {
        l = 2 * (hsize - 1) - l;
      }
      int DIFF = abs(IN[i * hsize + l] - IN[i * hsize + k]);
      bool IS_EDGE = (DIFF >= et);
      OUT[i * hsize + j] = DIFF;
      EDGE[i * hsize + j] = IS_EDGE;
    }
  }
}

// エッジ検出 (vvv == 0)
int DetectEdge0(int SIZE, int et, const int *HHEDGER, const int *VVEDGER,
                const int *HHEDGEG, const int *VVEDGEG, const int *HHEDGEB,
                const int *VVEDGEB, int16_t *OUT) {
  int EDGE_GASOSUU = 0;

  for (int i = 0; i < SIZE; i++) {
    if ((sqrt(Square((float)HHEDGER[i]) + Square((float)VVEDGER[i])) >=
         et * sqrt(2.0f)) ||
        (sqrt(Square((float)HHEDGEG[i]) + Square((float)VVEDGEG[i])) >=
         et * sqrt(2.0f)) ||
        (sqrt(Square((float)HHEDGEB[i]) + Square((float)VVEDGEB[i])) >=
         et * sqrt(2.0f))) {
      OUT[i] = kEdgePixel;
      EDGE_GASOSUU++;
    } else {
      OUT[i] = kRegularPixel;
    }
  }

  return EDGE_GASOSUU;
}

// エッジ検出 (vvv == 1)
int DetectEdge1(int SIZE, int et, const int *HHEDGER, const int *VVEDGER,
                const int *HHEDGEG, const int *VVEDGEG, const int *HHEDGEB,
                const int *VVEDGEB, int16_t *OUT) {
  int EDGE_GASOSUU = 0;

  for (int i = 0; i < SIZE; i++) {
    if (HHEDGER[i] >= et || VVEDGER[i] >= et || HHEDGEG[i] >= et ||
        VVEDGEG[i] >= et || HHEDGEB[i] >= et || VVEDGEB[i] >= et) {
      OUT[i] = kEdgePixel;
      EDGE_GASOSUU++;
    } else {
      OUT[i] = kRegularPixel;
    }
  }

  return EDGE_GASOSUU;
}

// エッジ検出 (vvv == 6)
int DetectEdge6(int SIZE, int et, const int *HHEDGER, const int *VVEDGER,
                const int *HHEDGEG, const int *VVEDGEG, const int *HHEDGEB,
                const int *VVEDGEB, int16_t *OUT) {
  int EDGE_GASOSUU = 0;

  for (int i = 0; i < SIZE; i++) {
    if (sqrt((Square((float)HHEDGER[i]) + Square((float)VVEDGER[i])) +
             (Square((float)HHEDGEG[i]) + Square((float)VVEDGEG[i])) +
             (Square((float)HHEDGEB[i]) + Square((float)VVEDGEB[i]))) >=
        et * sqrt(6.0f)) {
      OUT[i] = kEdgePixel;
      EDGE_GASOSUU++;
    } else {
      OUT[i] = kRegularPixel;
    }
  }

  return EDGE_GASOSUU;
}

// エッジ検出 (vvv == 2)
int DetectEdge2(int SIZE, const int *EDGER, const int *VEDGER, const int *EDGEG,
                const int *VEDGEG, const int *EDGEB, const int *VEDGEB,
                int16_t *OUT) {
  int EDGE_GASOSUU = 0;

  for (int i = 0; i < SIZE; i++) {
    if (EDGER[i] == 1 || EDGEG[i] == 1 || EDGEB[i] == 1 || VEDGER[i] == 1 ||
        VEDGEG[i] == 1 || VEDGEB[i] == 1) {
      OUT[i] = kEdgePixel;
      EDGE_GASOSUU++;
    } else {
      OUT[i] = kRegularPixel;
    }
  }

  return EDGE_GASOSUU;
}

// エッジ検出 (vvv == 3)
int DetectEdge3(int SIZE, int et, const int *HHEDGER, const int *VVEDGER,
                const int *HHEDGEG, const int *VVEDGEG, const int *HHEDGEB,
                const int *VVEDGEB, int16_t *OUT) {
  int EDGE_GASOSUU = 0;

  for (int i = 0; i < SIZE; i++) {
    if ((sqrt(Square((float)HHEDGER[i]) + Square((float)VVEDGER[i])) >=
         et * 4.0 * sqrt(2.0)) ||
        (sqrt((Square((float)HHEDGEG[i]) + Square((float)VVEDGEG[i])) >=
              et * 4.0 * sqrt(2.0)) ||
         (sqrt(Square((float)HHEDGEB[i]) + Square((float)VVEDGEB[i])) >=
          et * 4.0 * sqrt(2.0)))) {
      OUT[i] = kEdgePixel;
      EDGE_GASOSUU++;
    } else {
      OUT[i] = kRegularPixel;
    }
  }

  return EDGE_GASOSUU;
}

// エッジ検出 (vvv == 5)
int DetectEdge5(int SIZE, int et, const int *HHEDGER, const int *VVEDGER,
                const int *HHEDGEG, const int *VVEDGEG, const int *HHEDGEB,
                const int *VVEDGEB, int16_t *OUT) {
  int EDGE_GASOSUU = 0;

  for (int i = 0; i < SIZE; i++) {
    if (sqrt((Square((float)HHEDGER[i]) + Square((float)VVEDGER[i])) +
             (Square((float)HHEDGEG[i]) + Square((float)VVEDGEG[i])) +
             (Square((float)HHEDGEB[i]) + Square((float)VVEDGEB[i]))) >=
        et * 4.0 * sqrt(6.0)) {
      OUT[i] = kEdgePixel;
      EDGE_GASOSUU++;
    } else {
      OUT[i] = kRegularPixel;
    }
  }

  return EDGE_GASOSUU;
}

// エッジ検出 (vvv == 4)
int DetectEdge4(int SIZE, const int *HHEDGER, const int *VVEDGER,
                const int *HHEDGEG, const int *VVEDGEG, const int *HHEDGEB,
                const int *VVEDGEB, int16_t *OUT, double *EDGERASISAY) {
  for (int i = 0; i < SIZE; i++) {
    OUT[i] = kRegularPixel;
    double REDGE = /*sqrt*/ (Square((float)HHEDGER[i]) +
                             Square((float)VVEDGER[i])) /*/ 4.0 / sqrt(2.0)*/;
    double GEDGE = /*sqrt*/ (Square((float)HHEDGEG[i]) +
                             Square((float)VVEDGEG[i])) /*/ 4.0 / sqrt(2.0)*/;
    double BEDGE = /*sqrt*/ (Square((float)HHEDGEB[i]) +
                             Square((float)VVEDGEB[i])) /*/ 4.0 / sqrt(2.0)*/;
    //            EDGERASISAY[i] =
    //            0.29891*REDGE+0.58661*GEDGE+0.11448*BEDGE;
    EDGERASISAY[i] = (sqrt(REDGE + GEDGE + BEDGE)) / 4.0 / sqrt(6.0);
    //            if(EDGERASISAY[i] > EDGERASMAX){
    //                EDGERASMAX = EDGERASISAY[i];
    //            }
  }

  return 0;
}

// PIXTYPES[i] == PIXTYPE である画素の和
std::tuple<double, double, double> SumEquals(int SIZE, const int16_t *PIXTYPES,
                                             int16_t PIXTYPE, const double *v1,
                                             const double *v2,
                                             const double *v3) {
  double SUM1 = 0;
  double SUM2 = 0;
  double SUM3 = 0;

  for (int i = 0; i < SIZE; i++) {
    if (PIXTYPES[i] != PIXTYPE) {
      continue;
    }
    SUM1 += v1[i];
    SUM2 += v2[i];
    SUM3 += v3[i];
  }

  return std::make_tuple(SUM1, SUM2, SUM3);
}

// PIXTYPES[i] != PIXTYPE である画素の和
std::tuple<double, double, double>
SumNotEquals(int SIZE, const int16_t *PIXTYPES, int16_t PIXTYPE,
             const double *v1, const double *v2, const double *v3) {
  double SUM1 = 0;
  double SUM2 = 0;
  double SUM3 = 0;

  for (int i = 0; i < SIZE; i++) {
    if (PIXTYPES[i] == PIXTYPE) {
      continue;
    }
    SUM1 += v1[i];
    SUM2 += v2[i];
    SUM3 += v3[i];
  }

  return std::make_tuple(SUM1, SUM2, SUM3);
}

// (v - a)**2 の和が最小になるインデックスを探す
static int FindMinimumDistanceIndex(int SIZE, double v1, double v2, double v3,
                                    const uint8_t *a1, const uint8_t *a2,
                                    const uint8_t *a3) {
  int MIN_INDEX = 0;
  double MIN_VALUE = DBL_MAX;
  for (int i = 0; i < SIZE; i++) {
    double d = Square(v1 - a1[i]) + Square(v2 - a2[i]) + Square(v3 - a3[i]);
    if (d < MIN_VALUE) {
      MIN_VALUE = d;
      MIN_INDEX = i;
    }
  }
  return MIN_INDEX;
}
// (v - a)**2 の和が最小になる a の要素のポインタを返す
static const JYUSHIN_INDEX<double> *
FindMinimumDistance(int SIZE, const JYUSHIN_INDEX<double> &v,
                    const JYUSHIN_INDEX<double> *a) {
  const JYUSHIN_INDEX<double> *MIN_ITEM = nullptr;
  double MIN_VALUE = DBL_MAX;
  for (int i = 0; i < SIZE; i++) {
    double d =
        Square(v.Y - a[i].Y) + Square(v.U - a[i].U) + Square(v.V - a[i].V);
    if (d < MIN_VALUE) {
      MIN_VALUE = d;
      MIN_ITEM = &a[i];
    }
  }
  return MIN_ITEM;
}

// パレット番号に変換し OUT に出力する
void Dithering(int IMAGE_SIZE, const double *YIN, const double *UIN,
               const double *VIN, int PALETTE_SIZE, const double *Y_JYUSHIN,
               const double *U_JYUSHIN, const double *V_JYUSHIN, uint8_t *OUT) {
  for (int i = 0; i < IMAGE_SIZE; i++) {
    double Z[PALETTE_SIZE];
    for (int j = 0; j < PALETTE_SIZE; j++) {
      Z[j] = Square(YIN[i] - Y_JYUSHIN[j]) + Square(UIN[i] - U_JYUSHIN[j]) +
             Square(VIN[i] - V_JYUSHIN[j]);
    }
    ptrdiff_t MINZ_INDEX = std::min_element(&Z[0], &Z[PALETTE_SIZE]) - &Z[0];
    OUT[i] = static_cast<uint8_t>(MINZ_INDEX);
  }
}

// パレット番号に変換し OUT に出力する (edon == 1)
void Dithering1(int hsize, int vsize, int dither, double per,
                const uint8_t *RIN, const uint8_t *GIN, const uint8_t *BIN,
                int PALETTE_SIZE, const uint8_t *REDUCE_R,
                const uint8_t *REDUCE_G, const uint8_t *REDUCE_B,
                uint8_t *OUT) {
  int IMAGE_SIZE = hsize * vsize;
  std::vector<uint8_t> IIRIN(IMAGE_SIZE);
  std::vector<uint8_t> IIGIN(IMAGE_SIZE);
  std::vector<uint8_t> IIBIN(IMAGE_SIZE);
  for (int i = 0; i < IMAGE_SIZE; i++) {
    IIRIN[i] = RIN[i];
    IIGIN[i] = GIN[i];
    IIBIN[i] = BIN[i];
  }
  int mx = 0;
  if (dither == 0 || dither == 1) {
    mx = hsize + 2;
  } else if (dither == 2 || dither == 3) {
    mx = hsize + 4;
  }
  int ERROR_SIZE = 0;
  if (dither == 0 || dither == 1) {
    ERROR_SIZE = mx * 2;
  } else if (dither == 2 || dither == 3) {
    ERROR_SIZE = mx * 3;
  }
  std::vector<double> errorR(ERROR_SIZE);
  std::vector<double> errorG(ERROR_SIZE);
  std::vector<double> errorB(ERROR_SIZE);
  double r = 0.0;
  double g = 0.0;
  double b = 0.0;
  double re = 0.0;
  double ge = 0.0;
  double be = 0.0;
  int adr = 0;
  for (int y = 0; y < vsize; y++) {
    for (int x = 0; x < hsize; x++) {
      if (dither == 0 || dither == 1) {
        adr = x + 1;
      } else if (dither == 2 || dither == 3) {
        adr = x + 2;
      }
      // r= (IIRIN[x*vsize+y])+((errorR[adr])/16);
      // g= (IIGIN[x*vsize+y])+((errorG[adr])/16);
      // b= (IIBIN[x*vsize+y])+((errorB[adr])/16);
      if (dither == 0) {
        r = (double)(RIN[x + y * hsize]) +
            ((errorR[adr]) / (16.0 / (per + 0.0000001)) /*32.0*/);
        g = (double)(GIN[x + y * hsize]) +
            ((errorG[adr]) / (16.0 / (per + 0.0000001)) /*32.0*/);
        b = (double)(BIN[x + y * hsize]) +
            ((errorB[adr]) / (16.0 / (per + 0.0000001)) /*32.0*/);
      } else if (dither == 1) {
        r = (double)(RIN[x + y * hsize]) +
            ((errorR[adr]) / (4.0 / (per + 0.0000001)));
        g = (double)(GIN[x + y * hsize]) +
            ((errorG[adr]) / (4.0 / (per + 0.0000001)));
        b = (double)(BIN[x + y * hsize]) +
            ((errorB[adr]) / (4.0 / (per + 0.0000001)));
      } else if (dither == 2) {
        r = (double)(RIN[x + y * hsize]) +
            ((errorR[adr]) / (42.0 / (per + 0.0000001)));
        g = (double)(GIN[x + y * hsize]) +
            ((errorG[adr]) / (42.0 / (per + 0.0000001)));
        b = (double)(BIN[x + y * hsize]) +
            ((errorB[adr]) / (42.0 / (per + 0.0000001)));
      } else if (dither == 3) {
        r = (double)(RIN[x + y * hsize]) +
            ((errorR[adr]) / (48.0 / (per + 0.0000001)));
        g = (double)(GIN[x + y * hsize]) +
            ((errorG[adr]) / (48.0 / (per + 0.0000001)));
        b = (double)(BIN[x + y * hsize]) +
            ((errorB[adr]) / (48.0 / (per + 0.0000001)));
      }
      /* if(r > 255.0){
      r=255.0;
      }
      if(r < 0.0){
      r=0.0;
      }
      if(g > 255.0){
      g=255.0;
      }
      if(g < 0.0){
      g=0.0;
      }
      if(b > 255.0){
      b=255.0;
      }
      if(b < 0.0){
      b=0.0;
      }*/
      // r= (IIRIN[x*vsize+y]);//+((errorR[adr])/16);
      // g= (IIGIN[x*vsize+y]);//+((errorG[adr])/16);
      // b= (IIBIN[x*vsize+y]);//+((errorB[adr])/16);
      int bst = FindMinimumDistanceIndex(PALETTE_SIZE, r, g, b, &REDUCE_R[0],
                                         &REDUCE_G[0], &REDUCE_B[0]);
      // debug start
      // printf("gosa=%d\n",(int)est);
      // debug end
      re = r - (double)(REDUCE_R[bst]);
      ge = g - (double)(REDUCE_G[bst]);
      be = b - (double)(REDUCE_B[bst]);
      // (errorR[adr+1]) += re*7;
      // (errorG[adr+1]) += ge*7;
      // (errorB[adr+1]) += be*7;

      // (errorR[adr+mx-1]) += re*3;
      // (errorG[adr+mx-1]) += ge*3;
      // (errorB[adr+mx-1]) += be*3;

      // (errorR[adr+mx]) += re*5;
      // (errorG[adr+mx]) += ge*5;
      // (errorB[adr+mx]) += be*5;

      // (errorR[adr+mx+1]) += re;
      // (errorG[adr+mx+1]) += ge;
      // (errorB[adr+mx+1]) += be;
      if (dither == 0) {
        (errorR[adr + 1]) += re * 7.0;
        (errorG[adr + 1]) += ge * 7.0;
        (errorB[adr + 1]) += be * 7.0;

        (errorR[adr + mx - 1]) += re * 3.0;
        (errorG[adr + mx - 1]) += ge * 3.0;
        (errorB[adr + mx - 1]) += be * 3.0;

        (errorR[adr + mx]) += re * 5.0;
        (errorG[adr + mx]) += ge * 5.0;
        (errorB[adr + mx]) += be * 5.0;

        (errorR[adr + mx + 1]) += re;
        (errorG[adr + mx + 1]) += ge;
        (errorB[adr + mx + 1]) += be;
        // (errorR[adr+mx]) += re*1;
        // (errorG[adr+mx]) += ge*1;
        // (errorB[adr+mx]) += be*1;
      } else if (dither == 1) {
        (errorR[adr + 1]) += re * 2.0;
        (errorG[adr + 1]) += ge * 2.0;
        (errorB[adr + 1]) += be * 2.0;

        (errorR[adr + mx - 1]) += re;
        (errorG[adr + mx - 1]) += ge;
        (errorB[adr + mx - 1]) += be;

        (errorR[adr + mx]) += re;
        (errorG[adr + mx]) += ge;
        (errorB[adr + mx]) += be;

        // (errorR[adr+mx+1]) += re;
        // (errorG[adr+mx+1]) += ge;
        // (errorB[adr+mx+1]) += be;

      } else if (dither == 2) {
        (errorR[adr + 1]) += re * 8.0;
        (errorG[adr + 1]) += ge * 8.0;
        (errorB[adr + 1]) += be * 8.0;
        (errorR[adr + 2]) += re * 4.0;
        (errorG[adr + 2]) += ge * 4.0;
        (errorB[adr + 2]) += be * 4.0;
        (errorR[adr + mx - 2]) += re * 2.0;
        (errorG[adr + mx - 2]) += ge * 2.0;
        (errorB[adr + mx - 2]) += be * 2.0;
        (errorR[adr + mx - 1]) += re * 4.0;
        (errorG[adr + mx - 1]) += ge * 4.0;
        (errorB[adr + mx - 1]) += be * 4.0;
        (errorR[adr + mx]) += re * 8.0;
        (errorG[adr + mx]) += ge * 8.0;
        (errorB[adr + mx]) += be * 8.0;
        (errorR[adr + mx + 1]) += re * 4.0;
        (errorG[adr + mx + 1]) += ge * 4.0;
        (errorB[adr + mx + 1]) += be * 4.0;
        (errorR[adr + mx + 2]) += re * 2.0;
        (errorG[adr + mx + 2]) += ge * 2.0;
        (errorB[adr + mx + 2]) += be * 2.0;

        (errorR[adr + mx * 2 - 2]) += re;
        (errorG[adr + mx * 2 - 2]) += ge;
        (errorB[adr + mx * 2 - 2]) += be;
        (errorR[adr + mx * 2 - 1]) += re * 2.0;
        (errorG[adr + mx * 2 - 1]) += ge * 2.0;
        (errorB[adr + mx * 2 - 1]) += be * 2.0;
        (errorR[adr + mx * 2]) += re * 4.0;
        (errorG[adr + mx * 2]) += ge * 4.0;
        (errorB[adr + mx * 2]) += be * 4.0;
        (errorR[adr + mx * 2 + 1]) += re * 2.0;
        (errorG[adr + mx * 2 + 1]) += ge * 2.0;
        (errorB[adr + mx * 2 + 1]) += be * 2.0;
        (errorR[adr + mx * 2 + 2]) += re;
        (errorG[adr + mx * 2 + 2]) += ge;
        (errorB[adr + mx * 2 + 2]) += be;
      } else if (dither == 3) {
        (errorR[adr + 1]) += re * 7.0;
        (errorG[adr + 1]) += ge * 7.0;
        (errorB[adr + 1]) += be * 7.0;
        (errorR[adr + 2]) += re * 5.0;
        (errorG[adr + 2]) += ge * 5.0;
        (errorB[adr + 2]) += be * 5.0;
        (errorR[adr + mx - 2]) += re * 3.0;
        (errorG[adr + mx - 2]) += ge * 3.0;
        (errorB[adr + mx - 2]) += be * 3.0;
        (errorR[adr + mx - 1]) += re * 4.0;
        (errorG[adr + mx - 1]) += ge * 4.0;
        (errorB[adr + mx - 1]) += be * 4.0;
        (errorR[adr + mx]) += re * 7.0;
        (errorG[adr + mx]) += ge * 7.0;
        (errorB[adr + mx]) += be * 7.0;
        (errorR[adr + mx + 1]) += re * 5.0;
        (errorG[adr + mx + 1]) += ge * 5.0;
        (errorB[adr + mx + 1]) += be * 5.0;
        (errorR[adr + mx + 2]) += re * 3.0;
        (errorG[adr + mx + 2]) += ge * 3.0;
        (errorB[adr + mx + 2]) += be * 3.0;

        (errorR[adr + mx * 2 - 2]) += re;
        (errorG[adr + mx * 2 - 2]) += ge;
        (errorB[adr + mx * 2 - 2]) += be;
        (errorR[adr + mx * 2 - 1]) += re * 3.0;
        (errorG[adr + mx * 2 - 1]) += ge * 3.0;
        (errorB[adr + mx * 2 - 1]) += be * 3.0;
        (errorR[adr + mx * 2]) += re * 5.0;
        (errorG[adr + mx * 2]) += ge * 5.0;
        (errorB[adr + mx * 2]) += be * 5.0;
        (errorR[adr + mx * 2 + 1]) += re * 3.0;
        (errorG[adr + mx * 2 + 1]) += ge * 3.0;
        (errorB[adr + mx * 2 + 1]) += be * 3.0;
        (errorR[adr + mx * 2 + 2]) += re;
        (errorG[adr + mx * 2 + 2]) += ge;
        (errorB[adr + mx * 2 + 2]) += be;
      }

      (IIRIN[x + y * hsize]) = (REDUCE_R[bst]);
      (IIGIN[x + y * hsize]) = (REDUCE_G[bst]);
      (IIBIN[x + y * hsize]) = (REDUCE_B[bst]);
      // (IIRIN[x*vsize+y]) = (int)(REDUCE_R[OUT[x*vsize+y]]);
      // (IIGIN[x*vsize+y]) = (int)(REDUCE_G[OUT[x*vsize+y]]);
      // (IIBIN[x*vsize+y]) = (int)(REDUCE_B[OUT[x*vsize+y]]);
    }
    if (dither == 0 || dither == 1) {
      for (int j = 0; j < mx; j++) {
        errorR[j] = errorR[j + mx];
        errorG[j] = errorG[j + mx];
        errorB[j] = errorB[j + mx];
        errorR[j + mx] = 0.0;
        errorG[j + mx] = 0.0;
        errorB[j + mx] = 0.0;
      }
    } else if (dither == 2 || dither == 3) {
      for (int j = 0; j < mx; j++) {
        errorR[j] = errorR[j + mx];
        errorG[j] = errorG[j + mx];
        errorB[j] = errorB[j + mx];
        errorR[j + mx] = errorR[j + 2 * mx];
        errorG[j + mx] = errorG[j + 2 * mx];
        errorB[j + mx] = errorB[j + 2 * mx];
        errorR[j + 2 * mx] = 0.0;
        errorG[j + 2 * mx] = 0.0;
        errorB[j + 2 * mx] = 0.0;
      }
    }
  }

  // OUT[i]
  // IIRINに画像データがはいっている

  for (int i = 0; i < IMAGE_SIZE; i++) {
    for (int j = 0; j < PALETTE_SIZE; j++) {
      if ((IIRIN[i] == (REDUCE_R[j])) && (IIGIN[i] == (REDUCE_G[j])) &&
          (IIBIN[i] == (REDUCE_B[j]))) {
        OUT[i] = j;
      }
    }
  }
}

// lhs < rhs を V, Y, U の順に比較 (INDEX は比較しないことに注意)
template <typename T>
static inline bool VYULess(const JYUSHIN_INDEX<T> &lhs,
                           const JYUSHIN_INDEX<T> &rhs) {
  if (lhs.V == rhs.V) {
    if (lhs.Y == rhs.Y) {
      return lhs.U < rhs.U;
    }
    return lhs.Y < rhs.Y;
  }
  return lhs.V < rhs.V;
}

// Y, U, V が等価 (INDEX は比較しないことに注意)
template <typename T>
static inline bool YUVEqual(const JYUSHIN_INDEX<T> &lhs,
                            const JYUSHIN_INDEX<T> &rhs) {
  return (lhs.V == rhs.V) && (lhs.Y == rhs.Y) && (lhs.U == rhs.U);
}

// パレット番号に変換し OUT に出力する (edon == 2)
void Dithering2(int hsize, int vsize, int dither, double per, const double *VIN,
                const double *YIN, const double *UIN, int PALETTE_SIZE,
                const JYUSHIN_INDEX<double> *JYUSHIN, uint8_t *OUT) {
  int IMAGE_SIZE = hsize * vsize;

  std::vector<int> IIRIN(IMAGE_SIZE);
  std::vector<int> IIGIN(IMAGE_SIZE);
  std::vector<int> IIBIN(IMAGE_SIZE);
  for (int i = 0; i < IMAGE_SIZE; i++) {
    IIRIN[i] = static_cast<int>(VIN[i]);
    IIGIN[i] = static_cast<int>(YIN[i]);
    IIBIN[i] = static_cast<int>(UIN[i]);
  }
  // V, Y, U の順にソートされた JYUSHIN
  std::vector<JYUSHIN_INDEX<int32_t>> IJYUSHIN_VYU_SORTED(PALETTE_SIZE);
  for (int i = 0; i < PALETTE_SIZE; i++) {
    IJYUSHIN_VYU_SORTED[i] = JYUSHIN_INDEX<int32_t>{
        static_cast<int32_t>(JYUSHIN[i].Y), static_cast<int32_t>(JYUSHIN[i].U),
        static_cast<int32_t>(JYUSHIN[i].V), JYUSHIN[i].INDEX};
  }
  std::sort(IJYUSHIN_VYU_SORTED.begin(), IJYUSHIN_VYU_SORTED.end(),
            VYULess<int32_t>);

  int mx = 0;
  if (dither == 0 || dither == 1) {
    mx = hsize + 2;
  } else if (dither == 2 || dither == 3) {
    mx = hsize + 4;
  }
  int ERROR_SIZE = 0;
  if (dither == 0 || dither == 1) {
    ERROR_SIZE = mx * 2;
  } else if (dither == 2 || dither == 3) {
    ERROR_SIZE = mx * 3;
  }
  std::vector<double> errorR(ERROR_SIZE);
  std::vector<double> errorG(ERROR_SIZE);
  std::vector<double> errorB(ERROR_SIZE);
  double r = 0.0;
  double g = 0.0;
  double b = 0.0;
  double re = 0.0;
  double ge = 0.0;
  double be = 0.0;
  int adr = 0;
  for (int y = 0; y < vsize; y++) {
    for (int x = 0; x < hsize; x++) {
      if (dither == 0 || dither == 1) {
        adr = x + 1;
      } else if (dither == 2 || dither == 3) {
        adr = x + 2;
      }

      if (dither == 0) {
        r = (double)(VIN[x + y * hsize]) +
            ((errorR[adr]) / (16.0 / (per + 0.0000001)) /*32.0*/);
        g = (double)(YIN[x + y * hsize]) +
            ((errorG[adr]) / (16.0 / (per + 0.0000001)) /*32.0*/);
        b = (double)(UIN[x + y * hsize]) +
            ((errorB[adr]) / (16.0 / (per + 0.0000001)) /*32.0*/);
      } else if (dither == 1) {
        r = (double)(VIN[x + y * hsize]) +
            ((errorR[adr]) / (4.0 / (per + 0.0000001)));
        g = (double)(YIN[x + y * hsize]) +
            ((errorG[adr]) / (4.0 / (per + 0.0000001)));
        b = (double)(UIN[x + y * hsize]) +
            ((errorB[adr]) / (4.0 / (per + 0.0000001)));
      } else if (dither == 2) {
        r = (double)(VIN[x + y * hsize]) +
            ((errorR[adr]) / (42.0 / (per + 0.0000001)));
        g = (double)(YIN[x + y * hsize]) +
            ((errorG[adr]) / (42.0 / (per + 0.0000001)));
        b = (double)(UIN[x + y * hsize]) +
            ((errorB[adr]) / (42.0 / (per + 0.0000001)));
      } else if (dither == 3) {
        r = (double)(VIN[x + y * hsize]) +
            ((errorR[adr]) / (48.0 / (per + 0.0000001)));
        g = (double)(YIN[x + y * hsize]) +
            ((errorG[adr]) / (48.0 / (per + 0.0000001)));
        b = (double)(UIN[x + y * hsize]) +
            ((errorB[adr]) / (48.0 / (per + 0.0000001)));
      }

      // 関数の入り口で RGB = VYU にしているので YUV = GBR
      JYUSHIN_INDEX<double> YUV{g, b, r, 0};
      const JYUSHIN_INDEX<double> *bst =
          FindMinimumDistance(PALETTE_SIZE, YUV, &JYUSHIN[0]);

      re = r - bst->V;
      ge = g - bst->Y;
      be = b - bst->U;

      if (dither == 0) {
        (errorR[adr + 1]) += re * 7.0;
        (errorG[adr + 1]) += ge * 7.0;
        (errorB[adr + 1]) += be * 7.0;

        (errorR[adr + mx - 1]) += re * 3.0;
        (errorG[adr + mx - 1]) += ge * 3.0;
        (errorB[adr + mx - 1]) += be * 3.0;

        (errorR[adr + mx]) += re * 5.0;
        (errorG[adr + mx]) += ge * 5.0;
        (errorB[adr + mx]) += be * 5.0;

        (errorR[adr + mx + 1]) += re;
        (errorG[adr + mx + 1]) += ge;
        (errorB[adr + mx + 1]) += be;
        // (errorR[adr+mx]) += re*1;
        // (errorG[adr+mx]) += ge*1;
        // (errorB[adr+mx]) += be*1;
      } else if (dither == 1) {
        (errorR[adr + 1]) += re * 2.0;
        (errorG[adr + 1]) += ge * 2.0;
        (errorB[adr + 1]) += be * 2.0;

        (errorR[adr + mx - 1]) += re;
        (errorG[adr + mx - 1]) += ge;
        (errorB[adr + mx - 1]) += be;

        (errorR[adr + mx]) += re;
        (errorG[adr + mx]) += ge;
        (errorB[adr + mx]) += be;

        // (errorR[adr+mx+1]) += re;
        // (errorG[adr+mx+1]) += ge;
        // (errorB[adr+mx+1]) += be;

      } else if (dither == 2) {
        (errorR[adr + 1]) += re * 8.0;
        (errorG[adr + 1]) += ge * 8.0;
        (errorB[adr + 1]) += be * 8.0;
        (errorR[adr + 2]) += re * 4.0;
        (errorG[adr + 2]) += ge * 4.0;
        (errorB[adr + 2]) += be * 4.0;
        (errorR[adr + mx - 2]) += re * 2.0;
        (errorG[adr + mx - 2]) += ge * 2.0;
        (errorB[adr + mx - 2]) += be * 2.0;
        (errorR[adr + mx - 1]) += re * 4.0;
        (errorG[adr + mx - 1]) += ge * 4.0;
        (errorB[adr + mx - 1]) += be * 4.0;
        (errorR[adr + mx]) += re * 8.0;
        (errorG[adr + mx]) += ge * 8.0;
        (errorB[adr + mx]) += be * 8.0;
        (errorR[adr + mx + 1]) += re * 4.0;
        (errorG[adr + mx + 1]) += ge * 4.0;
        (errorB[adr + mx + 1]) += be * 4.0;
        (errorR[adr + mx + 2]) += re * 2.0;
        (errorG[adr + mx + 2]) += ge * 2.0;
        (errorB[adr + mx + 2]) += be * 2.0;

        (errorR[adr + mx * 2 - 2]) += re;
        (errorG[adr + mx * 2 - 2]) += ge;
        (errorB[adr + mx * 2 - 2]) += be;
        (errorR[adr + mx * 2 - 1]) += re * 2.0;
        (errorG[adr + mx * 2 - 1]) += ge * 2.0;
        (errorB[adr + mx * 2 - 1]) += be * 2.0;
        (errorR[adr + mx * 2]) += re * 4.0;
        (errorG[adr + mx * 2]) += ge * 4.0;
        (errorB[adr + mx * 2]) += be * 4.0;
        (errorR[adr + mx * 2 + 1]) += re * 2.0;
        (errorG[adr + mx * 2 + 1]) += ge * 2.0;
        (errorB[adr + mx * 2 + 1]) += be * 2.0;
        (errorR[adr + mx * 2 + 2]) += re;
        (errorG[adr + mx * 2 + 2]) += ge;
        (errorB[adr + mx * 2 + 2]) += be;
      } else if (dither == 3) {
        (errorR[adr + 1]) += re * 7.0;
        (errorG[adr + 1]) += ge * 7.0;
        (errorB[adr + 1]) += be * 7.0;
        (errorR[adr + 2]) += re * 5.0;
        (errorG[adr + 2]) += ge * 5.0;
        (errorB[adr + 2]) += be * 5.0;
        (errorR[adr + mx - 2]) += re * 3.0;
        (errorG[adr + mx - 2]) += ge * 3.0;
        (errorB[adr + mx - 2]) += be * 3.0;
        (errorR[adr + mx - 1]) += re * 4.0;
        (errorG[adr + mx - 1]) += ge * 4.0;
        (errorB[adr + mx - 1]) += be * 4.0;
        (errorR[adr + mx]) += re * 7.0;
        (errorG[adr + mx]) += ge * 7.0;
        (errorB[adr + mx]) += be * 7.0;
        (errorR[adr + mx + 1]) += re * 5.0;
        (errorG[adr + mx + 1]) += ge * 5.0;
        (errorB[adr + mx + 1]) += be * 5.0;
        (errorR[adr + mx + 2]) += re * 3.0;
        (errorG[adr + mx + 2]) += ge * 3.0;
        (errorB[adr + mx + 2]) += be * 3.0;

        (errorR[adr + mx * 2 - 2]) += re;
        (errorG[adr + mx * 2 - 2]) += ge;
        (errorB[adr + mx * 2 - 2]) += be;
        (errorR[adr + mx * 2 - 1]) += re * 3.0;
        (errorG[adr + mx * 2 - 1]) += ge * 3.0;
        (errorB[adr + mx * 2 - 1]) += be * 3.0;
        (errorR[adr + mx * 2]) += re * 5.0;
        (errorG[adr + mx * 2]) += ge * 5.0;
        (errorB[adr + mx * 2]) += be * 5.0;
        (errorR[adr + mx * 2 + 1]) += re * 3.0;
        (errorG[adr + mx * 2 + 1]) += ge * 3.0;
        (errorB[adr + mx * 2 + 1]) += be * 3.0;
        (errorR[adr + mx * 2 + 2]) += re;
        (errorG[adr + mx * 2 + 2]) += ge;
        (errorB[adr + mx * 2 + 2]) += be;
      }

      IIRIN[x + y * hsize] = static_cast<int>(bst->V);
      IIGIN[x + y * hsize] = static_cast<int>(bst->Y);
      IIBIN[x + y * hsize] = static_cast<int>(bst->U);
      // (IIRIN[x*vsize+y]) = (int)(REDUCE_R[OUT[x*vsize+y]]);
      // (IIGIN[x*vsize+y]) = (int)(REDUCE_G[OUT[x*vsize+y]]);
      // (IIBIN[x*vsize+y]) = (int)(REDUCE_B[OUT[x*vsize+y]]);
    }
    if (dither == 0 || dither == 1) {
      for (int j = 0; j < mx; j++) {
        errorR[j] = errorR[j + mx];
        errorG[j] = errorG[j + mx];
        errorB[j] = errorB[j + mx];
        errorR[j + mx] = 0.0;
        errorG[j + mx] = 0.0;
        errorB[j + mx] = 0.0;
      }
    } else if (dither == 2 || dither == 3) {
      for (int j = 0; j < mx; j++) {
        errorR[j] = errorR[j + mx];
        errorG[j] = errorG[j + mx];
        errorB[j] = errorB[j + mx];
        errorR[j + mx] = errorR[j + 2 * mx];
        errorG[j + mx] = errorG[j + 2 * mx];
        errorB[j + mx] = errorB[j + 2 * mx];
        errorR[j + 2 * mx] = 0.0;
        errorG[j + 2 * mx] = 0.0;
        errorB[j + 2 * mx] = 0.0;
      }
    }
  }

  // OUT[i]
  // IIRINに画像データがはいっている
  // 画素と等価な色を IJYUSHIN_VYU_SORTED から探す
  // IJYUSHIN_VYU_SORTED[n].INDEX にパレット番号が入っている
  const JYUSHIN_INDEX<int32_t> *BEGIN = &IJYUSHIN_VYU_SORTED[0];
  const JYUSHIN_INDEX<int32_t> *END = &IJYUSHIN_VYU_SORTED[PALETTE_SIZE];
  for (int i = 0; i < IMAGE_SIZE; i++) {
    // 関数の入り口で RGB = VYU にしているので YUV = GBR
    JYUSHIN_INDEX<int32_t> YUV{IIGIN[i], IIBIN[i], IIRIN[i], 0};
    auto pIJYUSHIN_VYU = std::lower_bound(BEGIN, END, YUV, VYULess<int32_t>);
    if (pIJYUSHIN_VYU == END || !YUVEqual(*pIJYUSHIN_VYU, YUV)) {
      continue; // 見つからなかった
    }
    OUT[i] = pIJYUSHIN_VYU->INDEX;
  }
}