#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tuple>
#include <vector>

#include <chrono>
struct StopWatch {
  StopWatch() { pre_ = std::chrono::high_resolution_clock::now(); }

  //前回のlap関数コールからの経過時間をmilli sec単位で返す
  double lap() {
    auto tmp = std::chrono::high_resolution_clock::now(); // 計測終了時刻を保存
    auto dur = tmp - pre_;
    pre_ = tmp;
    return std::chrono::duration_cast<std::chrono::nanoseconds>(dur).count() /
           1000000.0;
  }
  std::chrono::high_resolution_clock::time_point pre_;
};

#define RD 0
#define GD 0
#define BD 0
#define HASYORI 20
#define MHASYORI 250
#define IROSUU 256
#define YGAMMA 1.0
#define YMULT 6.0 // 1.0*255.0/pow(255.0,YGAMMA)*6.0
#define RMULT 3.0 // 1.0*3.0
#define BMULT 1.0 // 1.0*1.0
#define EDGETH 30 // 25

void JacobiSetupNotEquals(int SIZE, const int *EDGE, int EDGE_NUM,
                          const double *YIN, const double *UIN,
                          const double *VIN,
                          const std::tuple<double, double, double> &U_RGB,
                          double **A);
void JacobiSetupEquals(int SIZE, const int *EDGE, int EDGE_NUM,
                       const double *YIN, const double *UIN, const double *VIN,
                       const std::tuple<double, double, double> &U_RGB,
                       double **A);
int Jacobi(int n, int ct, double eps, double **A, double **A1, double **A2,
           double **X1, double **X2);
void LuvtoRGB(double L, double u, double v, double *R, double *G, double *B);
void RGBtoLuv(double R, double G, double B, double *L, double *u, double *v);
void LABtoRGB(double L, double a, double b, double *R, double *G, double *B);
void RGBtoLAB(double R, double G, double B, double *L, double *a, double *b);
int ohtsu(int NUM, const int *X);
int ohtsu2(int NUM, const double *X, const double *Y, const double *Z, int omh);
int media(int NUM, const int *X);
void modoshi(double V[3], double R, double G, double B, double *X, double *Y,
             double *Z);
void kaiten(double V[3], double X, double Y, double Z, double *R, double *G,
            double *B);
void SobelFilterHorizontal(int hsize, int vsize, const int *IN, int *OUT);
void SobelFilterVertical(int hsize, int vsize, const int *IN, int *OUT);
void PrewitAbsoluteFilterVertical(int hsize, int vsize, int et, const int *IN,
                                  int *OUT, int *EDGE);
void PrewitAbsoluteFilterHorizontal(int hsize, int vsize, int et, const int *IN,
                                    int *OUT, int *EDGE);
int DetectEdge0(int hsize, int vsize, int et, const int *HHEDGER,
                const int *VVEDGER, const int *HHEDGEG, const int *VVEDGEG,
                const int *HHEDGEB, const int *VVEDGEB, int *OUT);
int DetectEdge1(int hsize, int vsize, int et, const int *HHEDGER,
                const int *VVEDGER, const int *HHEDGEG, const int *VVEDGEG,
                const int *HHEDGEB, const int *VVEDGEB, int *OUT);
int DetectEdge2(int hsize, int vsize, const int *EDGER, const int *VEDGER,
                const int *EDGEG, const int *VEDGEG, const int *EDGEB,
                const int *VEDGEB, int *OUT);
int DetectEdge3(int hsize, int vsize, int et, const int *HHEDGER,
                const int *VVEDGER, const int *HHEDGEG, const int *VVEDGEG,
                const int *HHEDGEB, const int *VVEDGEB, int *OUT);
int DetectEdge4(int hsize, int vsize, const int *HHEDGER, const int *VVEDGER,
                const int *HHEDGEG, const int *VVEDGEG, const int *HHEDGEB,
                const int *VVEDGEB, int *OUT, double *EDGERASISAY);
int DetectEdge5(int hsize, int vsize, int et, const int *HHEDGER,
                const int *VVEDGER, const int *HHEDGEG, const int *VVEDGEG,
                const int *HHEDGEB, const int *VVEDGEB, int *OUT);
int DetectEdge6(int hsize, int vsize, int et, const int *HHEDGER,
                const int *VVEDGER, const int *HHEDGEG, const int *VVEDGEG,
                const int *HHEDGEB, const int *VVEDGEB, int *OUT);
std::tuple<double, double, double> SumEquals(int SIZE, const int *INDEX,
                                             int INDEX_NUM, const double *v1,
                                             const double *v2,
                                             const double *v3);
std::tuple<double, double, double> SumNotEquals(int SIZE, const int *INDEX,
                                                int INDEX_NUM, const double *v1,
                                                const double *v2,
                                                const double *v3);
int FindMinimumDistanceIndex(int SIZE, double v1, double v2, double v3,
                             const uint8_t *a1, const uint8_t *a2,
                             const uint8_t *a3);
int FindMinimumDistanceIndex(int SIZE, double v1, double v2, double v3,
                             const double *a1, const double *a2,
                             const double *a3);
void Dithering(int IMAGE_SIZE, const double *YIN, const double *UIN,
               const double *VIN, int PALET_SIZE, const double *Y_JYUSHIN,
               const double *U_JYUSHIN, const double *V_JYUSHIN, uint8_t *OUT);

struct palet //パレット構造体
{
  int U_R;
  int U_G;
  int U_B;
  int IG_R;
  int IG_G;
  int IG_B;
  int INDEXNUM;
  int INDEXNO;
  double MAXDISTANCE;
  double IGENMAX;
};

struct palet PT[2][IROSUU];

// U から T に変換して複製
template <typename T, typename U>
std::vector<T> DuplicateCast(size_t SIZE, const U *IN) {
  std::vector<T> OUT(SIZE);
  for (size_t i = 0; i < SIZE; i++) {
    OUT[i] = static_cast<T>(IN[i]);
  }

  return OUT;
}

void MedianCut(int hsize, int vsize, unsigned char *RIN, unsigned char *GIN,
               unsigned char *BIN, unsigned char *PALETGAZOU,
               unsigned char REDUCE_R[256], unsigned char REDUCE_G[256],
               unsigned char REDUCE_B[256], int vvv, int et, int edon, int kmon,
               double rmult, double gmult, double bmult, double gamma, int omh,
               int hasyori, int dither, int cs, double per, int div, double pw,
               double pw2, double IRO_THIETA, double aaaa, double rgamma,
               double bgamma, double ygamma, int bun, int lpf, double nbai,
               double lpfk) // unsigned char *ROUT, unsigned char *GOUT,
                            // unsigned char *BOUT)
{
  int DIVIDENUM;
  // int RDIV;
  // int GDIV;
  // int BDIV;
  int m, p;
  int MEN;
  int NMEN;
  int NUM;
  // unsigned char *PALETGAZOU;
  // unsigned char REDUCE_R[IROSUU];
  // unsigned char REDUCE_G[IROSUU];
  // unsigned char REDUCE_B[IROSUU];
  double SUM_R[IROSUU];
  double SUM_G[IROSUU];
  double SUM_B[IROSUU];
  double HEIKIN_R[IROSUU];
  double HEIKIN_G[IROSUU];
  double HEIKIN_B[IROSUU];
  int i1, ind;
  double IGENMAX;
  int Y;
  int *INDEX;
  double MAXINDEXNUM;
  double PI = atan(1.0) * 4.0;
  // RGBをYUVに変換
  double *YIN = (double *)malloc(sizeof(double) * hsize * vsize);
  double *UIN = (double *)malloc(sizeof(double) * hsize * vsize);
  double *VIN = (double *)malloc(sizeof(double) * hsize * vsize);
  double *RrIN = (double *)malloc(sizeof(double) * hsize * vsize);
  double *GgIN = (double *)malloc(sizeof(double) * hsize * vsize);
  double *BbIN = (double *)malloc(sizeof(double) * hsize * vsize);
  double *RRR33 = (double *)malloc(sizeof(double) * hsize * vsize);

  if (lpf == 1) {
    //    while(1);
    //        double lpfk = 0.75;
    double koi;
    koi = (1.0 - lpfk) / 2.0;
    for (int i = 0; i < vsize; i++) {
      for (int j = 0; j < hsize; j++) {
        if ((j >= 1) && (j <= (hsize - 2))) {
          YIN[j * vsize + i] = lpfk * (double)RIN[j * vsize + i] +
                               koi * (double)RIN[(j - 1) * vsize + i] +
                               koi * (double)RIN[(j + 1) * vsize + i];
          UIN[j * vsize + i] = lpfk * (double)GIN[j * vsize + i] +
                               koi * (double)GIN[(j - 1) * vsize + i] +
                               koi * (double)GIN[(j + 1) * vsize + i];
          VIN[j * vsize + i] = lpfk * (double)BIN[j * vsize + i] +
                               koi * (double)BIN[(j - 1) * vsize + i] +
                               koi * (double)BIN[(j + 1) * vsize + i];
        }
      }
      // if(i==1){while(1);}
    }
    //    while(1);
    for (int i = 0; i < hsize; i++) {
      for (int j = 0; j < vsize; j++) {
        if ((j >= 1) && (j <= (vsize - 2))) {
          RIN[i * vsize + j] =
              (unsigned char)(lpfk * YIN[i * vsize + j] +
                              koi * YIN[i * vsize + (j - 1)] +
                              koi * YIN[i * vsize + (j + 1)] + 0.5);
          GIN[i * vsize + j] =
              (unsigned char)(lpfk * UIN[i * vsize + j] +
                              koi * UIN[i * vsize + (j - 1)] +
                              koi * UIN[i * vsize + (j + 1)] + 0.5);
          BIN[i * vsize + j] =
              (unsigned char)(lpfk * VIN[i * vsize + j] +
                              koi * VIN[i * vsize + (j - 1)] +
                              koi * VIN[i * vsize + (j + 1)] + 0.5);
        }
      }
    }
  } // else {
  //    for(i=0;i<hsize*vsize;i++){
  //        RrIN[i] =
  //        (double)RIN[i];//+0.5*YIN[i*vsize+(j-1)]+0.25*YIN[i*vsize+(j+1)];
  //        GgIN[i] =
  //        (double)GIN[i];//+0.5*UIN[i*vsize+(j-1)]+0.25*UIN[i*vsize+(j+1)];
  //        BbIN[i] =
  //        (double)BIN[i];//+0.5*VIN[i*vsize+(j-1)]+0.25*VIN[i*vsize+(j+1)];
  //    }

  //}
  for (int i = 0; i < hsize * vsize; i++) {
    if (RIN[i] < hasyori) {
      RIN[i] = hasyori;
    }
    if (GIN[i] < hasyori) {
      GIN[i] = hasyori;
    }
    if (BIN[i] < hasyori) {
      BIN[i] = hasyori;
    }
    // if(RIN[i] > MHASYORI){
    // RIN[i] = MHASYORI;
    // }
    // if(GIN[i] > MHASYORI){
    // GIN[i] = MHASYORI;
    // }
    // if(BIN[i] > MHASYORI){
    // BIN[i] = MHASYORI;
    // }
  }
  for (int i = 0; i < hsize * vsize; i++) {
    RrIN[i] = 255.0 * pow((double)RIN[i] / 255.0, rgamma);
    GgIN[i] = 255.0 * pow((double)GIN[i] / 255.0, gamma);
    BbIN[i] = 255.0 * pow((double)BIN[i] / 255.0, bgamma);
  }

  IRO_THIETA *= PI / 180;
  if (cs == 0) {
    for (int i = 0; i < hsize * vsize; i++) {
      YIN[i] = *(
          GgIN +
          i); //(0.29891*(float)(RIN[i])+0.58661*(float)(GIN[i])+0.11448*(float)(BIN[i])
              ///*+ 0.5*/);
      UIN[i] = *(
          BbIN +
          i); //((-0.16874*(float)(RIN[i])-0.33126*(float)(GIN[i])+0.50000*(float)(BIN[i]))
              //+ 128.0);
      VIN[i] = *(
          RrIN +
          i); //(0.50000*(float)(RIN[i])-0.41869*(float)(GIN[i])-0.08131*(float)(BIN[i])
              //+ 128.0);
    }
  } else if (cs == 1) {
    for (int i = 0; i < hsize * vsize; i++) {
      RGBtoLAB(
          RrIN[i], GgIN[i], BbIN[i], VIN + i, YIN + i,
          UIN +
              i); //(0.29891*(float)(RIN[i])+0.58661*(float)(GIN[i])+0.11448*(float)(BIN[i])
                  //+ 0.5);
      // UIN[i] =
      // BIN[i];//((-0.16874*(float)(RIN[i])-0.33126*(float)(GIN[i])+0.50000*(float)(BIN[i]))
      //+ 128.5);
      // VIN[i] =
      // RIN[i];//(0.50000*(float)(RIN[i])-0.41869*(float)(GIN[i])-0.08131*(float)(BIN[i])
      //+ 128.5);
    }
  } else if (cs == 2) {
    for (int i = 0; i < hsize * vsize; i++) {
      RGBtoLuv(
          RrIN[i], GgIN[i], BbIN[i], VIN + i, YIN + i,
          UIN +
              i); //(0.29891*(float)(RIN[i])+0.58661*(float)(GIN[i])+0.11448*(float)(BIN[i])
                  //+ 0.5);
      // UIN[i] =
      // BIN[i];//((-0.16874*(float)(RIN[i])-0.33126*(float)(GIN[i])+0.50000*(float)(BIN[i]))
      //+ 128.5);
      // VIN[i] =
      // RIN[i];//(0.50000*(float)(RIN[i])-0.41869*(float)(GIN[i])-0.08131*(float)(BIN[i])
      //+ 128.5);
    }
  } else if (cs == 3 || cs == 9) {
    for (int i = 0; i < hsize * vsize; i++) {
      YIN[i] = (0.29891 * (float)(RrIN[i]) + 0.58661 * (float)(GgIN[i]) +
                0.11448 * (float)(BbIN[i]) /*+ 0.5*/);
      UIN[i] = ((-0.16874 * (float)(RrIN[i]) - 0.33126 * (float)(GgIN[i]) +
                 0.50000 * (float)(BbIN[i])) +
                128.0);
      VIN[i] = (0.50000 * (float)(RrIN[i]) - 0.41869 * (float)(GgIN[i]) -
                0.08131 * (float)(BbIN[i]) + 128.0);
    }
  } else if (cs == 4) {
    for (int i = 0; i < hsize * vsize; i++) {
      if ((double)(RIN[i]) / 255.0 >= 0.117647058823529411764705882352941) {
        RrIN[i] = 255.0 * pow((double)(RrIN[i]) / 255.0, 3.0 / 2.2);
      } else {
        RrIN[i] = 255.0 *
                  pow((double)(RrIN[i]) / 255.0 /**0.1429136476*/, 1.0 / 2.2) *
                  0.1429136476;
      }
      if ((double)(GIN[i]) / 255.0 >= 0.117647058823529411764705882352941) {
        GgIN[i] = 255.0 * pow((double)(GgIN[i]) / 255.0, 3.0 / 2.2);
      } else {
        GgIN[i] = 255.0 *
                  pow((double)(GgIN[i]) / 255.0 /**0.1429136476*/, 1.0 / 2.2) *
                  0.1429136476;
      }
      if ((double)(BIN[i]) / 255.0 >= 0.117647058823529411764705882352941) {
        BbIN[i] = 255.0 * pow((double)(BbIN[i]) / 255.0, 3.0 / 2.2);
      } else {
        BbIN[i] = 255.0 *
                  pow((double)(BbIN[i]) / 255.0 /*0.1429136476*/, 1.0 / 2.2) *
                  0.1429136476;
      }

      YIN[i] = (0.2989 * (float)(RrIN[i]) + 0.5866 * (float)(GgIN[i]) +
                0.1145 * (float)(BbIN[i]) /*+0.5*/);
      UIN[i] = ((-0.1350 * (float)(RrIN[i]) - 0.2650 * (float)(GgIN[i]) +
                 0.40000 * (float)(BbIN[i])));
      VIN[i] = (0.4000 * (float)(RrIN[i]) - 0.3346 * (float)(GgIN[i]) -
                0.0653 * (float)(BbIN[i]));
      //                n1= (double)r*0.2989+(double)g*0.5866+(double)b*0.1145;
      //                n2=-(double)r*0.1350-(double)g*0.2650+(double)b*0.4000;
      //                n3= (double)r*0.4000-(double)g*0.3346-(double)b*0.0653;
    }
  } else if (cs == 5) {
    for (int i = 0; i < hsize * vsize; i++) {
      RrIN[i] = (255.0 / (1.0 + exp(-aaaa * ((double)RrIN[i] - 127.5))) -
                 255.0 / (1.0 + exp(-aaaa * (-127.5)))) *
                255.0 / (255.0 - 2.0 * 255.0 / (1.0 + exp(-aaaa * (-127.5))));
      GgIN[i] = (255.0 / (1.0 + exp(-aaaa * ((double)GgIN[i] - 127.5))) -
                 255.0 / (1.0 + exp(-aaaa * (-127.5)))) *
                255.0 / (255.0 - 2.0 * 255.0 / (1.0 + exp(-aaaa * (-127.5))));
      BbIN[i] = (255.0 / (1.0 + exp(-aaaa * ((double)BbIN[i] - 127.5))) -
                 255.0 / (1.0 + exp(-aaaa * (-127.5)))) *
                255.0 / (255.0 - 2.0 * 255.0 / (1.0 + exp(-aaaa * (-127.5))));
      YIN[i] = (0.29891 * (float)(RrIN[i]) + 0.58661 * (float)(GgIN[i]) +
                0.11448 * (float)(BbIN[i]) /*+ 0.5*/);
      UIN[i] = ((-0.16874 * (float)(RrIN[i]) - 0.33126 * (float)(GgIN[i]) +
                 0.50000 * (float)(BbIN[i])) +
                128.0);
      VIN[i] = (0.50000 * (float)(RrIN[i]) - 0.41869 * (float)(GgIN[i]) -
                0.08131 * (float)(BbIN[i]) + 128.0);
    }
  } else if (cs == 6) {
    for (int i = 0; i < hsize * vsize; i++) {
      RrIN[i] = (255.0 / (1.0 + exp(-aaaa * ((double)RrIN[i] - 127.5))) -
                 255.0 / (1.0 + exp(-aaaa * (-127.5)))) *
                255.0 / (255.0 - 2.0 * 255.0 / (1.0 + exp(-aaaa * (-127.5))));
      GgIN[i] = (255.0 / (1.0 + exp(-aaaa * ((double)GgIN[i] - 127.5))) -
                 255.0 / (1.0 + exp(-aaaa * (-127.5)))) *
                255.0 / (255.0 - 2.0 * 255.0 / (1.0 + exp(-aaaa * (-127.5))));
      BbIN[i] = (255.0 / (1.0 + exp(-aaaa * ((double)BbIN[i] - 127.5))) -
                 255.0 / (1.0 + exp(-aaaa * (-127.5)))) *
                255.0 / (255.0 - 2.0 * 255.0 / (1.0 + exp(-aaaa * (-127.5))));

      YIN[i] = *(
          GgIN +
          i); //(0.29891*(float)(RIN[i])+0.58661*(float)(GIN[i])+0.11448*(float)(BIN[i])
              ///*+ 0.5*/);
      UIN[i] = *(
          BbIN +
          i); //((-0.16874*(float)(RIN[i])-0.33126*(float)(GIN[i])+0.50000*(float)(BIN[i]))
              //+ 128.0);
      VIN[i] = *(
          RrIN +
          i); //(0.50000*(float)(RIN[i])-0.41869*(float)(GIN[i])-0.08131*(float)(BIN[i])
              //+ 128.0);
    }
  } else if (cs == 7) {
    for (int i = 0; i < hsize * vsize; i++) {
      if ((double)(RIN[i]) / 255.0 >= 0.1176470588) {
        RrIN[i] = 255.0 * pow((double)(RrIN[i]) / 255.0, 3.0 / 2.2);
      } else {
        RrIN[i] = 255.0 *
                  pow((double)(RrIN[i]) / 255.0 /**0.1429136476*/, 1.0 / 2.2) *
                  0.1429136476;
      }
      if ((double)(GIN[i]) / 255.0 >= 0.1176470588) {
        GgIN[i] = 255.0 * pow((double)(GgIN[i]) / 255.0, 3.0 / 2.2);
      } else {
        GgIN[i] = 255.0 *
                  pow((double)(GgIN[i]) / 255.0 /**0.1429136476*/, 1.0 / 2.2) *
                  0.1429136476;
      }
      if ((double)(BIN[i]) / 255.0 >= 0.1176470588) {
        BbIN[i] = 255.0 * pow((double)(BbIN[i]) / 255.0, 3.0 / 2.2);
      } else {
        BbIN[i] = 255.0 *
                  pow((double)(BbIN[i]) / 255.0 /**0.1429136476*/, 1.0 / 2.2) *
                  0.1429136476;
      }

      YIN[i] = GgIN[i];
      UIN[i] = BbIN
          [i]; //((-0.1350*(float)(RrIN[i])-0.2650*(float)(GgIN[i])+0.40000*(float)(BbIN[i])));
      VIN[i] = RrIN
          [i]; //(0.4000*(float)(RrIN[i])-0.3346*(float)(GgIN[i])-0.0653*(float)(BbIN[i]));
    }
  } else if (cs == 8) {
    for (int i = 0; i < hsize * vsize; i++) {
      if ((double)(RIN[i]) / 255.0 >= 0.1176470588) {
        RrIN[i] = 255.0 * pow((double)(RrIN[i]) / 255.0, 3.0 / 2.2);
      } else {
        RrIN[i] = 255.0 *
                  pow((double)(RrIN[i]) / 255.0 /**0.1429136476*/, 1.0 / 2.2) *
                  0.1429136476;
      }
      if ((double)(GIN[i]) / 255.0 >= 0.1176470588) {
        GgIN[i] = 255.0 * pow((double)(GgIN[i]) / 255.0, 3.0 / 2.2);
      } else {
        GgIN[i] = 255.0 *
                  pow((double)(GgIN[i]) / 255.0 /**0.1429136476*/, 1.0 / 2.2) *
                  0.1429136476;
      }
      if ((double)(BIN[i]) / 255.0 >= 0.1176470588) {
        BbIN[i] = 255.0 * pow((double)(BbIN[i]) / 255.0, 3.0 / 2.2);
      } else {
        BbIN[i] = 255.0 *
                  pow((double)(BbIN[i]) / 255.0 /**0.1429136476*/, 1.0 / 2.2) *
                  0.1429136476;
      }

      // for(i=0;i<hsize*vsize;i++){
      YIN[i] = (0.29891 * (float)(RrIN[i]) + 0.58661 * (float)(GgIN[i]) +
                0.11448 * (float)(BbIN[i]) /*+ 0.5*/);
      UIN[i] = ((-0.16874 * (float)(RrIN[i]) - 0.33126 * (float)(GgIN[i]) +
                 0.50000 * (float)(BbIN[i])) +
                128.0);
      VIN[i] = (0.50000 * (float)(RrIN[i]) - 0.41869 * (float)(GgIN[i]) -
                0.08131 * (float)(BbIN[i]) + 128.0);
      //}

      // YIN[i] = GgIN[i];
      // UIN[i] =
      // BbIN[i];//((-0.1350*(float)(RrIN[i])-0.2650*(float)(GgIN[i])+0.40000*(float)(BbIN[i])));
      // VIN[i] =
      // RrIN[i];//(0.4000*(float)(RrIN[i])-0.3346*(float)(GgIN[i])-0.0653*(float)(BbIN[i]));
    }
  }
  double *KARIU, *KARIV;
  KARIU = (double *)malloc(sizeof(double) * hsize * vsize);
  KARIV = (double *)malloc(sizeof(double) * hsize * vsize);
  if (cs == 4 || cs == 9) {
    for (int i = 0; i < hsize * vsize; i++) {
      KARIU[i] = UIN[i] * /*cos(IRO_THIETA)*/ 0.941996495 +
                 VIN[i] * /*sin(IRO_THIETA)*/ 0.13632709;
      KARIV[i] = UIN[i] * /*sin(IRO_THIETA)*/ 0.118202579 +
                 VIN[i] * /*cos(IRO_THIETA)*/ 1.201426141;
      UIN[i] = KARIU[i];
      VIN[i] = KARIV[i];
    }
  } else {
    for (int i = 0; i < hsize * vsize; i++) {
      KARIU[i] = UIN[i] * cos(IRO_THIETA) - VIN[i] * sin(IRO_THIETA);
      KARIV[i] = UIN[i] * sin(IRO_THIETA) + VIN[i] * cos(IRO_THIETA);
      UIN[i] = KARIU[i];
      VIN[i] = KARIV[i];
    }
  }
  for (int i = 0; i < hsize * vsize; i++) {
    YIN[i] = 255.0 * pow(YIN[i] / 255.0, ygamma);
  }

  for (int i = 0; i < hsize * vsize; i++) {
    YIN[i] *= gmult;
  }
  for (int i = 0; i < hsize * vsize; i++) {
    VIN[i] *= rmult;
  }
  // for(i=0;i<hsize*vsize;i++){
  // VIN[i] = 255.0*rmult*pow(VIN[i]/(rmult*255.0),rgamma);
  //}
  for (int i = 0; i < hsize * vsize; i++) {
    UIN[i] *= bmult;
  }
  // for(i=0;i<hsize*vsize;i++){
  // UIN[i] = 255.0*bmult*pow(UIN[i]/(bmult*255.0),bgamma);
  //}
  // YUVは、0-255 たまに負や255を超える場合が発生するかも、あとで検証
  // debug start
  // for(i=0;i<hsize*vsize;i++){
  // fprintf(stderr," %d %d %d\n",YIN[i],UIN[i],VIN[i]);
  // }
  // while(1);
  // debug end
  /****************************************/
  /* */
  /* edge detect */
  /* */
  /****************************************/

  double EDGERASISAYD;
  std::vector<int> EDGER(hsize * vsize);
  std::vector<int> EDGEG(hsize * vsize);
  std::vector<int> EDGEB(hsize * vsize);
  std::vector<int> VEDGER(hsize * vsize);
  std::vector<int> VEDGEG(hsize * vsize);
  std::vector<int> VEDGEB(hsize * vsize);
  // unsigned char から int に変換して複製
  std::vector<int> IRIN = DuplicateCast<int>(hsize * vsize, &RIN[0]);
  std::vector<int> IGIN = DuplicateCast<int>(hsize * vsize, &GIN[0]);
  std::vector<int> IBIN = DuplicateCast<int>(hsize * vsize, &BIN[0]);

  std::vector<int> HHEDGER(hsize * vsize);
  std::vector<int> HHEDGEG(hsize * vsize);
  std::vector<int> HHEDGEB(hsize * vsize);
  std::vector<int> VVEDGER(hsize * vsize);
  std::vector<int> VVEDGEG(hsize * vsize);
  std::vector<int> VVEDGEB(hsize * vsize);

  std::vector<int> index3(hsize * vsize);

  int EDGE_GASOSUU = 0;

  std::vector<double> EDGERASISAY;

  if (vvv == 3 || vvv == 4 || vvv == 5) {
    SobelFilterHorizontal(hsize, vsize, &IRIN[0], &HHEDGER[0]);
    SobelFilterVertical(hsize, vsize, &IRIN[0], &VVEDGER[0]);

    SobelFilterHorizontal(hsize, vsize, &IGIN[0], &HHEDGEG[0]);
    SobelFilterVertical(hsize, vsize, &IGIN[0], &VVEDGEG[0]);

    SobelFilterHorizontal(hsize, vsize, &IBIN[0], &HHEDGEB[0]);
    SobelFilterVertical(hsize, vsize, &IBIN[0], &VVEDGEB[0]);

    if (vvv == 3) {
      EDGE_GASOSUU =
          DetectEdge3(hsize, vsize, et, &HHEDGER[0], &VVEDGER[0], &HHEDGEG[0],
                      &VVEDGEG[0], &HHEDGEB[0], &VVEDGEB[0], &index3[0]);
    } else if (vvv == 4) {
      EDGERASISAY = std::vector<double>(hsize * vsize);

      EDGE_GASOSUU = DetectEdge4(hsize, vsize, &HHEDGER[0], &VVEDGER[0],
                                 &HHEDGEG[0], &VVEDGEG[0], &HHEDGEB[0],
                                 &VVEDGEB[0], &index3[0], &EDGERASISAY[0]);
    } else if (vvv == 5) {
      EDGE_GASOSUU =
          DetectEdge5(hsize, vsize, et, &HHEDGER[0], &VVEDGER[0], &HHEDGEG[0],
                      &VVEDGEG[0], &HHEDGEB[0], &VVEDGEB[0], &index3[0]);
    }
  }

  if ((vvv >= 0 && vvv <= 2) || vvv == 6) {
    PrewitAbsoluteFilterVertical(hsize, vsize, et, &IRIN[0], &HHEDGER[0],
                                 &EDGER[0]);
    PrewitAbsoluteFilterHorizontal(hsize, vsize, et, &IRIN[0], &VVEDGER[0],
                                   &VEDGER[0]);

    PrewitAbsoluteFilterVertical(hsize, vsize, et, &IGIN[0], &HHEDGEG[0],
                                 &EDGEG[0]);
    PrewitAbsoluteFilterHorizontal(hsize, vsize, et, &IGIN[0], &VVEDGEG[0],
                                   &VEDGEG[0]);

    PrewitAbsoluteFilterVertical(hsize, vsize, et, &IBIN[0], &HHEDGEB[0],
                                 &EDGEB[0]);
    PrewitAbsoluteFilterHorizontal(hsize, vsize, et, &IBIN[0], &VVEDGEB[0],
                                   &VEDGEB[0]);

    if (vvv == 0) {
      EDGE_GASOSUU =
          DetectEdge0(hsize, vsize, et, &HHEDGER[0], &VVEDGER[0], &HHEDGEG[0],
                      &VVEDGEG[0], &HHEDGEB[0], &VVEDGEB[0], &index3[0]);
    } else if (vvv == 1) {
      EDGE_GASOSUU =
          DetectEdge1(hsize, vsize, et, &HHEDGER[0], &VVEDGER[0], &HHEDGEG[0],
                      &VVEDGEG[0], &HHEDGEB[0], &VVEDGEB[0], &index3[0]);
    } else if (vvv == 2) {
      EDGE_GASOSUU = DetectEdge2(hsize, vsize, &EDGER[0], &VEDGER[0], &EDGEG[0],
                                 &VEDGEG[0], &EDGEB[0], &VEDGEB[0], &index3[0]);
    } else if (vvv == 6) {
      EDGE_GASOSUU =
          DetectEdge6(hsize, vsize, et, &HHEDGER[0], &VVEDGER[0], &HHEDGEG[0],
                      &VVEDGEG[0], &HHEDGEB[0], &VVEDGEB[0], &index3[0]);
    }
  } // cnt if end

  // debug start]
  // if(vvv == 4){
  //    printf("EDGE rasisa mode\n");
  //} else {
  printf("EDGE gasosuu=%d\n", EDGE_GASOSUU);
  //}
  // debug end

  //    if(0/*omh==3*/){
  /*        int tempju;double tempuj;
          int *RGB;double *CONVY,*CONVU,*CONVV;
          RGB =(int*)malloc(sizeof(int)*hsize*vsize);
          CONVY =(double*)malloc(sizeof(double)*hsize*vsize);
          CONVU =(double*)malloc(sizeof(double)*hsize*vsize);
          CONVV =(double*)malloc(sizeof(double)*hsize*vsize);
          for(i=0,m=0;i<hsize*vsize;i++){
              if(index3[i] !=-1){
                  RGB[m] = (((int)RIN[i]) << 16) | (((int)GIN[i]) << 8) |
     ((int)BIN[i]); CONVY[m] = YIN[i]; CONVU[m] = UIN[i]; CONVV[m] = VIN[i];
                  m++;
              }
          }
          int *IROIRO;
          IROIRO = (int*)malloc(sizeof(int)*256*256*256);
          for(i=0;i<256*256*256;i++){
              IROIRO[i]=0;
          }
          for(i=0;i<hsize*vsize;i++){
          IROIRO[RGB[i]]++;
          }
          for(j=0,i=0;i<256*256*256;i++){
              if(IROIRO[i] !=0){
                  j++;
              }
          }
          printf("IROSUUSEIKAKU=%d\n",j);//while(1);
          int *IROINDEX;
          IROINDEX = (int*)malloc(sizeof(int)*hsize*vsize);
          for(i=0,m=0;i<hsize*vsize;i++){
              if(index3[i] !=-1){
                  IROINDEX[m] = i;
                  m++;
              }
          }
          for(i=0;i<m-1;i++){
              for(j=m-1;j>i;j--){
                  if(RGB[j-1] > RGB[j]){
                      tempju = RGB[j-1];
                      RGB[j-1] = RGB[j];
                      RGB[j] = tempju;
                      tempuj = CONVY[j-1];
                      CONVY[j-1] = CONVY[j];
                      CONVY[j] = tempuj;
                      tempuj = CONVU[j-1];
                      CONVU[j-1] = CONVU[j];
                      CONVU[j] = tempuj;
                      tempuj = CONVV[j-1];
                      CONVV[j-1] = CONVV[j];
                      CONVV[j] = tempuj;
                      tempju = IROINDEX[j-1];
                      IROINDEX[j-1] = IROINDEX[j];
                      IROINDEX[j] = tempju;
                  }
              }
          }
          for(i=0;i<m;i++){
              printf("%d %x\n",i,RGB[i]);
          }
          //while(1);
          //int samecolornum;
          int FLAGY;
          iro->RGB = RGB[0];
          iro->CONVY = CONVY[0];
          iro->CONVU = CONVU[0];
          iro->CONVV = CONVV[0];
          iro->INDEX = IROINDEX[0];
          i=0;j=1;k=1;
          while(i<m-1){
              i++;
              if(RGB[i-1] != RGB[i]){
                  (iro+k)->RGB = RGB[i];
                  (iro+k)->CONVY = CONVY[i];
                  (iro+k)->CONVU = CONVU[i];
                  (iro+k)->CONVV = CONVV[i];
                  (iro+k-1)->CNT = j;
                  (iro+k)->INDEX = IROINDEX[i];
                  j=1;k++;FLAGY=0;
              }else{
                  j++;FLAGY=1;
              }

          }
          if(FLAGY==0){
              (iro+k)->RGB = RGB[i];
              (iro+k)->CONVY = CONVY[i];
              (iro+k)->CONVU = CONVU[i];
              (iro+k)->CONVV = CONVV[i];
              (iro+k-1)->CNT = 1;
              (iro+k)->INDEX = IROINDEX[i];
          } else if(FLAGY==1){
              (iro+k-1)->CNT = j;
          }
          int IROSUUU;
          IROSUUU = k;
          printf("\n");//k ni iro no syuruisuu ga kioku sareteiru.
          for(i=0;i<k;i++){
              printf("%x %d\n",(iro+i)->RGB,(iro+i)->CNT);
          }//while(1);
          printf("IROSUU=%d\n",IROSUUU);
          //while(1);
      }

  */

  DIVIDENUM = 1;
  MEN = 0;
  // RDIV = RD;
  // GDIV = GD;
  // BDIV = BD;
  if (MEN == 1) {
    NMEN = 0;
  } else {
    NMEN = 1;
  }
  // MEN = 0,NMEN = 1
  auto U_RGB1 =
      SumNotEquals(hsize * vsize, &index3[0], -1, &YIN[0], &UIN[0], &VIN[0]);
  std::get<0>(U_RGB1) /= hsize * vsize - EDGE_GASOSUU;
  std::get<1>(U_RGB1) /= hsize * vsize - EDGE_GASOSUU;
  std::get<2>(U_RGB1) /= hsize * vsize - EDGE_GASOSUU;
  double U_R1 = std::get<0>(U_RGB1);
  double U_G1 = std::get<1>(U_RGB1);
  double U_B1 = std::get<2>(U_RGB1);
  // U_R = (U_R>>RDIV);
  // U_G = (U_G>>GDIV);
  // U_B = (U_B>>BDIV);
  // U_R,U_G,U_Bは平均
  int ct = 1000;
  double eps = 1.0e-10;
  int n = 3;
  double **A = new double *[n];
  double **A1 = new double *[n];
  double **A2 = new double *[n];
  double **X1 = new double *[n];
  double **X2 = new double *[n];
  for (i1 = 0; i1 < n; i1++) {
    A[i1] = new double[n];
    A1[i1] = new double[n];
    A2[i1] = new double[n];
    X1[i1] = new double[n];
    X2[i1] = new double[n];
  }
  JacobiSetupNotEquals(hsize * vsize, &index3[0], -1, &YIN[0], &UIN[0], &VIN[0],
                       U_RGB1, A);
  ind = Jacobi(3, ct, eps, A, A1, A2, X1, X2);

  if (ind > 0) {
    fprintf(stderr, "Jacobi syuusoku simasen!!\n");
    exit(0);
  }
  // if((A1[0][0] == A1[1][1]) || (A1[1][1] == A1[2][2]) || (A1[0][0] ==
  // A1[2][2]) ){ fprintf(stderr,"Igen Value duplication!!\n"); exit(0);
  // }
  // debug start
  // for(i=0;i<3;i++){
  // for(j=0;j<3;j++){
  // fprintf(stderr,"%f ",A1[i][j]);
  // }
  // fprintf(stderr,"\n");
  // }
  // debug end
  IGENMAX = -999999999999999999.9;
  for (int i = 0; i < 3; i++) {
    if (A1[i][i] > IGENMAX) {
      IGENMAX = A1[i][i];
      Y = i;
    }
  }
  // debug start
  // fprintf(stderr,"Y=%d\n",Y);
  // fprintf(stderr,"IGENVEC=%f %f %f\n",X1[0][Y],X1[1][Y],X1[2][Y]);
  // while(1);
  // debug end
  // Yは最大固有値の番号

  //(X1[0][Y],X1[1][Y],X1[2][Y]) 最大分散の固有ベクトル
  // ohtsu tuika start
  double V[3];
  double MAXX, MINN;
  int THRESH = 0;
  double DTHRESH;
  std::vector<double> RR(hsize * vsize);
  std::vector<double> GG(hsize * vsize);
  std::vector<double> BB(hsize * vsize);
  double *XXX = (double *)malloc(sizeof(double));
  double *YYY = (double *)malloc(sizeof(double));
  double *ZZZ = (double *)malloc(sizeof(double));
  double *RRR = (double *)malloc(sizeof(double) * hsize * vsize);
  double *GGG = (double *)malloc(sizeof(double) * hsize * vsize);
  double *BBB = (double *)malloc(sizeof(double) * hsize * vsize);
  for (int i = 0; i < hsize * vsize; i++) {
    if (index3[i] != -1) {
      RR[i] = YIN[i] - U_R1;
      GG[i] = UIN[i] - U_G1;
      BB[i] = VIN[i] - U_B1;
    }
  }
  V[0] = X1[0][Y];
  V[1] = X1[1][Y];
  V[2] = X1[2][Y];
  int *iRRRR = (int *)malloc(sizeof(int) * (hsize * vsize));
  for (int i = 0; i < hsize * vsize; i++) {
    kaiten(V, RR[i], GG[i], BB[i], RRR + i, GGG + i, BBB + i);
  }
  for (int i = 0, l = 0; i < hsize * vsize; i++) {
    if (index3[i] != -1) {
      iRRRR[l] = (int)(RRR[i] + 0.5);
      l++;
    }
  }

  double *RRRRRR = nullptr;
  double *GGGGGG = nullptr;
  double *BBBBBB = nullptr;
  int GASOSUU = hsize * vsize - EDGE_GASOSUU;
  if (omh == 3 || omh == 4) {

    RRRRRR = (double *)malloc(sizeof(double) * hsize * vsize);
    GGGGGG = (double *)malloc(sizeof(double) * hsize * vsize);
    BBBBBB = (double *)malloc(sizeof(double) * hsize * vsize);
    for (int i = 0, l = 0; i < hsize * vsize; i++) {
      if (index3[i] != -1) {
        RRRRRR[l] = RRR[i];
        GGGGGG[l] = GGG[i];
        BBBBBB[l] = BBB[i];
        l++;
      }
    }
    THRESH = ohtsu2(GASOSUU, RRRRRR, GGGGGG, BBBBBB, omh);
    printf("THRESH OHTSU2=%d\n", THRESH);
  }
  if (omh == 0) {
    THRESH = ohtsu(GASOSUU, iRRRR);
    printf("THRESH OHTSU=%d\n", THRESH);
  } else if (omh == 1) {
    THRESH = media(GASOSUU, iRRRR);
  } else if (omh == 2) {
    THRESH = 0;
  }
  // while(1);
  DTHRESH = (double)(THRESH);
  modoshi(V, DTHRESH, 0.0, 0.0, XXX, YYY, ZZZ);

  // MEN = 0,NMEN = 1
  PT[MEN][0].INDEXNO = 0;
  PT[MEN][1].INDEXNO = 1;

  INDEX = (int *)malloc(sizeof(int) * hsize * vsize);
  // int kk;
  // kk = k;
  for (int i = 0; i < hsize * vsize; i++) {
    INDEX[i] = -1;
  }
  int k2 = 0, l2 = 0;
  for (int i = 0; i < hsize * vsize; i++) {
    if (index3[i] != -1) {
      // ZZ = X1[0][Y]*((double)(((YIN[i]))) - ((double)(U_R)+(*XXX))) +
      // X1[1][Y]*((double)(((UIN[i]))) - ((double)(U_G)+(*YYY))) +
      // X1[2][Y]*((double)(((VIN[i]))) - ((double)(U_B)+(*ZZZ))) ;
      // if(ZZ>=0.0){
      if (RRR[i] >= DTHRESH) {
        INDEX[i] = PT[MEN][0].INDEXNO;
        k2++;
      } else {
        INDEX[i] = PT[MEN][1].INDEXNO;
        l2++;
      }
    }
  }
  PT[MEN][0].INDEXNUM = k2;
  PT[MEN][1].INDEXNUM = l2;
  // maxdistance tuika start
  // 0側のmaxdistanceを求める
  auto U_RGB = SumEquals(hsize * vsize, &INDEX[0], PT[MEN][0].INDEXNO, &YIN[0],
                         &UIN[0], &VIN[0]);

  double MAXD;
  double TEMP = 0.0;
  if (PT[MEN][0].INDEXNUM != 0) {
    std::get<0>(U_RGB) /= (double)PT[MEN][0].INDEXNUM;
    std::get<1>(U_RGB) /= (double)PT[MEN][0].INDEXNUM;
    std::get<2>(U_RGB) /= (double)PT[MEN][0].INDEXNUM;
    double U_R = std::get<0>(U_RGB);
    double U_G = std::get<1>(U_RGB);
    double U_B = std::get<2>(U_RGB);
    // U_R = (U_R>>RDIV);
    // U_G = (U_G>>GDIV);
    // U_B = (U_B>>BDIV);
    if (bun == 0) {
      MAXD = 0.0;
      for (int i = 0; i < hsize * vsize; i++) {
        if (INDEX[i] == PT[MEN][0].INDEXNO) {
          TEMP = (YIN[i] - U_R) * (YIN[i] - U_R) +
                 (UIN[i] - U_G) * (UIN[i] - U_G) +
                 (VIN[i] - U_B) * (VIN[i] - U_B);
          if (MAXD < TEMP) {
            MAXD = TEMP;
          }
        }
      }
      PT[MEN][0].MAXDISTANCE = MAXD;
    } else if (bun == 1) {
      if (vvv != 4) {
        if (div == 0) {
          TEMP = 0.0;
          for (int i = 0; i < hsize * vsize; i++) {
            if (INDEX[i] == PT[MEN][0].INDEXNO) {
              TEMP += (YIN[i] - U_R) * (YIN[i] - U_R) +
                      (UIN[i] - U_G) * (UIN[i] - U_G) +
                      (VIN[i] - U_B) * (VIN[i] - U_B);
              // if( MAXD < TEMP){
              // MAXD = TEMP;
              // }
            }
          }
        }
        if (div == 1) {
          TEMP = 0.0;
          for (int i = 0; i < hsize * vsize; i++) {
            if (INDEX[i] == PT[MEN][0].INDEXNO) {
              TEMP += pow(((YIN[i] - U_R) * (YIN[i] - U_R) +
                           (UIN[i] - U_G) * (UIN[i] - U_G) +
                           (VIN[i] - U_B) * (VIN[i] - U_B)),
                          pw / 2.0);
              // if( MAXD < TEMP){
              // MAXD = TEMP;
              // }
            }
          }
        }
      } // div !=4 if end

      if (vvv == 4) {
        if (div == 0) {
          TEMP = 0.0;
          for (int i = 0; i < hsize * vsize; i++) {
            if (INDEX[i] == PT[MEN][0].INDEXNO) {
              TEMP += ((YIN[i] - U_R) * (YIN[i] - U_R) +
                       (UIN[i] - U_G) * (UIN[i] - U_G) +
                       (VIN[i] - U_B) * (VIN[i] - U_B)) *
                      (255.0 - EDGERASISAY[i]) / 255.0;
              // if( MAXD < TEMP){
              // MAXD = TEMP;
              // }
            }
          }
        }
        if (div == 1) {
          TEMP = 0.0;
          for (int i = 0; i < hsize * vsize; i++) {
            if (INDEX[i] == PT[MEN][0].INDEXNO) {
              EDGERASISAYD = EDGERASISAY[i];
              EDGERASISAYD *= nbai;
              if (EDGERASISAYD > 255.0) {
                EDGERASISAYD = 255.0;
              }
              TEMP += pow(((YIN[i] - U_R) * (YIN[i] - U_R) +
                           (UIN[i] - U_G) * (UIN[i] - U_G) +
                           (VIN[i] - U_B) * (VIN[i] - U_B)),
                          pw / 2.0) *
                      pow((255.0 - EDGERASISAYD) / 255.0, pw2);
              // if( MAXD < TEMP){
              // MAXD = TEMP;
              // }
              // printf("EDGERASISAYD = %f\n",EDGERASISAYD);
            }
          }
        }
      } // div !=4 if end

      PT[MEN][0].MAXDISTANCE = TEMP;
    } else if (bun == 2 || bun == 3 || bun == 4 || bun == 5) {
      JacobiSetupEquals(hsize * vsize, &INDEX[0], PT[MEN][0].INDEXNO, &YIN[0],
                        &UIN[0], &VIN[0], U_RGB, A);
      ind = Jacobi(3, ct, eps, A, A1, A2, X1, X2);

      if (ind > 0) {
        fprintf(stderr, "Jacobi syuusoku simasen!!\n");
        exit(0);
      }
      // if((A1[0][0] == A1[1][1]) || (A1[1][1] == A1[2][2]) || (A1[0][0] ==
      // A1[2][2]) ){ fprintf(stderr,"Igen Value duplication!!\n"); exit(0);
      // }
      // debug start
      // for(i=0;i<3;i++){
      // for(j=0;j<3;j++){
      // fprintf(stderr,"%f ",A1[i][j]);
      // }
      // fprintf(stderr,"\n");
      // }
      // debug end
      IGENMAX = -999999999999999999.9;
      for (int i = 0; i < 3; i++) {
        if (A1[i][i] > IGENMAX) {
          IGENMAX = A1[i][i];
          Y = i;
        }
      }
      if (bun == 2) {
        PT[MEN][0].MAXDISTANCE = IGENMAX * (double)(PT[MEN][0].INDEXNUM);
      }
      if (bun == 3) {
        PT[MEN][0].MAXDISTANCE = IGENMAX;
      }
      if (bun == 4 || bun == 5) {

        for (int i = 0; i < hsize * vsize; i++) {
          if (INDEX[i] == PT[MEN][0].INDEXNO) {
            RR[i] = YIN[i] - U_R;
            GG[i] = UIN[i] - U_G;
            BB[i] = VIN[i] - U_B;
          }
        }
        V[0] = X1[0][Y];
        V[1] = X1[1][Y];
        V[2] = X1[2][Y];
        for (int i = 0; i < hsize * vsize; i++) {
          kaiten(V, RR[i], GG[i], BB[i], RRR + i, GGG + i, BBB + i);
        }
        // for(i=0,l=0;i<hsize*vsize;i++){
        // if(index3[i]!=-1){
        // iRRRR[l] = (int)(RRR[i]+0.5);
        // l++;
        //}
        //}

        // int GASOSUU = hsize*vsize-k;
        // if(omh == 3 || omh == 4){
        // double *RRRRRR,*GGGGGG,*BBBBBB;
        //            RRRRRR = (double*)malloc(sizeof(double)*hsize*vsize);
        // GGGGGG = (double*)malloc(sizeof(double)*hsize*vsize);
        // BBBBBB = (double*)malloc(sizeof(double)*hsize*vsize);
        int l2 = 0;
        for (int i = 0; i < hsize * vsize; i++) {
          if (INDEX[i] == PT[MEN][0].INDEXNO) {
            RRRRRR[l2] = RRR[i];
            GGGGGG[l2] = GGG[i];
            BBBBBB[l2] = BBB[i];
            l2++;
          }
        }
        MAXX = -99.9e64;
        MINN = 99.9e64;
        for (int i = 0; i < l2; i++) {
          if (MAXX < RRRRRR[i]) {
            MAXX = RRRRRR[i];
          }
          if (MINN > RRRRRR[i]) {
            MINN = RRRRRR[i];
          }
        }
        if (bun == 4) {
          PT[MEN][0].MAXDISTANCE = MAXX - MINN;
        }
        if (bun == 5) {
          PT[MEN][0].MAXDISTANCE = (MAXX - MINN) * (double)PT[MEN][0].INDEXNUM;
        }

      } // bun == 4 if end

    } // bun == 1 if end
  } else {
    PT[MEN][0].MAXDISTANCE = 0.0;
  }

  // 1側のmaxdistanceを求める
  auto U_RGB2 =
      SumEquals(hsize * vsize, &INDEX[0], PT[MEN][1].INDEXNO, YIN, UIN, VIN);

  if (PT[MEN][1].INDEXNUM != 0) {
    std::get<0>(U_RGB2) /= (double)PT[MEN][1].INDEXNUM;
    std::get<1>(U_RGB2) /= (double)PT[MEN][1].INDEXNUM;
    std::get<2>(U_RGB2) /= (double)PT[MEN][1].INDEXNUM;
    double U_R2 = std::get<0>(U_RGB2);
    double U_G2 = std::get<1>(U_RGB2);
    double U_B2 = std::get<2>(U_RGB2);
    // U_R = (U_R>>RDIV);
    // U_G = (U_G>>GDIV);
    // U_B = (U_B>>BDIV);
    if (bun == 0) {
      MAXD = 0.0;
      for (int i = 0; i < hsize * vsize; i++) {
        if (INDEX[i] == PT[MEN][1].INDEXNO) {
          TEMP = ((YIN[i] - U_R2) * (YIN[i] - U_R2) +
                  (UIN[i] - U_G2) * (UIN[i] - U_G2) +
                  (VIN[i] - U_B2) * (VIN[i] - U_B2));
          if (MAXD < TEMP) {
            MAXD = TEMP;
          }
        }
      }
      PT[MEN][1].MAXDISTANCE = MAXD;
    } else if (bun == 1) {
      if (vvv != 4) {
        if (div == 0) {
          TEMP = 0.0;
          for (int i = 0; i < hsize * vsize; i++) {
            if (INDEX[i] == PT[MEN][1].INDEXNO) {
              TEMP += ((YIN[i] - U_R2) * (YIN[i] - U_R2) +
                       (UIN[i] - U_G2) * (UIN[i] - U_G2) +
                       (VIN[i] - U_B2) * (VIN[i] - U_B2));
              // if( MAXD < TEMP){
              // MAXD = TEMP;
              // }
            }
          }
        }
        if (div == 1) {
          TEMP = 0.0;
          for (int i = 0; i < hsize * vsize; i++) {
            if (INDEX[i] == PT[MEN][1].INDEXNO) {
              TEMP += pow(((YIN[i] - U_R2) * (YIN[i] - U_R2) +
                           (UIN[i] - U_G2) * (UIN[i] - U_G2) +
                           (VIN[i] - U_B2) * (VIN[i] - U_B2)),
                          pw / 2.0);
              // if( MAXD < TEMP){
              // MAXD = TEMP;
              // }
            }
          }
        }
      }

      if (vvv == 4) {
        if (div == 0) {
          TEMP = 0.0;
          for (int i = 0; i < hsize * vsize; i++) {
            if (INDEX[i] == PT[MEN][1].INDEXNO) {
              TEMP += ((YIN[i] - U_R2) * (YIN[i] - U_R2) +
                       (UIN[i] - U_G2) * (UIN[i] - U_G2) +
                       (VIN[i] - U_B2) * (VIN[i] - U_B2)) *
                      (255.0 - EDGERASISAY[i]) / 255.0;
              // if( MAXD < TEMP){
              // MAXD = TEMP;
              // }
            }
          }
        }
        if (div == 1) {
          TEMP = 0.0;
          for (int i = 0; i < hsize * vsize; i++) {
            if (INDEX[i] == PT[MEN][1].INDEXNO) {
              EDGERASISAYD = EDGERASISAY[i];
              EDGERASISAYD *= nbai;
              if (EDGERASISAYD > 255.0) {
                EDGERASISAYD = 255.0;
              }

              TEMP += pow(((YIN[i] - U_R2) * (YIN[i] - U_R2) +
                           (UIN[i] - U_G2) * (UIN[i] - U_G2) +
                           (VIN[i] - U_B2) * (VIN[i] - U_B2)),
                          pw / 2.0) *
                      pow((255.0 - EDGERASISAYD) / 255.0, pw2);
              // if( MAXD < TEMP){
              // MAXD = TEMP;
              // }
            }
          }
        }
      }

      PT[MEN][1].MAXDISTANCE = TEMP;
    } else if (bun == 2 || bun == 3 || bun == 4 || bun == 5) {
      JacobiSetupEquals(hsize * vsize, &INDEX[0], PT[MEN][1].INDEXNO, &YIN[0],
                        &UIN[0], &VIN[0], U_RGB2, A);
      ind = Jacobi(3, ct, eps, A, A1, A2, X1, X2);

      if (ind > 0) {
        fprintf(stderr, "Jacobi syuusoku simasen!!\n");
        exit(0);
      }
      // if((A1[0][0] == A1[1][1]) || (A1[1][1] == A1[2][2]) || (A1[0][0] ==
      // A1[2][2]) ){ fprintf(stderr,"Igen Value duplication!!\n"); exit(0);
      // }
      // debug start
      // for(i=0;i<3;i++){
      // for(j=0;j<3;j++){
      // fprintf(stderr,"%f ",A1[i][j]);
      // }
      // fprintf(stderr,"\n");
      // }
      // debug end
      IGENMAX = -999999999999999999.9;
      for (int i = 0; i < 3; i++) {
        if (A1[i][i] > IGENMAX) {
          IGENMAX = A1[i][i];
          Y = i;
        }
      }
      if (bun == 2) {
        PT[MEN][1].MAXDISTANCE = IGENMAX * (double)(PT[MEN][0].INDEXNUM);
      }
      if (bun == 3) {
        PT[MEN][1].MAXDISTANCE = IGENMAX;
      }
      // while(1);
      if (bun == 4 || bun == 5) {

        for (int i = 0; i < hsize * vsize; i++) {
          if (INDEX[i] == PT[MEN][1].INDEXNO) {
            RR[i] = YIN[i] - U_R2;
            GG[i] = UIN[i] - U_G2;
            BB[i] = VIN[i] - U_B2;
          }
        }
        V[0] = X1[0][Y];
        V[1] = X1[1][Y];
        V[2] = X1[2][Y];
        for (int i = 0; i < hsize * vsize; i++) {
          kaiten(V, RR[i], GG[i], BB[i], RRR + i, GGG + i, BBB + i);
        }
        // for(i=0,l=0;i<hsize*vsize;i++){
        // if(index3[i]!=-1){
        // iRRRR[l] = (int)(RRR[i]+0.5);
        // l++;
        //}
        //}

        // int GASOSUU = hsize*vsize-k;
        // if(omh == 3 || omh == 4){
        // double *RRRRRR,*GGGGGG,*BBBBBB;
        //            RRRRRR = (double*)malloc(sizeof(double)*hsize*vsize);
        // GGGGGG = (double*)malloc(sizeof(double)*hsize*vsize);
        // BBBBBB = (double*)malloc(sizeof(double)*hsize*vsize);
        int l2 = 0;
        for (int i = 0; i < hsize * vsize; i++) {
          if (INDEX[i] == PT[MEN][1].INDEXNO) {
            RRRRRR[l2] = RRR[i];
            GGGGGG[l2] = GGG[i];
            BBBBBB[l2] = BBB[i];
            l2++;
          }
        }
        MAXX = -99.9e64;
        MINN = 99.9e64;
        for (int i = 0; i < l2; i++) {
          if (MAXX < RRRRRR[i]) {
            MAXX = RRRRRR[i];
          }
          if (MINN > RRRRRR[i]) {
            MINN = RRRRRR[i];
          }
        }
        if (bun == 4) {
          PT[MEN][1].MAXDISTANCE = MAXX - MINN;
        }
        if (bun == 5) {
          PT[MEN][1].MAXDISTANCE = (MAXX - MINN) * (double)PT[MEN][1].INDEXNUM;
        }
      }

    } // bun == 1 if end
  } else {
    PT[MEN][1].MAXDISTANCE = 0.0;
  }

  // max distance div end

  // ohtsu tuika end
  // debug start
  // fprintf(stderr,"%d %d\n",k,l);
  // while(1);
  // debug end
  // 2回目の分轄
  //分轄の中の画素数が多いものを分轄
  // if(PT[MEN][0].INDEXNUM >= PT[MEN][1].INDEXNUM){

  int j2 = 0;
  if ((bun == 1) || (bun == 2) || (bun == 3) || (bun == 4) || (bun == 5)) {
    if (PT[MEN][0].MAXDISTANCE /**PT[MEN][0].INDEXNUM*/ >=
        PT[MEN][1]
            .MAXDISTANCE /**PT[MEN][1].INDEXNUM*/) { // maxdistance div tuika
      j2 = 0;                                        // INDEXNO
    } else {
      j2 = 1;
    }
  } else if (bun == 0) {
    if (PT[MEN][0].MAXDISTANCE * PT[MEN][0].INDEXNUM >=
        PT[MEN][1].MAXDISTANCE * PT[MEN][1].INDEXNUM) { // maxdistance div tuika
      j2 = 0;                                           // INDEXNO
    } else {
      j2 = 1;
    }
  }

  auto U_RGB3 =
      SumEquals(hsize * vsize, &INDEX[0], PT[MEN][j2].INDEXNO, YIN, UIN, VIN);
  std::get<0>(U_RGB3) /= (double)PT[MEN][j2].INDEXNUM;
  std::get<1>(U_RGB3) /= (double)PT[MEN][j2].INDEXNUM;
  std::get<2>(U_RGB3) /= (double)PT[MEN][j2].INDEXNUM;
  double U_R3 = std::get<0>(U_RGB3);
  double U_G3 = std::get<1>(U_RGB3);
  double U_B3 = std::get<2>(U_RGB3);
  // U_R = (U_R>>RDIV);
  // U_G = (U_G>>GDIV);
  // U_B = (U_B>>BDIV);
  // debug start
  // fprintf(stderr,"2kaimeheikin %d %d %d\n",U_R,U_G,U_B);
  // debug end
  JacobiSetupEquals(hsize * vsize, &INDEX[0], PT[MEN][j2].INDEXNO, &YIN[0],
                    &UIN[0], &VIN[0], U_RGB3, A);
  ind = Jacobi(3, ct, eps, A, A1, A2, X1, X2);

  if (ind > 0) {
    fprintf(stderr, "Jacobi syuusoku simasen!!\n");
    exit(0);
  }
  // if((A1[0][0] == A1[1][1]) || (A1[1][1] == A1[2][2]) || (A1[0][0] ==
  // A1[2][2]) ){ fprintf(stderr,"Igen Value duplication!!\n"); exit(0);
  // }
  // debug start
  // for(i=0;i<3;i++){
  // for(m=0;m<3;m++){
  // fprintf(stderr,"%f ",A1[i][m]);
  // }
  // fprintf(stderr,"\n");
  // }
  // debug end
  IGENMAX = -999999999999999999.9;
  for (int i = 0; i < 3; i++) {
    if (A1[i][i] > IGENMAX) {
      IGENMAX = A1[i][i];
      Y = i;
    }
  }
  DIVIDENUM = 3;

  // ohtsu tuika start
  for (int i = 0; i < hsize * vsize; i++) {
    if (INDEX[i] == j2) {
      RR[i] = YIN[i] - U_R3;
      GG[i] = UIN[i] - U_G3;
      BB[i] = VIN[i] - U_B3;
    }
  }
  V[0] = X1[0][Y];
  V[1] = X1[1][Y];
  V[2] = X1[2][Y];
  for (int i = 0; i < hsize * vsize; i++) {
    if (INDEX[i] == j2) {
      kaiten(V, RR[i], GG[i], BB[i], RRR + i, GGG + i, BBB + i);
    }
  }
  GASOSUU = PT[MEN][j2].INDEXNUM;
  for (int i = 0; i < hsize * vsize; i++) {
    RRR33[i] = RRR[i];
  }
  for (int i = 0, k = 0; i < hsize * vsize; i++) {
    if (INDEX[i] == j2) {
      RR[k] = RRR[i];
      GG[k] = GGG[i];
      BB[k] = BBB[i];
      k++;
    }
  }
  for (int i = 0; i < hsize * vsize; i++) {
    iRRRR[i] = (int)(RR[i] + 0.5);
  }
  // ohtsu tuika end
  if (MEN == 1) {
    MEN = 0;
    NMEN = 1;
  } else {
    MEN = 1;
    NMEN = 0;
  }
  // ohtsu tuika start
  if (omh == 3 || omh == 4) {
    THRESH = ohtsu2(GASOSUU, &RR[0], &GG[0], &BB[0], omh);
  }
  if (omh == 0) {
    THRESH = ohtsu(GASOSUU, iRRRR);
  } else if (omh == 1) {
    THRESH = media(GASOSUU, iRRRR);
  } else if (omh == 2) {
    THRESH = 0;
  }
  DTHRESH = (double)(THRESH);
  modoshi(V, DTHRESH, 0.0, 0.0, XXX, YYY, ZZZ);

  int k3 = 0, l3 = 0;
  for (int i = 0; i < hsize * vsize; i++) {
    if (INDEX[i] == j2) {
      // ZZ = X1[0][Y]*((double)(((YIN[i]))) - ((double)(U_R)+(*XXX))) +
      // X1[1][Y]*((double)(((UIN[i]))) - ((double)(U_G)+(*YYY))) +
      // X1[2][Y]*((double)(((VIN[i]))) - ((double)(U_B)+(*ZZZ))) ;
      // if(ZZ>=0.0){
      if (RRR33[i] >= DTHRESH) {
        INDEX[i] = DIVIDENUM - 1;
        k3++;
      } else {
        ;
        l3++;
      }
    }
  }
  PT[MEN][0].INDEXNO = j2;
  PT[MEN][1].INDEXNO = DIVIDENUM - 1;
  PT[MEN][0].INDEXNUM = l3;
  PT[MEN][1].INDEXNUM = k3;

  // max distance div tuika

  // maxdistance tuika start
  // 0側のmaxdistanceを求める
  auto U_RGB4 = SumEquals(hsize * vsize, &INDEX[0], PT[MEN][0].INDEXNO, &YIN[0],
                          &UIN[0], &VIN[0]);
  if (PT[MEN][0].INDEXNUM != 0) {
    std::get<0>(U_RGB4) /= (double)PT[MEN][0].INDEXNUM;
    std::get<1>(U_RGB4) /= (double)PT[MEN][0].INDEXNUM;
    std::get<2>(U_RGB4) /= (double)PT[MEN][0].INDEXNUM;
    double U_R4 = std::get<0>(U_RGB4);
    double U_G4 = std::get<1>(U_RGB4);
    double U_B4 = std::get<2>(U_RGB4);
    // U_R = (U_R>>RDIV);
    // U_G = (U_G>>GDIV);
    // U_B = (U_B>>BDIV);
    if (bun == 0) {
      MAXD = 0.0;
      for (int i = 0; i < hsize * vsize; i++) {
        if (INDEX[i] == PT[MEN][0].INDEXNO) {
          TEMP = (YIN[i] - U_R4) * (YIN[i] - U_R4) +
                 (UIN[i] - U_G4) * (UIN[i] - U_G4) +
                 (VIN[i] - U_B4) * (VIN[i] - U_B4);
          if (MAXD < TEMP) {
            MAXD = TEMP;
          }
        }
      }
      PT[MEN][0].MAXDISTANCE = MAXD;
    } else if (bun == 1) {
      // int TEMP;
      if (vvv != 4) {
        if (div == 0) {
          TEMP = 0.0;
          for (int i = 0; i < hsize * vsize; i++) {
            if (INDEX[i] == PT[MEN][0].INDEXNO) {
              TEMP += (YIN[i] - U_R4) * (YIN[i] - U_R4) +
                      (UIN[i] - U_G4) * (UIN[i] - U_G4) +
                      (VIN[i] - U_B4) * (VIN[i] - U_B4);
              // if( MAXD < TEMP){
              // MAXD = TEMP;
              // }
            }
          }
        }
        if (div == 1) {
          TEMP = 0.0;
          for (int i = 0; i < hsize * vsize; i++) {
            if (INDEX[i] == PT[MEN][0].INDEXNO) {
              TEMP += pow(((YIN[i] - U_R4) * (YIN[i] - U_R4) +
                           (UIN[i] - U_G4) * (UIN[i] - U_G4) +
                           (VIN[i] - U_B4) * (VIN[i] - U_B4)),
                          pw / 2.0);
              // if( MAXD < TEMP){
              // MAXD = TEMP;
              // }
            }
          }
        }
      }

      if (vvv == 4) {
        if (div == 0) {
          TEMP = 0.0;
          for (int i = 0; i < hsize * vsize; i++) {
            if (INDEX[i] == PT[MEN][0].INDEXNO) {
              TEMP += ((YIN[i] - U_R4) * (YIN[i] - U_R4) +
                       (UIN[i] - U_G4) * (UIN[i] - U_G4) +
                       (VIN[i] - U_B4) * (VIN[i] - U_B4)) *
                      (255.0 - EDGERASISAY[i]) / 255.0;
              // if( MAXD < TEMP){
              // MAXD = TEMP;
              // }
            }
          }
        }
        if (div == 1) {
          TEMP = 0.0;
          for (int i = 0; i < hsize * vsize; i++) {
            if (INDEX[i] == PT[MEN][0].INDEXNO) {
              EDGERASISAYD = EDGERASISAY[i];
              EDGERASISAYD *= nbai;
              if (EDGERASISAYD > 255.0) {
                EDGERASISAYD = 255.0;
              }

              TEMP += pow(((YIN[i] - U_R4) * (YIN[i] - U_R4) +
                           (UIN[i] - U_G4) * (UIN[i] - U_G4) +
                           (VIN[i] - U_B4) * (VIN[i] - U_B4)),
                          pw / 2.0) *
                      pow((255.0 - EDGERASISAYD) / 255.0, pw2);
              // if( MAXD < TEMP){
              // MAXD = TEMP;
              // }
            }
          }
        }
      }

      PT[MEN][0].MAXDISTANCE = TEMP;
    } else if (bun == 2 || bun == 3 || bun == 4 || bun == 5) {
      JacobiSetupEquals(hsize * vsize, &INDEX[0], PT[MEN][0].INDEXNO, &YIN[0],
                        &UIN[0], &VIN[0], U_RGB4, A);
      ind = Jacobi(3, ct, eps, A, A1, A2, X1, X2);

      if (ind > 0) {
        fprintf(stderr, "Jacobi syuusoku simasen!!\n");
        exit(0);
      }
      // if((A1[0][0] == A1[1][1]) || (A1[1][1] == A1[2][2]) || (A1[0][0] ==
      // A1[2][2]) ){ fprintf(stderr,"Igen Value duplication!!\n"); exit(0);
      // }
      // debug start
      // for(i=0;i<3;i++){
      // for(m=0;m<3;m++){
      // fprintf(stderr,"%f ",A1[i][m]);
      // }
      // fprintf(stderr,"\n");
      // }
      // debug end
      IGENMAX = -999999999999999999.9;
      for (int i = 0; i < 3; i++) {
        if (A1[i][i] > IGENMAX) {
          IGENMAX = A1[i][i];
          Y = i;
        }
      }
      if (bun == 2) {
        PT[MEN][0].MAXDISTANCE = IGENMAX * (double)(PT[MEN][0].INDEXNUM);
      }
      if (bun == 3) {
        PT[MEN][0].MAXDISTANCE = IGENMAX;
      }
      if (bun == 4 || bun == 5) {

        for (int i = 0; i < hsize * vsize; i++) {
          if (INDEX[i] == PT[MEN][0].INDEXNO) {
            RR[i] = YIN[i] - U_R4;
            GG[i] = UIN[i] - U_G4;
            BB[i] = VIN[i] - U_B4;
          }
        }
        V[0] = X1[0][Y];
        V[1] = X1[1][Y];
        V[2] = X1[2][Y];
        for (int i = 0; i < hsize * vsize; i++) {
          kaiten(V, RR[i], GG[i], BB[i], RRR + i, GGG + i, BBB + i);
        }
        // for(i=0,l=0;i<hsize*vsize;i++){
        // if(index3[i]!=-1){
        // iRRRR[l] = (int)(RRR[i]+0.5);
        // l++;
        //}
        //}

        // int GASOSUU = hsize*vsize-k;
        // if(omh == 3 || omh == 4){
        // double *RRRRRR,*GGGGGG,*BBBBBB;
        //            RRRRRR = (double*)malloc(sizeof(double)*hsize*vsize);
        // GGGGGG = (double*)malloc(sizeof(double)*hsize*vsize);
        // BBBBBB = (double*)malloc(sizeof(double)*hsize*vsize);
        int l2 = 0;
        for (int i = 0; i < hsize * vsize; i++) {
          if (INDEX[i] == PT[MEN][0].INDEXNO) {
            RRRRRR[l2] = RRR[i];
            GGGGGG[l2] = GGG[i];
            BBBBBB[l2] = BBB[i];
            l2++;
          }
        }
        MAXX = -99.9e64;
        MINN = 99.9e64;
        for (int i = 0; i < l2; i++) {
          if (MAXX < RRRRRR[i]) {
            MAXX = RRRRRR[i];
          }
          if (MINN > RRRRRR[i]) {
            MINN = RRRRRR[i];
          }
        }
        if (bun == 4) {
          PT[MEN][0].MAXDISTANCE = MAXX - MINN;
        }
        if (bun == 5) {
          PT[MEN][0].MAXDISTANCE = (MAXX - MINN) * (double)PT[MEN][0].INDEXNUM;
        }
      }

    } // bun == 1 if end
  } else {
    PT[MEN][0].MAXDISTANCE = 0.0;
  }

  // 1側のmaxdistanceを求める
  auto U_RGB5 = SumEquals(hsize * vsize, &INDEX[0], PT[MEN][1].INDEXNO, &YIN[0],
                          &UIN[0], &VIN[0]);

  if (PT[MEN][1].INDEXNUM != 0) {
    std::get<0>(U_RGB5) /= (double)PT[MEN][1].INDEXNUM;
    std::get<1>(U_RGB5) /= (double)PT[MEN][1].INDEXNUM;
    std::get<2>(U_RGB5) /= (double)PT[MEN][1].INDEXNUM;
    double U_R5 = std::get<0>(U_RGB5);
    double U_G5 = std::get<1>(U_RGB5);
    double U_B5 = std::get<2>(U_RGB5);
    // U_R = (U_R>>RDIV);
    // U_G = (U_G>>GDIV);
    // U_B = (U_B>>BDIV);
    if (bun == 0) {
      MAXD = 0.0;
      for (int i = 0; i < hsize * vsize; i++) {
        if (INDEX[i] == PT[MEN][1].INDEXNO) {
          TEMP = (YIN[i] - U_R5) * (YIN[i] - U_R5) +
                 (UIN[i] - U_G5) * (UIN[i] - U_G5) +
                 (VIN[i] - U_B5) * (VIN[i] - U_B5);
          if (MAXD < TEMP) {
            MAXD = TEMP;
          }
        }
      }
      PT[MEN][1].MAXDISTANCE = MAXD;
    } else if (bun == 1) {
      if (vvv != 4) {
        if (div == 0) {
          TEMP = 0.0;
          for (int i = 0; i < hsize * vsize; i++) {
            if (INDEX[i] == PT[MEN][1].INDEXNO) {
              TEMP += (YIN[i] - U_R5) * (YIN[i] - U_R5) +
                      (UIN[i] - U_G5) * (UIN[i] - U_G5) +
                      (VIN[i] - U_B5) * (VIN[i] - U_B5);
              // if( MAXD < TEMP){
              // MAXD = TEMP;
              // }
            }
          }
        }
        if (div == 1) {
          TEMP = 0.0;
          for (int i = 0; i < hsize * vsize; i++) {
            if (INDEX[i] == PT[MEN][1].INDEXNO) {
              TEMP += pow(((YIN[i] - U_R5) * (YIN[i] - U_R5) +
                           (UIN[i] - U_G5) * (UIN[i] - U_G5) +
                           (VIN[i] - U_B5) * (VIN[i] - U_B5)),
                          pw / 2.0);
              // if( MAXD < TEMP){
              // MAXD = TEMP;
              // }
            }
          }
        }
      }

      if (vvv == 4) {
        if (div == 0) {
          TEMP = 0.0;
          for (int i = 0; i < hsize * vsize; i++) {
            if (INDEX[i] == PT[MEN][1].INDEXNO) {
              TEMP += ((YIN[i] - U_R5) * (YIN[i] - U_R5) +
                       (UIN[i] - U_G5) * (UIN[i] - U_G5) +
                       (VIN[i] - U_B5) * (VIN[i] - U_B5)) *
                      (255.0 - EDGERASISAY[i]) / 255.0;
              // if( MAXD < TEMP){
              // MAXD = TEMP;
              // }
            }
          }
        }
        if (div == 1) {
          TEMP = 0.0;
          for (int i = 0; i < hsize * vsize; i++) {
            if (INDEX[i] == PT[MEN][1].INDEXNO) {
              EDGERASISAYD = EDGERASISAY[i];
              EDGERASISAYD *= nbai;
              if (EDGERASISAYD > 255.0) {
                EDGERASISAYD = 255.0;
              }

              TEMP += pow(((YIN[i] - U_R5) * (YIN[i] - U_R5) +
                           (UIN[i] - U_G5) * (UIN[i] - U_G5) +
                           (VIN[i] - U_B5) * (VIN[i] - U_B5)),
                          pw / 2.0) *
                      pow((255.0 - EDGERASISAYD) / 255.0, pw2);
              // if( MAXD < TEMP){
              // MAXD = TEMP;
              // }
            }
          }
        }
      }

      PT[MEN][1].MAXDISTANCE = TEMP;
    } else if (bun == 2 || bun == 3 || bun == 4 || bun == 5) {
      JacobiSetupEquals(hsize * vsize, &INDEX[0], PT[MEN][1].INDEXNO, &YIN[0],
                        &UIN[0], &VIN[0], U_RGB5, A);
      ind = Jacobi(3, ct, eps, A, A1, A2, X1, X2);

      if (ind > 0) {
        fprintf(stderr, "Jacobi syuusoku simasen!!\n");
        exit(0);
      }
      // if((A1[0][0] == A1[1][1]) || (A1[1][1] == A1[2][2]) || (A1[0][0] ==
      // A1[2][2]) ){ fprintf(stderr,"Igen Value duplication!!\n"); exit(0);
      // }
      // debug start
      // for(i=0;i<3;i++){
      // for(m=0;m<3;m++){
      // fprintf(stderr,"%f ",A1[i][m]);
      // }
      // fprintf(stderr,"\n");
      // }
      // debug end

      IGENMAX = -999999999999999999.9;
      for (int i = 0; i < 3; i++) {
        if (A1[i][i] > IGENMAX) {
          IGENMAX = A1[i][i];
          Y = i;
        }
      }
      if (bun == 2) {
        PT[MEN][1].MAXDISTANCE = IGENMAX * (double)(PT[MEN][1].INDEXNUM);
      }
      if (bun == 3) {
        PT[MEN][1].MAXDISTANCE = IGENMAX;
      }
      if (bun == 4 || bun == 5) {
        for (int i = 0; i < hsize * vsize; i++) {
          if (INDEX[i] == PT[MEN][1].INDEXNO) {
            RR[i] = YIN[i] - U_R5;
            GG[i] = UIN[i] - U_G5;
            BB[i] = VIN[i] - U_B5;
          }
        }
        V[0] = X1[0][Y];
        V[1] = X1[1][Y];
        V[2] = X1[2][Y];
        for (int i = 0; i < hsize * vsize; i++) {
          kaiten(V, RR[i], GG[i], BB[i], RRR + i, GGG + i, BBB + i);
        }
        // for(i=0,l=0;i<hsize*vsize;i++){
        // if(index3[i]!=-1){
        // iRRRR[l] = (int)(RRR[i]+0.5);
        // l++;
        //}
        //}

        // int GASOSUU = hsize*vsize-k;
        // if(omh == 3 || omh == 4){
        // double *RRRRRR,*GGGGGG,*BBBBBB;
        //            RRRRRR = (double*)malloc(sizeof(double)*hsize*vsize);
        // GGGGGG = (double*)malloc(sizeof(double)*hsize*vsize);
        // BBBBBB = (double*)malloc(sizeof(double)*hsize*vsize);
        int l2 = 0;
        for (int i = 0; i < hsize * vsize; i++) {
          if (INDEX[i] == PT[MEN][1].INDEXNO) {
            RRRRRR[l2] = RRR[i];
            GGGGGG[l2] = GGG[i];
            BBBBBB[l2] = BBB[i];
            l2++;
          }
        }
        MAXX = -99.9e64;
        MINN = 99.9e64;
        for (int i = 0; i < l2; i++) {
          if (MAXX < RRRRRR[i]) {
            MAXX = RRRRRR[i];
          }
          if (MINN > RRRRRR[i]) {
            MINN = RRRRRR[i];
          }
        }
        if (bun == 4) {
          PT[MEN][1].MAXDISTANCE = MAXX - MINN;
        }
        if (bun == 5) {
          PT[MEN][1].MAXDISTANCE = (MAXX - MINN) * (double)PT[MEN][1].INDEXNUM;
        }
      }

    } // bun == 1 if end
  } else {
    PT[MEN][1].MAXDISTANCE = 0.0;
  }
  // max distance div end

  if (j2 == 0) {
    for (int i = 2; i <= DIVIDENUM - 1; i++) {
      PT[MEN][i].INDEXNO = PT[NMEN][i - 1].INDEXNO;
      PT[MEN][i].INDEXNUM = PT[NMEN][i - 1].INDEXNUM;
      PT[MEN][i].MAXDISTANCE =
          PT[NMEN][i - 1].MAXDISTANCE; // max distance div tuika
    }
  }
  if ((j2 >= 1) && (j2 <= (DIVIDENUM - 3))) {
    for (int i = 2; i <= j2 + 1; i++) {
      //// i i-2
      PT[MEN][i].INDEXNO = PT[NMEN][i - 2].INDEXNO;
      PT[MEN][i].INDEXNUM = PT[NMEN][i - 2].INDEXNUM;
      PT[MEN][i].MAXDISTANCE =
          PT[NMEN][i - 2].MAXDISTANCE; // max distance div tuika
    }

    for (int i = j2 + 2; i <= DIVIDENUM - 1; i++) {
      //// i i-1
      PT[MEN][i].INDEXNO = PT[NMEN][i - 1].INDEXNO;
      PT[MEN][i].INDEXNUM = PT[NMEN][i - 1].INDEXNUM;
      PT[MEN][i].MAXDISTANCE =
          PT[NMEN][i - 1].MAXDISTANCE; // max distance div tuika
    }
  }
  if (j2 == (DIVIDENUM - 2)) {
    for (int i = 2; i <= j2 + 1; i++) {
      //// i i-2
      PT[MEN][i].INDEXNO = PT[NMEN][i - 2].INDEXNO;
      PT[MEN][i].INDEXNUM = PT[NMEN][i - 2].INDEXNUM;
      PT[MEN][i].MAXDISTANCE =
          PT[NMEN][i - 2].MAXDISTANCE; // max distance div tuika
    }
  }

  // 9/9 kokokara
  while (DIVIDENUM < 256) {
    StopWatch sw_before;
    //最大画素数のブロックをさがし、そのINDEXNOを調べる。
    MAXINDEXNUM = -1.0;

    if ((bun == 1) || (bun == 2) || (bun == 3) || (bun == 4) || (bun == 5)) {
      for (int i = 0; i < DIVIDENUM; i++) {
        // if(PT[MEN][i].INDEXNUM > MAXINDEXNUM)
        // MAXINDEXNUM = PT[MEN][i].INDEXNUM;
        if (PT[MEN][i].MAXDISTANCE /**PT[MEN][i].INDEXNUM*/ > MAXINDEXNUM) {
          MAXINDEXNUM =
              PT[MEN][i]
                  .MAXDISTANCE /**PT[MEN][i].INDEXNUM*/; // max distance tuika
          NUM = i;
        }
      }
    } else if (bun == 0) {
      for (int i = 0; i < DIVIDENUM; i++) {
        // if(PT[MEN][i].INDEXNUM > MAXINDEXNUM)
        // MAXINDEXNUM = PT[MEN][i].INDEXNUM;
        if (PT[MEN][i].MAXDISTANCE * PT[MEN][i].INDEXNUM > MAXINDEXNUM) {
          MAXINDEXNUM = PT[MEN][i].MAXDISTANCE *
                        PT[MEN][i].INDEXNUM; // max distance tuika
          NUM = i;
        }
      }
    }

    // debug start
    if (DIVIDENUM == 41) {
      fprintf(stderr, "NUM= %d\n", NUM);
      fprintf(stderr, "PT[MEN][%d].INDEXNO=%d\n", NUM, PT[MEN][NUM].INDEXNO);
    }
    // debug end
    //次に分轄するブロックNO----- NUM
    //次に分轄するINDEXNO-------- PT[MEN][NUM].INDEXNO
    auto U_RGB6 = SumEquals(hsize * vsize, &INDEX[0], PT[MEN][NUM].INDEXNO,
                            &YIN[0], &UIN[0], &VIN[0]);
    std::get<0>(U_RGB6) /= (double)PT[MEN][NUM].INDEXNUM;
    std::get<1>(U_RGB6) /= (double)PT[MEN][NUM].INDEXNUM;
    std::get<2>(U_RGB6) /= (double)PT[MEN][NUM].INDEXNUM;
    double U_R6 = std::get<0>(U_RGB6);
    double U_G6 = std::get<1>(U_RGB6);
    double U_B6 = std::get<2>(U_RGB6);
    // U_R = (U_R>>RDIV);
    // U_G = (U_G>>GDIV);
    // U_B = (U_B>>BDIV);
    // debug start
    // if(DIVIDENUM == 41){
    // fprintf(stderr,"AVE= %d %d %d\n",U_R,U_G,U_B);
    // fprintf(stderr,"NUM= %d\n",NUM);
    // }
    // debug end
    JacobiSetupEquals(hsize * vsize, &INDEX[0], PT[MEN][NUM].INDEXNO, &YIN[0],
                      &UIN[0], &VIN[0], U_RGB6, A);
    // debug start
    if (DIVIDENUM == 41) {
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          fprintf(stderr, "%f ", A[i][j]);
        }
        fprintf(stderr, "\n");
      }
    }
    // debug end
    ind = Jacobi(3, ct, eps, A, A1, A2, X1, X2);

    if (ind > 0) {
      fprintf(stderr, "Jacobi syuusoku simasen!!\n");
      exit(0);
    }
    // if((A1[0][0] == A1[1][1]) || (A1[1][1] == A1[2][2]) || (A1[0][0] ==
    // A1[2][2]) ){ fprintf(stderr,"Igen Value duplication!!\n"); exit(0);
    // }
    // debug start
    // if(DIVIDENUM == 4){
    // for(i=0;i<3;i++){
    // for(m=0;m<3;m++){
    // fprintf(stderr,"%f ",A1[i][m]);
    // }
    // fprintf(stderr,"\n");
    // }
    //}
    // debug end
    IGENMAX = -999999999999999999.9;
    for (int i = 0; i < 3; i++) {
      if (A1[i][i] > IGENMAX) {
        IGENMAX = A1[i][i];
        Y = i;
      }
    }
    DIVIDENUM++;

    // ohtsu tuika start
    for (int i = 0; i < hsize * vsize; i++) {
      if (INDEX[i] == PT[MEN][NUM].INDEXNO) {
        RR[i] = YIN[i] - U_R6;
        GG[i] = UIN[i] - U_G6;
        BB[i] = VIN[i] - U_B6;
      }
    }
    V[0] = X1[0][Y];
    V[1] = X1[1][Y];
    V[2] = X1[2][Y];
    for (int i = 0; i < hsize * vsize; i++) {
      if (INDEX[i] == PT[MEN][NUM].INDEXNO) {
        kaiten(V, RR[i], GG[i], BB[i], RRR + i, GGG + i, BBB + i);
      }
    }
    GASOSUU = PT[MEN][NUM].INDEXNUM;
    for (int i = 0; i < hsize * vsize; i++) {
      RRR33[i] = RRR[i];
    }
    {
      int k = 0;
      for (int i = 0; i < hsize * vsize; i++) {
        if (INDEX[i] == PT[MEN][NUM].INDEXNO) {
          RR[k] = RRR[i];
          GG[k] = GGG[i];
          BB[k] = BBB[i];
          k++;
        }
      }
    }
    for (int i = 0; i < hsize * vsize; i++) {
      iRRRR[i] = (int)(RR[i] + 0.5);
    }
    // ohtsu tuika end
    // if(MEN == 1){
    // MEN = 0;
    // NMEN = 1;
    //// } else {
    // MEN = 1;
    // NMEN = 0;
    // }
    // ohtsu tuika start
    printf("before %8.4f ms\n", sw_before.lap());
    if (omh == 3 || omh == 4) {
      printf("hoge4\n");
      StopWatch sw_ohtsu2;
      THRESH = ohtsu2(GASOSUU, &RR[0], &GG[0], &BB[0], omh);
      printf("ohtsu2 %8.4f ms\n", sw_ohtsu2.lap());
    }
    if (omh == 0) {
      THRESH = ohtsu(GASOSUU, iRRRR);
    } else if (omh == 1) {
      THRESH = media(GASOSUU, iRRRR);
    } else if (omh == 2) {
      THRESH = 0;
    }
    StopWatch sw_after;
    DTHRESH = (double)(THRESH);
    modoshi(V, DTHRESH, 0.0, 0.0, XXX, YYY, ZZZ);

    int k2 = 0, l2 = 0;
    for (int i = 0; i < hsize * vsize; i++) {
      if (INDEX[i] == PT[MEN][NUM].INDEXNO) {
        // ZZ = X1[0][Y]*((double)(((YIN[i]))) - ((double)(U_R)+(*XXX))) +
        // X1[1][Y]*((double)(((UIN[i]))) - ((double)(U_G)+(*YYY))) +
        // X1[2][Y]*((double)(((VIN[i]))) - ((double)(U_B)+(*ZZZ))) ;
        // if(ZZ>=0.0){
        if (RRR33[i] >= DTHRESH) {
          INDEX[i] = DIVIDENUM - 1;
          k2++;
        } else {
          ;
          l2++;
        }
      }
    }
    if (MEN == 1) {
      MEN = 0;
      NMEN = 1;
    } else {
      MEN = 1;
      NMEN = 0;
    }
    PT[MEN][0].INDEXNO = PT[NMEN][NUM].INDEXNO;
    PT[MEN][1].INDEXNO = DIVIDENUM - 1;
    PT[MEN][0].INDEXNUM = l2;
    PT[MEN][1].INDEXNUM = k2;
    // debug start
    if (DIVIDENUM == 5 || DIVIDENUM == 6) {
      fprintf(stderr, "INDEXNO%d %d\n", PT[MEN][0].INDEXNO, l2);
      fprintf(stderr, "INDEXNO%d %d\n", PT[MEN][1].INDEXNO, k2);
    }
    // while(1);
    // debug end

    // max distance divide tuika
    // maxdistance tuika start
    // 0側のmaxdistanceを求める
    auto U_RGB7 = SumEquals(hsize * vsize, &INDEX[0], PT[MEN][0].INDEXNO,
                            &YIN[0], &UIN[0], &VIN[0]);

    if (PT[MEN][0].INDEXNUM != 0) {
      std::get<0>(U_RGB7) /= PT[MEN][0].INDEXNUM;
      std::get<1>(U_RGB7) /= PT[MEN][0].INDEXNUM;
      std::get<2>(U_RGB7) /= PT[MEN][0].INDEXNUM;
      double U_R7 = std::get<0>(U_RGB7);
      double U_G7 = std::get<1>(U_RGB7);
      double U_B7 = std::get<2>(U_RGB7);
      // U_R = (U_R>>RDIV);
      // U_G = (U_G>>GDIV);
      // U_B = (U_B>>BDIV);
      if (bun == 0) {
        MAXD = 0.0;
        for (int i = 0; i < hsize * vsize; i++) {
          if (INDEX[i] == PT[MEN][0].INDEXNO) {
            TEMP = (YIN[i] - U_R7) * (YIN[i] - U_R7) +
                   (UIN[i] - U_G7) * (UIN[i] - U_G7) +
                   (VIN[i] - U_B7) * (VIN[i] - U_B7);
            if (MAXD < TEMP) {
              MAXD = TEMP;
            }
          }
        }
        PT[MEN][0].MAXDISTANCE = MAXD;
      } else if (bun == 1) {
        // int TEMP;
        if (vvv != 4) {
          if (div == 0) {
            TEMP = 0.0;
            for (int i = 0; i < hsize * vsize; i++) {
              if (INDEX[i] == PT[MEN][0].INDEXNO) {
                TEMP += (YIN[i] - U_R7) * (YIN[i] - U_R7) +
                        (UIN[i] - U_G7) * (UIN[i] - U_G7) +
                        (VIN[i] - U_B7) * (VIN[i] - U_B7);
                // if( MAXD < TEMP){
                // MAXD = TEMP;
                // }
              }
            }
          }
          if (div == 1) {
            TEMP = 0.0;
            for (int i = 0; i < hsize * vsize; i++) {
              if (INDEX[i] == PT[MEN][0].INDEXNO) {
                TEMP += pow(((YIN[i] - U_R7) * (YIN[i] - U_R7) +
                             (UIN[i] - U_G7) * (UIN[i] - U_G7) +
                             (VIN[i] - U_B7) * (VIN[i] - U_B7)),
                            pw / 2.0);
                // if( MAXD < TEMP){
                // MAXD = TEMP;
                // }
              }
            }
          }
        }

        if (vvv == 4) {
          if (div == 0) {
            TEMP = 0.0;
            for (int i = 0; i < hsize * vsize; i++) {
              if (INDEX[i] == PT[MEN][0].INDEXNO) {
                TEMP += ((YIN[i] - U_R7) * (YIN[i] - U_R7) +
                         (UIN[i] - U_G7) * (UIN[i] - U_G7) +
                         (VIN[i] - U_B7) * (VIN[i] - U_B7)) *
                        (255.0 - EDGERASISAY[i]) / 255.0;
                // if( MAXD < TEMP){
                // MAXD = TEMP;
                // }
              }
            }
          }
          if (div == 1) {
            TEMP = 0.0;
            for (int i = 0; i < hsize * vsize; i++) {
              if (INDEX[i] == PT[MEN][0].INDEXNO) {
                EDGERASISAYD = EDGERASISAY[i];
                EDGERASISAYD *= nbai;
                if (EDGERASISAYD > 255.0) {
                  EDGERASISAYD = 255.0;
                }

                TEMP += pow(((YIN[i] - U_R7) * (YIN[i] - U_R7) +
                             (UIN[i] - U_G7) * (UIN[i] - U_G7) +
                             (VIN[i] - U_B7) * (VIN[i] - U_B7)),
                            pw / 2.0) *
                        pow((255.0 - EDGERASISAYD) / 255.0, pw2);
                // if( MAXD < TEMP){
                // MAXD = TEMP;
                // }
              }
            }
          }
        }

        PT[MEN][0].MAXDISTANCE = TEMP;
      } else if (bun == 2 || bun == 3 || bun == 4 || bun == 5) {
        JacobiSetupEquals(hsize * vsize, &INDEX[0], PT[MEN][0].INDEXNO, &YIN[0],
                          &UIN[0], &VIN[0], U_RGB7, A);
        ind = Jacobi(3, ct, eps, A, A1, A2, X1, X2);

        if (ind > 0) {
          fprintf(stderr, "Jacobi syuusoku simasen!!\n");
          exit(0);
        }
        // if((A1[0][0] == A1[1][1]) || (A1[1][1] == A1[2][2]) || (A1[0][0] ==
        // A1[2][2]) ){ fprintf(stderr,"Igen Value duplication!!\n"); exit(0);
        // }
        // debug start
        // if(DIVIDENUM == 4){
        // for(i=0;i<3;i++){
        // for(m=0;m<3;m++){
        // fprintf(stderr,"%f ",A1[i][m]);
        // }
        // fprintf(stderr,"\n");
        // }
        //}
        // debug end
        IGENMAX = -999999999999999999.9;
        for (int i = 0; i < 3; i++) {
          if (A1[i][i] > IGENMAX) {
            IGENMAX = A1[i][i];
            Y = i;
          }
        }
        if (bun == 2) {
          PT[MEN][0].MAXDISTANCE = IGENMAX * (double)PT[MEN][0].INDEXNUM;
        }
        if (bun == 3) {
          PT[MEN][0].MAXDISTANCE = IGENMAX;
        }
        if (bun == 4 || bun == 5) {
          for (int i = 0; i < hsize * vsize; i++) {
            if (INDEX[i] == PT[MEN][0].INDEXNO) {
              RR[i] = YIN[i] - U_R7;
              GG[i] = UIN[i] - U_G7;
              BB[i] = VIN[i] - U_B7;
            }
          }
          V[0] = X1[0][Y];
          V[1] = X1[1][Y];
          V[2] = X1[2][Y];
          for (int i = 0; i < hsize * vsize; i++) {
            kaiten(V, RR[i], GG[i], BB[i], RRR + i, GGG + i, BBB + i);
          }
          // for(i=0,l=0;i<hsize*vsize;i++){
          // if(index3[i]!=-1){
          // iRRRR[l] = (int)(RRR[i]+0.5);
          // l++;
          //}
          //}

          // int GASOSUU = hsize*vsize-k;
          // if(omh == 3 || omh == 4){
          // double *RRRRRR,*GGGGGG,*BBBBBB;
          //            RRRRRR = (double*)malloc(sizeof(double)*hsize*vsize);
          // GGGGGG = (double*)malloc(sizeof(double)*hsize*vsize);
          // BBBBBB = (double*)malloc(sizeof(double)*hsize*vsize);
          int l2 = 0;
          for (int i = 0; i < hsize * vsize; i++) {
            if (INDEX[i] == PT[MEN][0].INDEXNO) {
              RRRRRR[l2] = RRR[i];
              GGGGGG[l2] = GGG[i];
              BBBBBB[l2] = BBB[i];
              l2++;
            }
          }
          MAXX = -99.9e64;
          MINN = 99.9e64;
          for (int i = 0; i < l2; i++) {
            if (MAXX < RRRRRR[i]) {
              MAXX = RRRRRR[i];
            }
            if (MINN > RRRRRR[i]) {
              MINN = RRRRRR[i];
            }
          }
          if (bun == 4) {
            PT[MEN][0].MAXDISTANCE = MAXX - MINN;
          }
          if (bun == 5) {
            PT[MEN][0].MAXDISTANCE =
                (MAXX - MINN) * (double)PT[MEN][0].INDEXNUM;
          }
        }

      } // bun == 1 if end
      printf("after %8.4f ms\n", sw_after.lap());
    } else {
      PT[MEN][0].MAXDISTANCE = 0.0;
    }

    // 1側のmaxdistanceを求める
    auto U_RGB8 = SumEquals(hsize * vsize, &INDEX[0], PT[MEN][1].INDEXNO,
                            &YIN[0], &UIN[0], &VIN[0]);
    if (PT[MEN][1].INDEXNUM != 0) {
      std::get<0>(U_RGB8) /= PT[MEN][1].INDEXNUM;
      std::get<1>(U_RGB8) /= PT[MEN][1].INDEXNUM;
      std::get<2>(U_RGB8) /= PT[MEN][1].INDEXNUM;
      double U_R8 = std::get<0>(U_RGB8);
      double U_G8 = std::get<1>(U_RGB8);
      double U_B8 = std::get<2>(U_RGB8);
      // U_R = (U_R>>RDIV);
      // U_G = (U_G>>GDIV);
      // U_B = (U_B>>BDIV);
      if (bun == 0) {
        MAXD = 0.0;
        for (int i = 0; i < hsize * vsize; i++) {
          if (INDEX[i] == PT[MEN][1].INDEXNO) {
            TEMP = (YIN[i] - U_R8) * (YIN[i] - U_R8) +
                   (UIN[i] - U_G8) * (UIN[i] - U_G8) +
                   (VIN[i] - U_B8) * (VIN[i] - U_B8);
            if (MAXD < TEMP) {
              MAXD = TEMP;
            }
          }
        }
        PT[MEN][1].MAXDISTANCE = MAXD;
      } else if (bun == 1) {
        //    if(div != 4){   //2019.9.15  bug syuusei
        if (vvv != 4) {
          if (div == 0) {
            TEMP = 0.0;
            for (int i = 0; i < hsize * vsize; i++) {
              if (INDEX[i] == PT[MEN][1].INDEXNO) {
                TEMP += (YIN[i] - U_R8) * (YIN[i] - U_R8) +
                        (UIN[i] - U_G8) * (UIN[i] - U_G8) +
                        (VIN[i] - U_B8) * (VIN[i] - U_B8);
                // if( MAXD < TEMP){
                // MAXD = TEMP;
                // }
              }
            }
          }
          if (div == 1) {
            TEMP = 0.0;
            for (int i = 0; i < hsize * vsize; i++) {
              if (INDEX[i] == PT[MEN][1].INDEXNO) {
                TEMP += pow(((YIN[i] - U_R8) * (YIN[i] - U_R8) +
                             (UIN[i] - U_G8) * (UIN[i] - U_G8) +
                             (VIN[i] - U_B8) * (VIN[i] - U_B8)),
                            pw / 2.0);
                // if( MAXD < TEMP){
                // MAXD = TEMP;
                // }
              }
            }
          }
        }
        if (vvv == 4) {
          if (div == 0) {
            TEMP = 0.0;
            for (int i = 0; i < hsize * vsize; i++) {
              if (INDEX[i] == PT[MEN][1].INDEXNO) {
                TEMP += ((YIN[i] - U_R8) * (YIN[i] - U_R8) +
                         (UIN[i] - U_G8) * (UIN[i] - U_G8) +
                         (VIN[i] - U_B8) * (VIN[i] - U_B8)) *
                        (255.0 - EDGERASISAY[i]) / 255.0;
                // if( MAXD < TEMP){
                // MAXD = TEMP;
                // }
              }
            }
          }
          if (div == 1) {
            TEMP = 0.0;
            for (int i = 0; i < hsize * vsize; i++) {
              if (INDEX[i] == PT[MEN][1].INDEXNO) {
                EDGERASISAYD = EDGERASISAY[i];
                EDGERASISAYD *= nbai;
                if (EDGERASISAYD > 255.0) {
                  EDGERASISAYD = 255.0;
                }

                TEMP += pow(((YIN[i] - U_R8) * (YIN[i] - U_R8) +
                             (UIN[i] - U_G8) * (UIN[i] - U_G8) +
                             (VIN[i] - U_B8) * (VIN[i] - U_B8)),
                            pw / 2.0) *
                        pow((255.0 - EDGERASISAYD) / 255.0, pw2);
                // if( MAXD < TEMP){
                // MAXD = TEMP;
                // }
              }
            }
          }
        }

        PT[MEN][1].MAXDISTANCE = TEMP;
      } else if (bun == 2 || bun == 3 || bun == 4 || bun == 5) {
        JacobiSetupEquals(hsize * vsize, &INDEX[0], PT[MEN][1].INDEXNO, &YIN[0],
                          &UIN[0], &VIN[0], U_RGB8, A);
        ind = Jacobi(3, ct, eps, A, A1, A2, X1, X2);

        if (ind > 0) {
          fprintf(stderr, "Jacobi syuusoku simasen!!\n");
          exit(0);
        }
        // if((A1[0][0] == A1[1][1]) || (A1[1][1] == A1[2][2]) || (A1[0][0] ==
        // A1[2][2]) ){ fprintf(stderr,"Igen Value duplication!!\n"); exit(0);
        // }
        // debug start
        // if(DIVIDENUM == 4){
        // for(i=0;i<3;i++){
        // for(m=0;m<3;m++){
        // fprintf(stderr,"%f ",A1[i][m]);
        // }
        // fprintf(stderr,"\n");
        // }
        //}
        // debug end
        IGENMAX = -999999999999999999.9;
        for (int i = 0; i < 3; i++) {
          if (A1[i][i] > IGENMAX) {
            IGENMAX = A1[i][i];
            Y = i;
          }
        }
        if (bun == 2) {
          PT[MEN][1].MAXDISTANCE = IGENMAX * (double)PT[MEN][1].INDEXNUM;
        }
        if (bun == 3) {
          PT[MEN][1].MAXDISTANCE = IGENMAX;
        }
        if (bun == 4 || bun == 5) {
          for (int i = 0; i < hsize * vsize; i++) {
            if (INDEX[i] == PT[MEN][1].INDEXNO) {
              RR[i] = YIN[i] - U_R8;
              GG[i] = UIN[i] - U_G8;
              BB[i] = VIN[i] - U_B8;
            }
          }
          V[0] = X1[0][Y];
          V[1] = X1[1][Y];
          V[2] = X1[2][Y];
          for (int i = 0; i < hsize * vsize; i++) {
            kaiten(V, RR[i], GG[i], BB[i], RRR + i, GGG + i, BBB + i);
          }
          // for(i=0,l=0;i<hsize*vsize;i++){
          // if(index3[i]!=-1){
          // iRRRR[l] = (int)(RRR[i]+0.5);
          // l++;
          //}
          //}

          // int GASOSUU = hsize*vsize-k;
          // if(omh == 3 || omh == 4){
          // double *RRRRRR,*GGGGGG,*BBBBBB;
          //            RRRRRR = (double*)malloc(sizeof(double)*hsize*vsize);
          // GGGGGG = (double*)malloc(sizeof(double)*hsize*vsize);
          // BBBBBB = (double*)malloc(sizeof(double)*hsize*vsize);
          int l2 = 0;
          for (int i = 0; i < hsize * vsize; i++) {
            if (INDEX[i] == PT[MEN][1].INDEXNO) {
              RRRRRR[l2] = RRR[i];
              GGGGGG[l2] = GGG[i];
              BBBBBB[l2] = BBB[i];
              l2++;
            }
          }
          MAXX = -99.9e64;
          MINN = 99.9e64;
          for (int i = 0; i < l2; i++) {
            if (MAXX < RRRRRR[i]) {
              MAXX = RRRRRR[i];
            }
            if (MINN > RRRRRR[i]) {
              MINN = RRRRRR[i];
            }
          }
          if (bun == 4) {
            PT[MEN][1].MAXDISTANCE = MAXX - MINN;
          }
          if (bun == 5) {
            PT[MEN][1].MAXDISTANCE =
                (MAXX - MINN) * (double)PT[MEN][1].INDEXNUM;
          }
        }

      } // bun == 1 if end
    } else {
      PT[MEN][1].MAXDISTANCE = 0.0;
    }
    // max distance div end

    //分割しなかったボックスのパラメータを違うMEN*に保存
    if (NUM == 0) {
      for (int i = 2; i <= DIVIDENUM - 1; i++) {
        PT[MEN][i].INDEXNO = PT[NMEN][i - 1].INDEXNO;
        PT[MEN][i].INDEXNUM = PT[NMEN][i - 1].INDEXNUM;
        PT[MEN][i].MAXDISTANCE = PT[NMEN][i - 1].MAXDISTANCE;
      }
    }
    if ((NUM >= 1) && (NUM <= (DIVIDENUM - 3))) {
      for (int i = 2; i <= NUM + 1; i++) {
        // i i-2
        PT[MEN][i].INDEXNO = PT[NMEN][i - 2].INDEXNO;
        PT[MEN][i].INDEXNUM = PT[NMEN][i - 2].INDEXNUM;
        PT[MEN][i].MAXDISTANCE = PT[NMEN][i - 2].MAXDISTANCE;
      }
      for (int i = NUM + 2; i <= DIVIDENUM - 1; i++) {
        // i i-1
        PT[MEN][i].INDEXNO = PT[NMEN][i - 1].INDEXNO;
        PT[MEN][i].INDEXNUM = PT[NMEN][i - 1].INDEXNUM;
        PT[MEN][i].MAXDISTANCE = PT[NMEN][i - 1].MAXDISTANCE;
      }
    }
    if (NUM == (DIVIDENUM - 2)) {
      for (int i = 2; i <= NUM + 1; i++) {
        // i i-2
        PT[MEN][i].INDEXNO = PT[NMEN][i - 2].INDEXNO;
        PT[MEN][i].INDEXNUM = PT[NMEN][i - 2].INDEXNUM;
        PT[MEN][i].MAXDISTANCE = PT[NMEN][i - 2].MAXDISTANCE;
      }
    }
    // debug start
    if (DIVIDENUM == 256) {
      for (int i = 0; i < DIVIDENUM; i++) {
        fprintf(stderr, "%d %d %d\n", i, PT[MEN][i].INDEXNO,
                PT[MEN][i].INDEXNUM);
      }
    }
    // debug end
  } // while end
  // debug start
  // for(i=0;i<256;i++){
  // fprintf(stderr,"%d %d %d\n",i,PT[MEN][i].INDEXNO,PT[MEN][i].INDEXNUM);
  // }
  // debug end

  // 256個のパラメータ空間の重心をだしパレットとして登録

  double Y_JYUSHIN[IROSUU];
  double U_JYUSHIN[IROSUU];
  double V_JYUSHIN[IROSUU];
  for (int i = 0; i < IROSUU; i++) {
    Y_JYUSHIN[i] = 0.0;
    U_JYUSHIN[i] = 0.0;
    V_JYUSHIN[i] = 0.0;
  }
  int w;
  w = 0;
  for (int i = 0; i < hsize * vsize; i++) {

    if (index3[i] != -1) {
      for (int j = 0; j < IROSUU; j++) {
        if (INDEX[i] == j) {
          Y_JYUSHIN[j] += YIN[i];
          U_JYUSHIN[j] += UIN[i];
          V_JYUSHIN[j] += VIN[i];
        }
      }
      w++;
    }
  }
  printf("edge denai gasosuu = %d\n", w);
  int y;
  y = 0;
  for (int i = 0; i < IROSUU; i++) {
    y += PT[MEN][i].INDEXNUM;
  }
  printf("palet sou gasosuu=%d\n", y);

  for (m = 0; m < IROSUU; m++) {
    int k = PT[MEN][m].INDEXNO;
    if (PT[MEN][m].INDEXNUM != 0) {
      Y_JYUSHIN[k] /= (double)PT[MEN][m].INDEXNUM;
      U_JYUSHIN[k] /= (double)PT[MEN][m].INDEXNUM;
      V_JYUSHIN[k] /= (double)PT[MEN][m].INDEXNUM;
    } else {
      Y_JYUSHIN[k] = 255.0;
      U_JYUSHIN[k] = 255.0; // 128.0;
      V_JYUSHIN[k] = 255.0; // 128.0;
    }
  }
  // debug start
  // for(i=0;i<IROSUU;i++){
  // fprintf(stderr,"PALET%d %d %d\n",Y_JYUSHIN[i],U_JYUSHIN[i],V_JYUSHIN[i]);
  // }
  // while(1);
  // debug end
  // for(i=0;i<IROSUU;i++){
  // REDUCE_R[i] = (unsigned char)(R_JYUSHIN[i]);
  // REDUCE_G[i] = (unsigned char)(G_JYUSHIN[i]);
  // REDUCE_B[i] = (unsigned char)(B_JYUSHIN[i]);
  // }
  // PALETGAZOU = (unsigned char*)malloc(sizeof(unsigned char)*hsize*vsize);
  // int q;
  int s;
  p = 0;
  s = 0;
  if (kmon == 1) {
    while (1) {
      p++;
      fprintf(stderr, "%d iter\n", p);
      for (int i = 0; i < hsize * vsize; i++) {
        if (index3[i] != -1) {
          double Z[IROSUU];
          for (int j = 0; j < IROSUU; j++) {
            Z[j] = (YIN[i] - Y_JYUSHIN[j]) * (YIN[i] - Y_JYUSHIN[j]) +
                   (UIN[i] - U_JYUSHIN[j]) * (UIN[i] - U_JYUSHIN[j]) +
                   (VIN[i] - V_JYUSHIN[j]) * (VIN[i] - V_JYUSHIN[j]);
          }
          double MINZ = DBL_MAX;
          int X = 0;
          for (int j = 0; j < IROSUU; j++) {
            if (Z[j] < MINZ) {
              MINZ = Z[j];
              X = j;
            }
          }
          PALETGAZOU[i] = X;
        }
      }
      // debug start
      int KARI[IROSUU];
      if (1) {
        for (int i = 0; i < IROSUU; i++) {
          KARI[i] = 0;
        }
        for (int i = 0; i < hsize * vsize; i++) {
          KARI[PALETGAZOU[i]]++;
        }
        for (int i = 0; i < IROSUU; i++) {
          fprintf(stderr, "i=%d %d\n", i, KARI[i]);
        }
      }
      // debug end

      // kmeans法
      int l2 = 0;
      for (int i = 0; i < IROSUU; i++) {
        SUM_R[i] = 0.0;
        SUM_G[i] = 0.0;
        SUM_B[i] = 0.0;
      }
      for (int j = 0; j < IROSUU; j++) {
        int k = 0;
        for (int i = 0; i < hsize * vsize; i++) {
          if (index3[i] != -1) {
            if ((PALETGAZOU[i]) == j) {
              SUM_R[j] += (YIN[i]);
              SUM_G[j] += (UIN[i]);
              SUM_B[j] += (VIN[i]);
              k++;
            }
          }
        }
        if (k != 0) {
          HEIKIN_R[j] = SUM_R[j] / (double)k;
          HEIKIN_G[j] = SUM_G[j] / (double)k;
          HEIKIN_B[j] = SUM_B[j] / (double)k;
        } else if (k == 0) {
          l2++;
          HEIKIN_R[j] = (double)((rand()) % ((int)(gmult * 256)));
          HEIKIN_G[j] = (double)((rand()) % ((int)(bmult * 256)));
          HEIKIN_B[j] = (double)((rand()) % ((int)(rmult * 256)));
          // debug start
          printf("RAND de dasita palet = %f %f %f\n", HEIKIN_R[j], HEIKIN_G[j],
                 HEIKIN_B[j]);
          // debug end
        }
      }
      if (l2 == 0) {
        if (s == 0) {
          // q = p + 50;
          s = 1;
        }
      }

      bool ALL_SAME = true;
      for (int i = 0; ALL_SAME && i < 256; i++) {
        ALL_SAME = ALL_SAME && ((int)HEIKIN_R[i] == (int)Y_JYUSHIN[i]);
      }
      for (int i = 0; ALL_SAME && i < 256; i++) {
        ALL_SAME = ALL_SAME && ((int)HEIKIN_G[i] == (int)U_JYUSHIN[i]);
      }
      for (int i = 0; ALL_SAME && i < 256; i++) {
        ALL_SAME = ALL_SAME && ((int)HEIKIN_B[i] == (int)V_JYUSHIN[i]);
      }
      if (ALL_SAME) {
        fprintf(stderr, "%d iter\n", p);
        break;
      }

      for (int i = 0; i < IROSUU; i++) {
        Y_JYUSHIN[i] = HEIKIN_R[i];
        U_JYUSHIN[i] = HEIKIN_G[i];
        V_JYUSHIN[i] = HEIKIN_B[i];
      }
      // if(p==100){
      // break;
      // }
      // if(l==0 && q==p){
      // break;
      // }
    } // kmean loop end
  }
  printf("hoge3\n");

  Dithering(hsize * vsize, &YIN[0], &UIN[0], &VIN[0], IROSUU, &Y_JYUSHIN[0],
            &U_JYUSHIN[0], &V_JYUSHIN[0], &PALETGAZOU[0]);

  double *Y_JYUSHIN1;
  double *U_JYUSHIN1;
  double *V_JYUSHIN1;
  Y_JYUSHIN1 = (double *)malloc(sizeof(double) * IROSUU);
  U_JYUSHIN1 = (double *)malloc(sizeof(double) * IROSUU);
  V_JYUSHIN1 = (double *)malloc(sizeof(double) * IROSUU);
  for (int i = 0; i < IROSUU; i++) {
    Y_JYUSHIN1[i] = Y_JYUSHIN[i];
    U_JYUSHIN1[i] = U_JYUSHIN[i];
    V_JYUSHIN1[i] = V_JYUSHIN[i];
  }

  // YUVパレットが完成
  //これをＲＧＢに変換
  int IREDUCE_R[IROSUU];
  int IREDUCE_G[IROSUU];
  int IREDUCE_B[IROSUU];
  for (int i = 0; i < IROSUU; i++) {
    // U_JYUSHIN[i] = 255.0*bmult*pow(U_JYUSHIN[i]/(bmult*255.0),1.0/bgamma);
    // V_JYUSHIN[i] = 255.0*rmult*pow(V_JYUSHIN[i]/(rmult*255.0),1.0/rgamma);
  }
  if (cs == 4 || cs == 9) {
    for (int i = 0; i < IROSUU; i++) {
      KARIU[i] =
          U_JYUSHIN[i] / ((float)(bmult)) * /*cos(IRO_THIETA)*/ 1.076908585 -
          V_JYUSHIN[i] / ((float)(rmult)) * /*sin(IRO_THIETA)*/ 0.122197951;
      KARIV[i] =
          -U_JYUSHIN[i] / ((float)(bmult)) * /*sin(IRO_THIETA)*/ 0.105951892 +
          V_JYUSHIN[i] / ((float)(rmult)) * /*cos(IRO_THIETA)*/ 0.844366606;
      U_JYUSHIN[i] = KARIU[i];
      V_JYUSHIN[i] = KARIV[i];
      Y_JYUSHIN[i] /= gmult;
    }
  } else {
    for (int i = 0; i < IROSUU; i++) {
      KARIU[i] = U_JYUSHIN[i] / ((float)(bmult)) * cos(IRO_THIETA) +
                 V_JYUSHIN[i] / ((float)(rmult)) * sin(IRO_THIETA);
      KARIV[i] = -U_JYUSHIN[i] / ((float)(bmult)) * sin(IRO_THIETA) +
                 V_JYUSHIN[i] / ((float)(rmult)) * cos(IRO_THIETA);
      U_JYUSHIN[i] = KARIU[i];
      V_JYUSHIN[i] = KARIV[i];
      Y_JYUSHIN[i] /= gmult;
    }
  }
  for (int i = 0; i < IROSUU; i++) {
    Y_JYUSHIN[i] = 255.0 * pow(Y_JYUSHIN[i] / 255.0, 1.0 / ygamma);
  }
  double KAARII;
  double *FREDUCE_R, *FREDUCE_G, *FREDUCE_B;
  FREDUCE_R = (double *)malloc(sizeof(double) * IROSUU);
  FREDUCE_G = (double *)malloc(sizeof(double) * IROSUU);
  FREDUCE_B = (double *)malloc(sizeof(double) * IROSUU);

  if (cs == 0) {
    for (int i = 0; i < IROSUU; i++) {
      FREDUCE_R[i] =
          255.0 *
          pow(((float)V_JYUSHIN[i]) / 255.0,
              1.0 /
                  rgamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))+1.40200*((float)(V_JYUSHIN[i])/((float)(RMULT))-128.0))+0.5);
      FREDUCE_G[i] =
          255.0 *
          pow(((float)Y_JYUSHIN[i]) / 255.0,
              1.0 /
                  gamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))-0.34414*((float)(U_JYUSHIN[i])/((float)(BMULT))-128.0)-0.71414*((float)(V_JYUSHIN[i])/((float)(RMULT))-128.0))+0.5);
      FREDUCE_B[i] =
          255.0 *
          pow(((float)U_JYUSHIN[i]) / 255.0,
              1.0 /
                  bgamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))+1.77200*((float)(U_JYUSHIN[i]/((float)(BMULT)))-128.0))+0.5);
      IREDUCE_R[i] = (int)(FREDUCE_R[i] + 0.5);
      IREDUCE_G[i] = (int)(FREDUCE_G[i] + 0.5);
      IREDUCE_B[i] = (int)(FREDUCE_B[i] + 0.5);
    }
  } else if (cs == 1) {
    for (int i = 0; i < IROSUU; i++) {
      LABtoRGB((float)V_JYUSHIN[i], (float)Y_JYUSHIN[i], (float)U_JYUSHIN[i],
               FREDUCE_R + i, FREDUCE_G + i, FREDUCE_B + i);
      FREDUCE_R[i] =
          255.0 *
          pow(((float)FREDUCE_R[i]) / 255.0,
              1.0 /
                  rgamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))+1.40200*((float)(V_JYUSHIN[i])/((float)(RMULT))-128.0))+0.5);
      FREDUCE_G[i] =
          255.0 *
          pow(((float)FREDUCE_G[i]) / 255.0,
              1.0 /
                  gamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))-0.34414*((float)(U_JYUSHIN[i])/((float)(BMULT))-128.0)-0.71414*((float)(V_JYUSHIN[i])/((float)(RMULT))-128.0))+0.5);
      FREDUCE_B[i] =
          255.0 *
          pow(((float)FREDUCE_B[i]) / 255.0,
              1.0 /
                  bgamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))+1.77200*((float)(U_JYUSHIN[i]/((float)(BMULT)))-128.0))+0.5);
      IREDUCE_R[i] = (int)(FREDUCE_R[i] + 0.5);
      IREDUCE_G[i] = (int)(FREDUCE_G[i] + 0.5);
      IREDUCE_B[i] = (int)(FREDUCE_B[i] + 0.5);
    }
  } else if (cs == 2) {
    for (int i = 0; i < IROSUU; i++) {
      LuvtoRGB((float)V_JYUSHIN[i], (float)Y_JYUSHIN[i], (float)U_JYUSHIN[i],
               FREDUCE_R + i, FREDUCE_G + i, FREDUCE_B + i);
      FREDUCE_R[i] =
          255.0 *
          pow(((float)FREDUCE_R[i]) / 255.0,
              1.0 /
                  rgamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))+1.40200*((float)(V_JYUSHIN[i])/((float)(RMULT))-128.0))+0.5);
      FREDUCE_G[i] =
          255.0 *
          pow(((float)FREDUCE_G[i]) / 255.0,
              1.0 /
                  gamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))-0.34414*((float)(U_JYUSHIN[i])/((float)(BMULT))-128.0)-0.71414*((float)(V_JYUSHIN[i])/((float)(RMULT))-128.0))+0.5);
      FREDUCE_B[i] =
          255.0 *
          pow(((float)FREDUCE_B[i]) / 255.0,
              1.0 /
                  bgamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))+1.77200*((float)(U_JYUSHIN[i]/((float)(BMULT)))-128.0))+0.5);
      IREDUCE_R[i] = (int)(FREDUCE_R[i] + 0.5);
      IREDUCE_G[i] = (int)(FREDUCE_G[i] + 0.5);
      IREDUCE_B[i] = (int)(FREDUCE_B[i] + 0.5);
    }
  } else if (cs == 3 || cs == 9) {
    for (int i = 0; i < IROSUU; i++) {
      // IREDUCE_R[i] =
      // (int)((((float)(Y_JYUSHIN[i]))+1.40200*((float)(V_JYUSHIN[i])-128.0))+0.5);
      // IREDUCE_G[i] =
      // (int)((((float)(Y_JYUSHIN[i]))-0.34414*((float)(U_JYUSHIN[i])-128.0)-0.71414*((float)(V_JYUSHIN[i])-128.0))+0.5);
      // IREDUCE_B[i] =
      // (int)((((float)(Y_JYUSHIN[i]))+1.77200*((float)(U_JYUSHIN[i])-128.0))+0.5);
      FREDUCE_R[i] = ((((float)(Y_JYUSHIN[i])) +
                       1.4018452447105 * ((float)(V_JYUSHIN[i]) - 128.0)));
      FREDUCE_G[i] = ((((float)(Y_JYUSHIN[i])) -
                       0.34511790321824 * ((float)(U_JYUSHIN[i]) - 128.0) -
                       0.71429038997809 * ((float)(V_JYUSHIN[i]) - 128.0)));
      FREDUCE_B[i] = ((((float)(Y_JYUSHIN[i])) +
                       1.7710177314704 * ((float)(U_JYUSHIN[i]) - 128.0)));
      FREDUCE_R[i] =
          255.0 *
          pow(((float)FREDUCE_R[i]) / 255.0,
              1.0 /
                  rgamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))+1.40200*((float)(V_JYUSHIN[i])/((float)(RMULT))-128.0))+0.5);
      FREDUCE_G[i] =
          255.0 *
          pow(((float)FREDUCE_G[i]) / 255.0,
              1.0 /
                  gamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))-0.34414*((float)(U_JYUSHIN[i])/((float)(BMULT))-128.0)-0.71414*((float)(V_JYUSHIN[i])/((float)(RMULT))-128.0))+0.5);
      FREDUCE_B[i] =
          255.0 *
          pow(((float)FREDUCE_B[i]) / 255.0,
              1.0 /
                  bgamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))+1.77200*((float)(U_JYUSHIN[i]/((float)(BMULT)))-128.0))+0.5);
      IREDUCE_R[i] = (int)(FREDUCE_R[i] + 0.5);
      IREDUCE_G[i] = (int)(FREDUCE_G[i] + 0.5);
      IREDUCE_B[i] = (int)(FREDUCE_B[i] + 0.5);
    }
    //} else if(cs == 4){
    // for(i=0;i<IROSUU;i++){
    // IREDUCE_R[i] =
    // (int)((((float)(Y_JYUSHIN[i])*0.9998247134019967)-0.000044452970991771824*((float)(U_JYUSHIN[i]))+1.7528659800326482*((float)(V_JYUSHIN[i])))+0.5);
    // IREDUCE_G[i] =
    // (int)((((float)(Y_JYUSHIN[i])*1.0000893144198046)-0.4320813565838843*((float)(U_JYUSHIN[i]))-0.89314419804726*((float)(V_JYUSHIN[i])))+0.5);
    // IREDUCE_B[i] =
    // (int)((((float)(Y_JYUSHIN[i])*1.0000000115762946)+2.213731098385467*((float)(U_JYUSHIN[i]))-0.00011576294529099052*((float)(V_JYUSHIN[i])))+0.5);
    ////            RR=
    /// n1*0.9998247134019967+n2*-0.000044452970991771824+n3* 1.7528659800326482;
    ////            GG= n1*1.0000893144198046+n2*-0.4320813565838843+
    /// n3*-0.89314419804726; /            BB=
    /// n1*1.0000000115762946+n2* 2.213731098385467+ n3*-0.00011576294529099052;
    //}
    //}
  } else if (cs == 4) {
    for (int i = 0; i < IROSUU; i++) {
      KAARII = ((((float)(Y_JYUSHIN[i]) * 0.9998247134019967) -
                 0.000044452970991771824 * ((float)(U_JYUSHIN[i])) +
                 1.7528659800326482 * ((float)(V_JYUSHIN[i]))));
      if (KAARII / 255.0 >= 0.0540269630587776405948631399435028) {
        KAARII = pow(KAARII / 255.0, 2.2 / 3.0) * 255.0;
      } else {
        KAARII =
            pow(KAARII / 255.0 / 0.142913647595774056881018286010431, 2.2) *
            255.0;
      }
      FREDUCE_R[i] = KAARII;
      //                r3=colortrim2(R/255/16384);
      //            if (r3>=0.0540269630587776405948631399435028)
      //                rr=Math.Pow(r3,2.2/3.0)*255;
      //            else
      //                rr=Math.Pow(r3/0.142913647595774056881018286010431,2.2)*255;
      KAARII = ((((float)(Y_JYUSHIN[i]) * 1.0000893144198046) -
                 0.4320813565838843 * ((float)(U_JYUSHIN[i])) -
                 0.89314419804726 * ((float)(V_JYUSHIN[i]))));
      if (KAARII / 255.0 >= 0.0540269630587776405948631399435028) {
        KAARII = pow(KAARII / 255.0, 2.2 / 3.0) * 255.0;
      } else {
        KAARII =
            pow(KAARII / 255.0 / 0.142913647595774056881018286010431, 2.2) *
            255.0;
      }
      FREDUCE_G[i] = KAARII;
      KAARII = ((((float)(Y_JYUSHIN[i]) * 1.0000000115762946) +
                 2.213731098385467 * ((float)(U_JYUSHIN[i])) -
                 0.00011576294529099052 * ((float)(V_JYUSHIN[i]))));
      if (KAARII / 255.0 >= 0.0540269630587776405948631399435028) {
        KAARII = pow(KAARII / 255.0, 2.2 / 3.0) * 255.0;
      } else {
        KAARII =
            pow(KAARII / 255.0 / 0.142913647595774056881018286010431, 2.2) *
            255.0;
      }
      FREDUCE_B[i] = KAARII;

      FREDUCE_R[i] =
          255.0 *
          pow(((float)FREDUCE_R[i]) / 255.0,
              1.0 /
                  rgamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))+1.40200*((float)(V_JYUSHIN[i])/((float)(RMULT))-128.0))+0.5);
      FREDUCE_G[i] =
          255.0 *
          pow(((float)FREDUCE_G[i]) / 255.0,
              1.0 /
                  gamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))-0.34414*((float)(U_JYUSHIN[i])/((float)(BMULT))-128.0)-0.71414*((float)(V_JYUSHIN[i])/((float)(RMULT))-128.0))+0.5);
      FREDUCE_B[i] =
          255.0 *
          pow(((float)FREDUCE_B[i]) / 255.0,
              1.0 /
                  bgamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))+1.77200*((float)(U_JYUSHIN[i]/((float)(BMULT)))-128.0))+0.5);
      IREDUCE_R[i] = (int)(FREDUCE_R[i] + 0.5);
      IREDUCE_G[i] = (int)(FREDUCE_G[i] + 0.5);
      IREDUCE_B[i] = (int)(FREDUCE_B[i] + 0.5);
    }
  } else if (cs == 5) {
    for (int i = 0; i < IROSUU; i++) {
      // IREDUCE_R[i] =
      // (int)((((float)(Y_JYUSHIN[i]))+1.40200*((float)(V_JYUSHIN[i])-128.0))+0.5);
      // IREDUCE_G[i] =
      // (int)((((float)(Y_JYUSHIN[i]))-0.34414*((float)(U_JYUSHIN[i])-128.0)-0.71414*((float)(V_JYUSHIN[i])-128.0))+0.5);
      // IREDUCE_B[i] =
      // (int)((((float)(Y_JYUSHIN[i]))+1.77200*((float)(U_JYUSHIN[i])-128.0))+0.5);
      //    IREDUCE_R[i] =
      //    (int)((((float)(Y_JYUSHIN[i]))+1.4018452447105*((float)(V_JYUSHIN[i])-128.0))+0.5);
      //    IREDUCE_G[i] =
      //    (int)((((float)(Y_JYUSHIN[i]))-0.34511790321824*((float)(U_JYUSHIN[i])-128.0)-0.71429038997809*((float)(V_JYUSHIN[i])-128.0))+0.5);
      //    IREDUCE_B[i] =
      //    (int)((((float)(Y_JYUSHIN[i]))+1.7710177314704*((float)(U_JYUSHIN[i])-128.0))+0.5);
      //        RrIN[i] =
      //        (255.0/(1.0+exp(-aaaa*((double)RIN[i]-127.5)))-255.0/(1.0+exp(-aaaa*(-127.5))))*255.0/(255.0-2.0*255.0/(1.0+exp(-aaaa*(-127.5))));
      KAARII =
          ((((float)(Y_JYUSHIN[i])) +
            /*1.4018452447105*/ 1.40200 * ((float)(V_JYUSHIN[i]) - 128.0)));
      KAARII =
          logl(255.0 /
                   ((KAARII /
                     (255.0 /
                      (255.0 - 2.0 * 255.0 / (1.0 + exp(-aaaa * (-127.5)))))) +
                    255.0 / (1.0 + exp(-aaaa * (-127.5)))) -
               1.0) /
              (-aaaa) +
          127.5;
      FREDUCE_R[i] = (KAARII);
      KAARII =
          ((((float)(Y_JYUSHIN[i])) -
            /*0.34511790321824*/ 0.34414 * ((float)(U_JYUSHIN[i]) - 128.0) -
            /*0.71429038997809*/ 0.71414 * ((float)(V_JYUSHIN[i]) - 128.0)));
      KAARII =
          logl(255.0 /
                   ((KAARII /
                     (255.0 /
                      (255.0 - 2.0 * 255.0 / (1.0 + exp(-aaaa * (-127.5)))))) +
                    255.0 / (1.0 + exp(-aaaa * (-127.5)))) -
               1.0) /
              (-aaaa) +
          127.5;
      FREDUCE_G[i] = (KAARII);
      KAARII =
          ((((float)(Y_JYUSHIN[i])) +
            /*1.7710177314704*/ 1.77200 * ((float)(U_JYUSHIN[i]) - 128.0)) +
           0.5);
      KAARII =
          logl(255.0 /
                   ((KAARII /
                     (255.0 /
                      (255.0 - 2.0 * 255.0 / (1.0 + exp(-aaaa * (-127.5)))))) +
                    255.0 / (1.0 + exp(-aaaa * (-127.5)))) -
               1.0) /
              (-aaaa) +
          127.5;
      FREDUCE_B[i] = (KAARII);
      FREDUCE_R[i] =
          255.0 *
          pow(((float)FREDUCE_R[i]) / 255.0,
              1.0 /
                  rgamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))+1.40200*((float)(V_JYUSHIN[i])/((float)(RMULT))-128.0))+0.5);
      FREDUCE_G[i] =
          255.0 *
          pow(((float)FREDUCE_G[i]) / 255.0,
              1.0 /
                  gamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))-0.34414*((float)(U_JYUSHIN[i])/((float)(BMULT))-128.0)-0.71414*((float)(V_JYUSHIN[i])/((float)(RMULT))-128.0))+0.5);
      FREDUCE_B[i] =
          255.0 *
          pow(((float)FREDUCE_B[i]) / 255.0,
              1.0 /
                  bgamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))+1.77200*((float)(U_JYUSHIN[i]/((float)(BMULT)))-128.0))+0.5);
      IREDUCE_R[i] = (int)(FREDUCE_R[i] + 0.5);
      IREDUCE_G[i] = (int)(FREDUCE_G[i] + 0.5);
      IREDUCE_B[i] = (int)(FREDUCE_B[i] + 0.5);
    }
  } else if (cs == 6) {
    for (int i = 0; i < IROSUU; i++) {
      KAARII = V_JYUSHIN[i];
      KAARII =
          logl(255.0 /
                   ((KAARII /
                     (255.0 /
                      (255.0 - 2.0 * 255.0 / (1.0 + exp(-aaaa * (-127.5)))))) +
                    255.0 / (1.0 + exp(-aaaa * (-127.5)))) -
               1.0) /
              (-aaaa) +
          127.5;
      FREDUCE_R[i] =
          ((float)
               KAARII); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))+1.40200*((float)(V_JYUSHIN[i])/((float)(RMULT))-128.0))+0.5);
      KAARII = Y_JYUSHIN[i];
      KAARII =
          logl(255.0 /
                   ((KAARII /
                     (255.0 /
                      (255.0 - 2.0 * 255.0 / (1.0 + exp(-aaaa * (-127.5)))))) +
                    255.0 / (1.0 + exp(-aaaa * (-127.5)))) -
               1.0) /
              (-aaaa) +
          127.5;
      FREDUCE_G[i] =
          ((float)
               KAARII); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))-0.34414*((float)(U_JYUSHIN[i])/((float)(BMULT))-128.0)-0.71414*((float)(V_JYUSHIN[i])/((float)(RMULT))-128.0))+0.5);
      KAARII = U_JYUSHIN[i];
      KAARII =
          logl(255.0 /
                   ((KAARII /
                     (255.0 /
                      (255.0 - 2.0 * 255.0 / (1.0 + exp(-aaaa * (-127.5)))))) +
                    255.0 / (1.0 + exp(-aaaa * (-127.5)))) -
               1.0) /
              (-aaaa) +
          127.5;
      FREDUCE_B[i] =
          ((float)
               KAARII); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))+1.77200*((float)(U_JYUSHIN[i]/((float)(BMULT)))-128.0))+0.5);
      FREDUCE_R[i] =
          255.0 *
          pow(((float)FREDUCE_R[i]) / 255.0,
              1.0 /
                  rgamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))+1.40200*((float)(V_JYUSHIN[i])/((float)(RMULT))-128.0))+0.5);
      FREDUCE_G[i] =
          255.0 *
          pow(((float)FREDUCE_G[i]) / 255.0,
              1.0 /
                  gamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))-0.34414*((float)(U_JYUSHIN[i])/((float)(BMULT))-128.0)-0.71414*((float)(V_JYUSHIN[i])/((float)(RMULT))-128.0))+0.5);
      FREDUCE_B[i] =
          255.0 *
          pow(((float)FREDUCE_B[i]) / 255.0,
              1.0 /
                  bgamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))+1.77200*((float)(U_JYUSHIN[i]/((float)(BMULT)))-128.0))+0.5);
      IREDUCE_R[i] = (int)(FREDUCE_R[i] + 0.5);
      IREDUCE_G[i] = (int)(FREDUCE_G[i] + 0.5);
      IREDUCE_B[i] = (int)(FREDUCE_B[i] + 0.5);
    }

  } else if (cs == 7) {
    for (int i = 0; i < IROSUU; i++) {
      KAARII = V_JYUSHIN[i];
      if (KAARII / 255.0 >= 0.0540269630587776405948631399435028) {
        KAARII = pow(KAARII / 255.0, 2.2 / 3.0) * 255.0;
      } else {
        KAARII =
            pow(KAARII / 255.0 / 0.142913647595774056881018286010431, 2.2) *
            255.0;
      }
      FREDUCE_R[i] = KAARII;
      //                r3=colortrim2(R/255/16384);
      //            if (r3>=0.0540269630587776405948631399435028)
      //                rr=Math.Pow(r3,2.2/3.0)*255;
      //            else
      //                rr=Math.Pow(r3/0.142913647595774056881018286010431,2.2)*255;
      KAARII = Y_JYUSHIN
          [i]; //)*1.0000893144198046)-0.4320813565838843*((float)(U_JYUSHIN[i]))-0.89314419804726*((float)(V_JYUSHIN[i]))));
      if (KAARII / 255.0 >= 0.0540269630587776405948631399435028) {
        KAARII = pow(KAARII / 255.0, 2.2 / 3.0) * 255.0;
      } else {
        KAARII =
            pow(KAARII / 255.0 / 0.142913647595774056881018286010431, 2.2) *
            255.0;
      }
      FREDUCE_G[i] = KAARII;
      KAARII =
          U_JYUSHIN[i]; //))-0.00011576294529099052*((float)(V_JYUSHIN[i]))));
      if (KAARII / 255.0 >= 0.0540269630587776405948631399435028) {
        KAARII = pow(KAARII / 255.0, 2.2 / 3.0) * 255.0;
      } else {
        KAARII =
            pow(KAARII / 255.0 / 0.142913647595774056881018286010431, 2.2) *
            255.0;
      }
      FREDUCE_B[i] = KAARII;
      FREDUCE_R[i] =
          255.0 *
          pow(((float)FREDUCE_R[i]) / 255.0,
              1.0 /
                  rgamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))+1.40200*((float)(V_JYUSHIN[i])/((float)(RMULT))-128.0))+0.5);
      FREDUCE_G[i] =
          255.0 *
          pow(((float)FREDUCE_G[i]) / 255.0,
              1.0 /
                  gamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))-0.34414*((float)(U_JYUSHIN[i])/((float)(BMULT))-128.0)-0.71414*((float)(V_JYUSHIN[i])/((float)(RMULT))-128.0))+0.5);
      FREDUCE_B[i] =
          255.0 *
          pow(((float)FREDUCE_B[i]) / 255.0,
              1.0 /
                  bgamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))+1.77200*((float)(U_JYUSHIN[i]/((float)(BMULT)))-128.0))+0.5);
      IREDUCE_R[i] = (int)(FREDUCE_R[i] + 0.5);
      IREDUCE_G[i] = (int)(FREDUCE_G[i] + 0.5);
      IREDUCE_B[i] = (int)(FREDUCE_B[i] + 0.5);
    }
  } else if (cs == 8) {
    for (int i = 0; i < IROSUU; i++) {
      // IREDUCE_R[i] =
      // (int)((((float)(Y_JYUSHIN[i]))+1.40200*((float)(V_JYUSHIN[i])-128.0))+0.5);
      // IREDUCE_G[i] =
      // (int)((((float)(Y_JYUSHIN[i]))-0.34414*((float)(U_JYUSHIN[i])-128.0)-0.71414*((float)(V_JYUSHIN[i])-128.0))+0.5);
      // IREDUCE_B[i] =
      // (int)((((float)(Y_JYUSHIN[i]))+1.77200*((float)(U_JYUSHIN[i])-128.0))+0.5);
      FREDUCE_R[i] = ((((float)(Y_JYUSHIN[i])) +
                       1.4018452447105 * ((float)(V_JYUSHIN[i]) - 128.0)));
      FREDUCE_G[i] = ((((float)(Y_JYUSHIN[i])) -
                       0.34511790321824 * ((float)(U_JYUSHIN[i]) - 128.0) -
                       0.71429038997809 * ((float)(V_JYUSHIN[i]) - 128.0)));
      FREDUCE_B[i] = ((((float)(Y_JYUSHIN[i])) +
                       1.7710177314704 * ((float)(U_JYUSHIN[i]) - 128.0)));
      KAARII = FREDUCE_R[i];
      if (KAARII / 255.0 >= 0.0540269630587776405948631399435028) {
        KAARII = pow(KAARII / 255.0, 2.2 / 3.0) * 255.0;
      } else {
        KAARII =
            pow(KAARII / 255.0 / 0.142913647595774056881018286010431, 2.2) *
            255.0;
      }
      FREDUCE_R[i] = KAARII;
      //                r3=colortrim2(R/255/16384);
      //            if (r3>=0.0540269630587776405948631399435028)
      //                rr=Math.Pow(r3,2.2/3.0)*255;
      //            else
      //                rr=Math.Pow(r3/0.142913647595774056881018286010431,2.2)*255;
      KAARII = FREDUCE_G
          [i]; //)*1.0000893144198046)-0.4320813565838843*((float)(U_JYUSHIN[i]))-0.89314419804726*((float)(V_JYUSHIN[i]))));
      if (KAARII / 255.0 >= 0.0540269630587776405948631399435028) {
        KAARII = pow(KAARII / 255.0, 2.2 / 3.0) * 255.0;
      } else {
        KAARII =
            pow(KAARII / 255.0 / 0.142913647595774056881018286010431, 2.2) *
            255.0;
      }
      FREDUCE_G[i] = KAARII;
      KAARII =
          FREDUCE_B[i]; //))-0.00011576294529099052*((float)(V_JYUSHIN[i]))));
      if (KAARII / 255.0 >= 0.0540269630587776405948631399435028) {
        KAARII = pow(KAARII / 255.0, 2.2 / 3.0) * 255.0;
      } else {
        KAARII =
            pow(KAARII / 255.0 / 0.142913647595774056881018286010431, 2.2) *
            255.0;
      }
      FREDUCE_B[i] = KAARII;
      FREDUCE_R[i] =
          255.0 *
          pow(((float)FREDUCE_R[i]) / 255.0,
              1.0 /
                  rgamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))+1.40200*((float)(V_JYUSHIN[i])/((float)(RMULT))-128.0))+0.5);
      FREDUCE_G[i] =
          255.0 *
          pow(((float)FREDUCE_G[i]) / 255.0,
              1.0 /
                  gamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))-0.34414*((float)(U_JYUSHIN[i])/((float)(BMULT))-128.0)-0.71414*((float)(V_JYUSHIN[i])/((float)(RMULT))-128.0))+0.5);
      FREDUCE_B[i] =
          255.0 *
          pow(((float)FREDUCE_B[i]) / 255.0,
              1.0 /
                  bgamma); //(int)((((float)(Y_JYUSHIN[i])/((float)(YMULT)))+1.77200*((float)(U_JYUSHIN[i]/((float)(BMULT)))-128.0))+0.5);
      IREDUCE_R[i] = (int)(FREDUCE_R[i] + 0.5);
      IREDUCE_G[i] = (int)(FREDUCE_G[i] + 0.5);
      IREDUCE_B[i] = (int)(FREDUCE_B[i] + 0.5);
    }
  }

  for (int i = 0; i < IROSUU; i++) {
    if (IREDUCE_R[i] > 255) {
      REDUCE_R[i] = 255;
    } else if (IREDUCE_R[i] < 0) {
      REDUCE_R[i] = 0;
    } else {
      REDUCE_R[i] = (unsigned char)(IREDUCE_R[i]);
    }

    if (IREDUCE_G[i] > 255) {
      REDUCE_G[i] = 255;
    } else if (IREDUCE_G[i] < 0) {
      REDUCE_G[i] = 0;
    } else {
      REDUCE_G[i] = (unsigned char)(IREDUCE_G[i]);
    }

    if (IREDUCE_B[i] > 255) {
      REDUCE_B[i] = 255;
    } else if (IREDUCE_B[i] < 0) {
      REDUCE_B[i] = 0;
    } else {
      REDUCE_B[i] = (unsigned char)(IREDUCE_B[i]);
    }
  }
  // debug start
  // for(i=0;i<IROSUU;i++){
  // fprintf(stderr,"RGBPALET%d %d %d\n",REDUCE_R[i],REDUCE_G[i],REDUCE_B[i]);
  // }
  // debug end

  //出力画像をつくる
  // for(i=0;i<hsize*vsize;i++){
  // ROUT[i] = REDUCE_R[PALETGAZOU[i]];
  // GOUT[i] = REDUCE_G[PALETGAZOU[i]];
  // BOUT[i] = REDUCE_B[PALETGAZOU[i]];
  // }

  // Y_JYUSHIN[i] REDUCE_R[i]

  if (edon == 1) {
    unsigned char *IIRIN;
    unsigned char *IIGIN;
    unsigned char *IIBIN;
    IIRIN = (unsigned char *)malloc(sizeof(unsigned char) * hsize * vsize);
    IIGIN = (unsigned char *)malloc(sizeof(unsigned char) * hsize * vsize);
    IIBIN = (unsigned char *)malloc(sizeof(unsigned char) * hsize * vsize);
    for (int i = 0; i < hsize * vsize; i++) {
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
    double *errorR = nullptr;
    double *errorG = nullptr;
    double *errorB = nullptr;
    if (dither == 0 || dither == 1) {
      errorR = (double *)malloc(sizeof(double) * mx * 2);
      errorG = (double *)malloc(sizeof(double) * mx * 2);
      errorB = (double *)malloc(sizeof(double) * mx * 2);
    } else if (dither == 2 || dither == 3) {
      errorR = (double *)malloc(sizeof(double) * mx * 3);
      errorG = (double *)malloc(sizeof(double) * mx * 3);
      errorB = (double *)malloc(sizeof(double) * mx * 3);
    }
    double r = 0.0;
    double g = 0.0;
    double b = 0.0;
    double re = 0.0;
    double ge = 0.0;
    double be = 0.0;
    int adr = 0;
    if (dither == 0 || dither == 1) {
      for (int i = 0; i < mx * 2; i++) {
        errorR[i] = 0.0;
        errorG[i] = 0.0;
        errorB[i] = 0.0;
      }
    } else if (dither == 2 || dither == 3) {
      for (int i = 0; i < mx * 3; i++) {
        errorR[i] = 0.0;
        errorG[i] = 0.0;
        errorB[i] = 0.0;
      }
    }
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
          r = (double)(RIN[x * vsize + y]) +
              ((errorR[adr]) / (16.0 / (per + 0.0000001)) /*32.0*/);
          g = (double)(GIN[x * vsize + y]) +
              ((errorG[adr]) / (16.0 / (per + 0.0000001)) /*32.0*/);
          b = (double)(BIN[x * vsize + y]) +
              ((errorB[adr]) / (16.0 / (per + 0.0000001)) /*32.0*/);
        } else if (dither == 1) {
          r = (double)(RIN[x * vsize + y]) +
              ((errorR[adr]) / (4.0 / (per + 0.0000001)));
          g = (double)(GIN[x * vsize + y]) +
              ((errorG[adr]) / (4.0 / (per + 0.0000001)));
          b = (double)(BIN[x * vsize + y]) +
              ((errorB[adr]) / (4.0 / (per + 0.0000001)));
        } else if (dither == 2) {
          r = (double)(RIN[x * vsize + y]) +
              ((errorR[adr]) / (42.0 / (per + 0.0000001)));
          g = (double)(GIN[x * vsize + y]) +
              ((errorG[adr]) / (42.0 / (per + 0.0000001)));
          b = (double)(BIN[x * vsize + y]) +
              ((errorB[adr]) / (42.0 / (per + 0.0000001)));
        } else if (dither == 3) {
          r = (double)(RIN[x * vsize + y]) +
              ((errorR[adr]) / (48.0 / (per + 0.0000001)));
          g = (double)(GIN[x * vsize + y]) +
              ((errorG[adr]) / (48.0 / (per + 0.0000001)));
          b = (double)(BIN[x * vsize + y]) +
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
        int bst = FindMinimumDistanceIndex(IROSUU, r, g, b, &REDUCE_R[0],
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

        (IIRIN[x * vsize + y]) = (REDUCE_R[bst]);
        (IIGIN[x * vsize + y]) = (REDUCE_G[bst]);
        (IIBIN[x * vsize + y]) = (REDUCE_B[bst]);
        // (IIRIN[x*vsize+y]) = (int)(REDUCE_R[PALETGAZOU[x*vsize+y]]);
        // (IIGIN[x*vsize+y]) = (int)(REDUCE_G[PALETGAZOU[x*vsize+y]]);
        // (IIBIN[x*vsize+y]) = (int)(REDUCE_B[PALETGAZOU[x*vsize+y]]);
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

    // PALETGAZOU[i]
    // IIRINに画像データがはいっている

    for (int i = 0; i < hsize * vsize; i++) {
      for (int j = 0; j < IROSUU; j++) {
        if ((IIRIN[i] == (REDUCE_R[j])) && (IIGIN[i] == (REDUCE_G[j])) &&
            (IIBIN[i] == (REDUCE_B[j]))) {
          PALETGAZOU[i] = j;
        }
      }
    }
  } // edon if end
  // free(PALETGAZOU);

  // U B,V R
  // Y_JYUSHIN[],REDUCE_G[]

  if (edon == 2) {
    double *IIRIN;
    double *IIGIN;
    double *IIBIN;
    IIRIN = (double *)malloc(sizeof(double) * hsize * vsize);
    IIGIN = (double *)malloc(sizeof(double) * hsize * vsize);
    IIBIN = (double *)malloc(sizeof(double) * hsize * vsize);
    for (int i = 0; i < hsize * vsize; i++) {
      IIRIN[i] = (VIN[i]);
      IIGIN[i] = (YIN[i]);
      IIBIN[i] = (UIN[i]);
    }
    int mx = 0;
    if (dither == 0 || dither == 1) {
      mx = hsize + 2;
    } else if (dither == 2 || dither == 3) {
      mx = hsize + 4;
    }
    double *errorR = nullptr;
    double *errorG = nullptr;
    double *errorB = nullptr;
    if (dither == 0 || dither == 1) {
      errorR = (double *)malloc(sizeof(double) * mx * 2);
      errorG = (double *)malloc(sizeof(double) * mx * 2);
      errorB = (double *)malloc(sizeof(double) * mx * 2);
    } else if (dither == 2 || dither == 3) {
      errorR = (double *)malloc(sizeof(double) * mx * 3);
      errorG = (double *)malloc(sizeof(double) * mx * 3);
      errorB = (double *)malloc(sizeof(double) * mx * 3);
    }
    double r = 0.0;
    double g = 0.0;
    double b = 0.0;
    double re = 0.0;
    double ge = 0.0;
    double be = 0.0;
    int adr = 0;
    if (dither == 0 || dither == 1) {
      for (int i = 0; i < mx * 2; i++) {
        errorR[i] = 0.0;
        errorG[i] = 0.0;
        errorB[i] = 0.0;
      }
    } else if (dither == 2 || dither == 3) {
      for (int i = 0; i < mx * 3; i++) {
        errorR[i] = 0.0;
        errorG[i] = 0.0;
        errorB[i] = 0.0;
      }
    }
    printf("hoge1\n");
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
          r = (double)(VIN[x * vsize + y]) +
              ((errorR[adr]) / (16.0 / (per + 0.0000001)) /*32.0*/);
          g = (double)(YIN[x * vsize + y]) +
              ((errorG[adr]) / (16.0 / (per + 0.0000001)) /*32.0*/);
          b = (double)(UIN[x * vsize + y]) +
              ((errorB[adr]) / (16.0 / (per + 0.0000001)) /*32.0*/);
        } else if (dither == 1) {
          r = (double)(VIN[x * vsize + y]) +
              ((errorR[adr]) / (4.0 / (per + 0.0000001)));
          g = (double)(YIN[x * vsize + y]) +
              ((errorG[adr]) / (4.0 / (per + 0.0000001)));
          b = (double)(UIN[x * vsize + y]) +
              ((errorB[adr]) / (4.0 / (per + 0.0000001)));
        } else if (dither == 2) {
          r = (double)(VIN[x * vsize + y]) +
              ((errorR[adr]) / (42.0 / (per + 0.0000001)));
          g = (double)(YIN[x * vsize + y]) +
              ((errorG[adr]) / (42.0 / (per + 0.0000001)));
          b = (double)(UIN[x * vsize + y]) +
              ((errorB[adr]) / (42.0 / (per + 0.0000001)));
        } else if (dither == 3) {
          r = (double)(VIN[x * vsize + y]) +
              ((errorR[adr]) / (48.0 / (per + 0.0000001)));
          g = (double)(YIN[x * vsize + y]) +
              ((errorG[adr]) / (48.0 / (per + 0.0000001)));
          b = (double)(UIN[x * vsize + y]) +
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
        int bst = FindMinimumDistanceIndex(IROSUU, r, g, b, &V_JYUSHIN1[0],
                                           &Y_JYUSHIN1[0], &U_JYUSHIN1[0]);
        // debug start
        // printf("gosa=%d\n",(int)est);
        // debug end
        re = r - (double)(V_JYUSHIN1[bst]);
        ge = g - (double)(Y_JYUSHIN1[bst]);
        be = b - (double)(U_JYUSHIN1[bst]);
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

        (IIRIN[x * vsize + y]) = (V_JYUSHIN1[bst]);
        (IIGIN[x * vsize + y]) = (Y_JYUSHIN1[bst]);
        (IIBIN[x * vsize + y]) = (U_JYUSHIN1[bst]);
        // (IIRIN[x*vsize+y]) = (int)(REDUCE_R[PALETGAZOU[x*vsize+y]]);
        // (IIGIN[x*vsize+y]) = (int)(REDUCE_G[PALETGAZOU[x*vsize+y]]);
        // (IIBIN[x*vsize+y]) = (int)(REDUCE_B[PALETGAZOU[x*vsize+y]]);
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

    // PALETGAZOU[i]
    // IIRINに画像データがはいっている

    for (int i = 0; i < hsize * vsize; i++) {
      for (int j = 0; j < IROSUU; j++) {
        if (((int)IIRIN[i] == (int)V_JYUSHIN1[j]) &&
            ((int)IIGIN[i] == (int)Y_JYUSHIN1[j]) &&
            ((int)IIBIN[i] == (int)U_JYUSHIN1[j])) {
          PALETGAZOU[i] = j;
          break;
        }
      }
    }
  } // edon if end

} // median cut end
