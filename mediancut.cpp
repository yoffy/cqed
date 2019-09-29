#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
  int m, n, p;
  int MEN;
  int NMEN;
  int NUM;
  // unsigned char *PALETGAZOU;
  // unsigned char REDUCE_R[IROSUU];
  // unsigned char REDUCE_G[IROSUU];
  // unsigned char REDUCE_B[IROSUU];
  double MINZ;
  double Z[IROSUU];
  int X;
  double SUM_R[IROSUU];
  double SUM_G[IROSUU];
  double SUM_B[IROSUU];
  double HEIKIN_R[IROSUU];
  double HEIKIN_G[IROSUU];
  double HEIKIN_B[IROSUU];
  double **A, **A1, **A2, **X1, **X2, eps;
  int i1, ind, ct;
  double TMP_RR;
  double TMP_GG;
  double TMP_BB;
  double TMP_RG;
  double TMP_RB;
  double TMP_GB;
  double IGENMAX;
  int Y;
  int *INDEX;
  double U_R;
  double U_G;
  double U_B;
  double MAXINDEXNUM;
  double *YIN;
  double *UIN;
  double *VIN;
  double *RrIN;
  double *GgIN;
  double *BbIN;
  double *RRR33;
  double PI = atan(1.0) * 4.0;
  // RGBをYUVに変換
  YIN = (double *)malloc(sizeof(double) * hsize * vsize);
  UIN = (double *)malloc(sizeof(double) * hsize * vsize);
  VIN = (double *)malloc(sizeof(double) * hsize * vsize);
  RrIN = (double *)malloc(sizeof(double) * hsize * vsize);
  GgIN = (double *)malloc(sizeof(double) * hsize * vsize);
  BbIN = (double *)malloc(sizeof(double) * hsize * vsize);
  RRR33 = (double *)malloc(sizeof(double) * hsize * vsize);

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
              (unsigned char)(lpfk * (double)YIN[i * vsize + j] +
                              koi * (double)YIN[i * vsize + (j - 1)] +
                              koi * (double)YIN[i * vsize + (j + 1)] + 0.5);
          GIN[i * vsize + j] =
              (unsigned char)(lpfk * (double)UIN[i * vsize + j] +
                              koi * (double)UIN[i * vsize + (j - 1)] +
                              koi * (double)UIN[i * vsize + (j + 1)] + 0.5);
          BIN[i * vsize + j] =
              (unsigned char)(lpfk * (double)VIN[i * vsize + j] +
                              koi * (double)VIN[i * vsize + (j - 1)] +
                              koi * (double)VIN[i * vsize + (j + 1)] + 0.5);
        }
      }
    }
  } // else {
  //    for(i=0;i<hsize*vsize;i++){
  //        RrIN[i] =
  //        (double)RIN[i];//+0.5*(double)YIN[i*vsize+(j-1)]+0.25*(double)YIN[i*vsize+(j+1)];
  //        GgIN[i] =
  //        (double)GIN[i];//+0.5*(double)UIN[i*vsize+(j-1)]+0.25*(double)UIN[i*vsize+(j+1)];
  //        BbIN[i] =
  //        (double)BIN[i];//+0.5*(double)VIN[i*vsize+(j-1)]+0.25*(double)VIN[i*vsize+(j+1)];
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
      *(YIN + i) = *(
          GgIN +
          i); //(0.29891*(float)(*(RIN+i))+0.58661*(float)(*(GIN+i))+0.11448*(float)(*(BIN+i))
              ///*+ 0.5*/);
      *(UIN + i) = *(
          BbIN +
          i); //((-0.16874*(float)(*(RIN+i))-0.33126*(float)(*(GIN+i))+0.50000*(float)(*(BIN+i)))
              //+ 128.0);
      *(VIN + i) = *(
          RrIN +
          i); //(0.50000*(float)(*(RIN+i))-0.41869*(float)(*(GIN+i))-0.08131*(float)(*(BIN+i))
              //+ 128.0);
    }
  } else if (cs == 1) {
    for (int i = 0; i < hsize * vsize; i++) {
      RGBtoLAB(
          RrIN[i], GgIN[i], BbIN[i], VIN + i, YIN + i,
          UIN +
              i); //(0.29891*(float)(*(RIN+i))+0.58661*(float)(*(GIN+i))+0.11448*(float)(*(BIN+i))
                  //+ 0.5);
      //*(UIN+i) =
      //*(BIN+i);//((-0.16874*(float)(*(RIN+i))-0.33126*(float)(*(GIN+i))+0.50000*(float)(*(BIN+i)))
      //+ 128.5);
      //*(VIN+i) =
      //*(RIN+i);//(0.50000*(float)(*(RIN+i))-0.41869*(float)(*(GIN+i))-0.08131*(float)(*(BIN+i))
      //+ 128.5);
    }
  } else if (cs == 2) {
    for (int i = 0; i < hsize * vsize; i++) {
      RGBtoLuv(
          RrIN[i], GgIN[i], BbIN[i], VIN + i, YIN + i,
          UIN +
              i); //(0.29891*(float)(*(RIN+i))+0.58661*(float)(*(GIN+i))+0.11448*(float)(*(BIN+i))
                  //+ 0.5);
      //*(UIN+i) =
      //*(BIN+i);//((-0.16874*(float)(*(RIN+i))-0.33126*(float)(*(GIN+i))+0.50000*(float)(*(BIN+i)))
      //+ 128.5);
      //*(VIN+i) =
      //*(RIN+i);//(0.50000*(float)(*(RIN+i))-0.41869*(float)(*(GIN+i))-0.08131*(float)(*(BIN+i))
      //+ 128.5);
    }
  } else if (cs == 3 || cs == 9) {
    for (int i = 0; i < hsize * vsize; i++) {
      *(YIN + i) =
          (0.29891 * (float)(*(RrIN + i)) + 0.58661 * (float)(*(GgIN + i)) +
           0.11448 * (float)(*(BbIN + i)) /*+ 0.5*/);
      *(UIN + i) =
          ((-0.16874 * (float)(*(RrIN + i)) - 0.33126 * (float)(*(GgIN + i)) +
            0.50000 * (float)(*(BbIN + i))) +
           128.0);
      *(VIN + i) =
          (0.50000 * (float)(*(RrIN + i)) - 0.41869 * (float)(*(GgIN + i)) -
           0.08131 * (float)(*(BbIN + i)) + 128.0);
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

      *(YIN + i) =
          (0.2989 * (float)(*(RrIN + i)) + 0.5866 * (float)(*(GgIN + i)) +
           0.1145 * (float)(*(BbIN + i)) /*+0.5*/);
      *(UIN + i) =
          ((-0.1350 * (float)(*(RrIN + i)) - 0.2650 * (float)(*(GgIN + i)) +
            0.40000 * (float)(*(BbIN + i))));
      *(VIN + i) =
          (0.4000 * (float)(*(RrIN + i)) - 0.3346 * (float)(*(GgIN + i)) -
           0.0653 * (float)(*(BbIN + i)));
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
      *(YIN + i) =
          (0.29891 * (float)(*(RrIN + i)) + 0.58661 * (float)(*(GgIN + i)) +
           0.11448 * (float)(*(BbIN + i)) /*+ 0.5*/);
      *(UIN + i) =
          ((-0.16874 * (float)(*(RrIN + i)) - 0.33126 * (float)(*(GgIN + i)) +
            0.50000 * (float)(*(BbIN + i))) +
           128.0);
      *(VIN + i) =
          (0.50000 * (float)(*(RrIN + i)) - 0.41869 * (float)(*(GgIN + i)) -
           0.08131 * (float)(*(BbIN + i)) + 128.0);
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

      *(YIN + i) = *(
          GgIN +
          i); //(0.29891*(float)(*(RIN+i))+0.58661*(float)(*(GIN+i))+0.11448*(float)(*(BIN+i))
              ///*+ 0.5*/);
      *(UIN + i) = *(
          BbIN +
          i); //((-0.16874*(float)(*(RIN+i))-0.33126*(float)(*(GIN+i))+0.50000*(float)(*(BIN+i)))
              //+ 128.0);
      *(VIN + i) = *(
          RrIN +
          i); //(0.50000*(float)(*(RIN+i))-0.41869*(float)(*(GIN+i))-0.08131*(float)(*(BIN+i))
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

      *(YIN + i) = GgIN[i];
      *(UIN + i) = BbIN
          [i]; //((-0.1350*(float)(*(RrIN+i))-0.2650*(float)(*(GgIN+i))+0.40000*(float)(*(BbIN+i))));
      *(VIN + i) = RrIN
          [i]; //(0.4000*(float)(*(RrIN+i))-0.3346*(float)(*(GgIN+i))-0.0653*(float)(*(BbIN+i)));
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
      *(YIN + i) =
          (0.29891 * (float)(*(RrIN + i)) + 0.58661 * (float)(*(GgIN + i)) +
           0.11448 * (float)(*(BbIN + i)) /*+ 0.5*/);
      *(UIN + i) =
          ((-0.16874 * (float)(*(RrIN + i)) - 0.33126 * (float)(*(GgIN + i)) +
            0.50000 * (float)(*(BbIN + i))) +
           128.0);
      *(VIN + i) =
          (0.50000 * (float)(*(RrIN + i)) - 0.41869 * (float)(*(GgIN + i)) -
           0.08131 * (float)(*(BbIN + i)) + 128.0);
      //}

      //*(YIN+i) = GgIN[i];
      //*(UIN+i) =
      // BbIN[i];//((-0.1350*(float)(*(RrIN+i))-0.2650*(float)(*(GgIN+i))+0.40000*(float)(*(BbIN+i))));
      //*(VIN+i) =
      // RrIN[i];//(0.4000*(float)(*(RrIN+i))-0.3346*(float)(*(GgIN+i))-0.0653*(float)(*(BbIN+i)));
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
    *(YIN + i) = 255.0 * pow(*(YIN + i) / 255.0, ygamma);
  }

  for (int i = 0; i < hsize * vsize; i++) {
    *(YIN + i) *= gmult;
  }
  for (int i = 0; i < hsize * vsize; i++) {
    *(VIN + i) *= rmult;
  }
  // for(i=0;i<hsize*vsize;i++){
  //*(VIN+i) = 255.0*rmult*pow(*(VIN+i)/(rmult*255.0),rgamma);
  //}
  for (int i = 0; i < hsize * vsize; i++) {
    *(UIN + i) *= bmult;
  }
  // for(i=0;i<hsize*vsize;i++){
  //*(UIN+i) = 255.0*bmult*pow(*(UIN+i)/(bmult*255.0),bgamma);
  //}
  // YUVは、0-255 たまに負や255を超える場合が発生するかも、あとで検証
  // debug start
  // for(i=0;i<hsize*vsize;i++){
  // fprintf(stderr," %d %d %d\n",*(YIN+i),*(UIN+i),*(VIN+i));
  // }
  // while(1);
  // debug end
  /****************************************/
  /* */
  /* edge detect */
  /* */
  /****************************************/

  int *IRIN;
  int *IGIN;
  int *IBIN;
  int *HHEDGER;
  int *HHEDGEG;
  int *HHEDGEB;
  int *VVEDGER;
  int *VVEDGEG;
  int *VVEDGEB;
  double EDGERASISAYD;
  IRIN = (int *)malloc(sizeof(int) * hsize * vsize);
  IGIN = (int *)malloc(sizeof(int) * hsize * vsize);
  IBIN = (int *)malloc(sizeof(int) * hsize * vsize);
  int *EDGER, *EDGEG, *EDGEB;
  EDGER = (int *)malloc(sizeof(int) * hsize * vsize);
  EDGEG = (int *)malloc(sizeof(int) * hsize * vsize);
  EDGEB = (int *)malloc(sizeof(int) * hsize * vsize);
  HHEDGER = (int *)malloc(sizeof(int) * hsize * vsize);
  HHEDGEG = (int *)malloc(sizeof(int) * hsize * vsize);
  HHEDGEB = (int *)malloc(sizeof(int) * hsize * vsize);
  VVEDGER = (int *)malloc(sizeof(int) * hsize * vsize);
  VVEDGEG = (int *)malloc(sizeof(int) * hsize * vsize);
  VVEDGEB = (int *)malloc(sizeof(int) * hsize * vsize);
  int *VEDGER, *VEDGEG, *VEDGEB;
  VEDGER = (int *)malloc(sizeof(int) * hsize * vsize);
  VEDGEG = (int *)malloc(sizeof(int) * hsize * vsize);
  VEDGEB = (int *)malloc(sizeof(int) * hsize * vsize);
  for (int i = 0; i < hsize * vsize; i++) {
    *(IRIN + i) = (int)(*(RIN + i));
    *(IGIN + i) = (int)(*(GIN + i));
    *(IBIN + i) = (int)(*(BIN + i));
  }

  int *index3;
  index3 = (int *)malloc(sizeof(int) * hsize * vsize);
  // double EDGERASMAX;
  for (int i = 0; i < hsize * vsize; i++) {
    HHEDGER[i] = 0;
    HHEDGEG[i] = 0;
    HHEDGEB[i] = 0;
    VVEDGER[i] = 0;
    VVEDGEG[i] = 0;
    VVEDGEB[i] = 0;
  }

  int SOBEL[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
  int SOBEL2[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};

  int EDGE_GASOSUU = 0;

  if (vvv == 3 || vvv == 4 || vvv == 5) {
    for (int i = 0; i < vsize; i++) {   // suityoku
      for (int j = 0; j < hsize; j++) { // suihei
        for (int k = 0; k < 3; k++) {   // suityoku
          for (int l = 0; l < 3; l++) { // suihei
            m = i + k - 1;              //-1~1suiityoku
            n = j + l - 1;              //-1~1suihei
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
            HHEDGER[j * vsize + i] +=
                SOBEL2[3 * k + l] * (int)RIN[n * vsize + m];
            HHEDGEG[j * vsize + i] +=
                SOBEL2[3 * k + l] * (int)GIN[n * vsize + m];
            HHEDGEB[j * vsize + i] +=
                SOBEL2[3 * k + l] * (int)BIN[n * vsize + m];
            VVEDGER[j * vsize + i] +=
                SOBEL[3 * k + l] * (int)RIN[n * vsize + m];
            VVEDGEG[j * vsize + i] +=
                SOBEL[3 * k + l] * (int)GIN[n * vsize + m];
            VVEDGEB[j * vsize + i] +=
                SOBEL[3 * k + l] * (int)BIN[n * vsize + m];
            // printf("%d\n",RIN[n*vsize+m]);
          }
        }
      }
    }
  }

  if ((vvv >= 0 && vvv <= 2) || vvv == 6) {
    for (int i = 0; i < vsize; i++) {
      for (int j = 0; j < hsize; j++) {
        int k = j - 1;
        int l = j + 1;
        if (k < 0) {
          k = -k;
        }
        if (k > (hsize - 1)) {
          k = 2 * (hsize - 1) - k;
        }
        if (l < 0) {
          l = -l;
        }
        if (l > (hsize - 1)) {
          l = 2 * (hsize - 1) - l;
        }
        HHEDGER[j * vsize + i] = abs(IRIN[k * vsize + i] - IRIN[l * vsize + i]);
        // EDGER[j*vsize+i]=1;
        // } else {
        // EDGER[j*vsize+i]=0;
        // }
        HHEDGEG[j * vsize + i] = abs(IGIN[k * vsize + i] - IGIN[l * vsize + i]);
        // EDGEG[j*vsize+i]=1;
        // } else {
        // EDGEG[j*vsize+i]=0;
        // }
        HHEDGEB[j * vsize + i] = abs(IBIN[k * vsize + i] - IBIN[l * vsize + i]);
        // EDGEB[j*vsize+i]=1;
        // } else {
        // EDGEB[j*vsize+i]=0;
        // }
        if ((abs(IRIN[k * vsize + i] - IRIN[l * vsize + i])) >= et) {
          EDGER[j * vsize + i] = 1;
        } else {
          EDGER[j * vsize + i] = 0;
        }
        if ((abs(IGIN[k * vsize + i] - IGIN[l * vsize + i])) >= et) {
          EDGEG[j * vsize + i] = 1;
        } else {
          EDGEG[j * vsize + i] = 0;
        }
        if ((abs(IBIN[k * vsize + i] - IBIN[l * vsize + i])) >= et) {
          EDGEB[j * vsize + i] = 1;
        } else {
          EDGEB[j * vsize + i] = 0;
        }
      }
    }

    for (int i = 0; i < hsize; i++) {
      for (int j = 0; j < vsize; j++) {
        int k = j - 1;
        int l = j + 1;
        if (k < 0) {
          k = -k;
        }
        if (l > (vsize - 1)) {
          l = 2 * (vsize - 1) - l;
        }
        VVEDGER[i * vsize + j] = abs(IRIN[i * vsize + l] - IRIN[i * vsize + k]);
        // VEDGER[i*vsize+j] = 1;
        // } else {
        // VEDGER[i*vsize+j] = 0;
        // }
        VVEDGEG[i * vsize + j] = abs(IGIN[i * vsize + l] - IGIN[i * vsize + k]);
        // VEDGEG[i*vsize+j] = 1;
        // } else {
        // VEDGEG[i*vsize+j] = 0;
        // }
        VVEDGEB[i * vsize + j] = abs(IBIN[i * vsize + l] - IBIN[i * vsize + k]);
        // VEDGEB[i*vsize+j] = 1;
        // } else {
        // VEDGEB[i*vsize+j] = 0;
        // }
        if (abs(IRIN[i * vsize + l] - IRIN[i * vsize + k]) >= et) {
          VEDGER[i * vsize + j] = 1;
        } else {
          VEDGER[i * vsize + j] = 0;
        }
        if (abs(IGIN[i * vsize + l] - IGIN[i * vsize + k]) >= et) {
          VEDGEG[i * vsize + j] = 1;
        } else {
          VEDGEG[i * vsize + j] = 0;
        }
        if (abs(IBIN[i * vsize + l] - IBIN[i * vsize + k]) >= et) {
          VEDGEB[i * vsize + j] = 1;
        } else {
          VEDGEB[i * vsize + j] = 0;
        }
      }
    }

    for (int i = 0; i < vsize; i++) {
      for (int j = 0; j < hsize; j++) {
        if (vvv == 0) {
          if ((sqrt((float)HHEDGER[j * vsize + i] *
                        (float)HHEDGER[j * vsize + i] +
                    (float)VVEDGER[j * vsize + i] *
                        (float)VVEDGER[j * vsize + i]) >= et * sqrt(2.0)) ||
              (sqrt((float)HHEDGEG[j * vsize + i] *
                        (float)HHEDGEG[j * vsize + i] +
                    (float)VVEDGEG[j * vsize + i] *
                        (float)VVEDGEG[j * vsize + i]) >= et * sqrt(2.0)) ||
              (sqrt((float)HHEDGEB[j * vsize + i] *
                        (float)HHEDGEB[j * vsize + i] +
                    (float)VVEDGEB[j * vsize + i] *
                        (float)VVEDGEB[j * vsize + i]) >= et * sqrt(2.0))) {
            *(index3 + j * vsize + i) = -1;
            EDGE_GASOSUU++;
          } else {
            *(index3 + j * vsize + i) = 1;
          }
        }
        if (vvv == 6) {
          if (sqrt((((float)HHEDGER[j * vsize + i] *
                         (float)HHEDGER[j * vsize + i] +
                     (float)VVEDGER[j * vsize + i] *
                         (float)VVEDGER[j * vsize + i])) +
                   (((float)HHEDGEG[j * vsize + i] *
                         (float)HHEDGEG[j * vsize + i] +
                     (float)VVEDGEG[j * vsize + i] *
                         (float)VVEDGEG[j * vsize + i])) +
                   (((float)HHEDGEB[j * vsize + i] *
                         (float)HHEDGEB[j * vsize + i] +
                     (float)VVEDGEB[j * vsize + i] *
                         (float)VVEDGEB[j * vsize + i]))) >= et * sqrt(6.0)) {
            *(index3 + j * vsize + i) = -1;
            EDGE_GASOSUU++;
          } else {
            *(index3 + j * vsize + i) = 1;
          }
        }
        if (vvv == 1) {
          if (HHEDGER[j * vsize + i] >= et || VVEDGER[j * vsize + i] >= et ||
              HHEDGEG[j * vsize + i] >= et || VVEDGEG[j * vsize + i] >= et ||
              HHEDGEB[j * vsize + i] >= et || VVEDGEB[j * vsize + i] >= et) {
            *(index3 + j * vsize + i) = -1;
            EDGE_GASOSUU++;
          } else {
            *(index3 + j * vsize + i) = 1;
          }
        }
        if (vvv == 2) {
          if (EDGER[j * vsize + i] == 1 || EDGEG[j * vsize + i] == 1 ||
              EDGEB[j * vsize + i] == 1 || VEDGER[j * vsize + i] == 1 ||
              VEDGEG[j * vsize + i] == 1 || VEDGEB[j * vsize + i] == 1) {
            *(index3 + j * vsize + i) = -1;
            EDGE_GASOSUU++;
          } else {
            *(index3 + j * vsize + i) = 1;
          }
        }
      }
    }

  } // cnt if end

  if (vvv == 3 || vvv == 5) {
    for (int i = 0; i < vsize; i++) {
      for (int j = 0; j < hsize; j++) {
        // printf("%f\n",sqrt((float)HHEDGER[j*vsize+i]*(float)HHEDGER[j*vsize+i]
        // + (float)VVEDGER[j*vsize+i] * (float)VVEDGER[j*vsize+i]));
        // printf("%f\n",sqrt((float)HHEDGEG[j*vsize+i]*(float)HHEDGEG[j*vsize+i]
        // + (float)VVEDGEG[j*vsize+i] * (float)VVEDGEG[j*vsize+i]));
        // printf("%f\n",sqrt((float)HHEDGEB[j*vsize+i]*(float)HHEDGEB[j*vsize+i]
        // + (float)VVEDGEB[j*vsize+i] * (float)VVEDGEB[j*vsize+i]));

        if (vvv == 3) {
          if ((sqrt((float)HHEDGER[j * vsize + i] *
                        (float)HHEDGER[j * vsize + i] +
                    (float)VVEDGER[j * vsize + i] *
                        (float)VVEDGER[j * vsize + i]) >=
               et * 4.0 * sqrt(2.0)) ||
              (sqrt((float)HHEDGEG[j * vsize + i] *
                        (float)HHEDGEG[j * vsize + i] +
                    (float)VVEDGEG[j * vsize + i] *
                        (float)VVEDGEG[j * vsize + i]) >=
               et * 4.0 * sqrt(2.0)) ||
              (sqrt((float)HHEDGEB[j * vsize + i] *
                        (float)HHEDGEB[j * vsize + i] +
                    (float)VVEDGEB[j * vsize + i] *
                        (float)VVEDGEB[j * vsize + i]) >=
               et * 4.0 * sqrt(2.0))) {
            *(index3 + j * vsize + i) = -1;
            EDGE_GASOSUU++;
          } else {
            *(index3 + j * vsize + i) = 1;
          }
        }

        if (vvv == 5) {
          //            printf("%f\n",sqrt(((float)HHEDGER[j*vsize+i]*(float)HHEDGER[j*vsize+i]
          //            + (float)VVEDGER[j*vsize+i] * (float)VVEDGER[j*vsize+i])
          //                +((float)HHEDGEG[j*vsize+i]*(float)HHEDGEG[j*vsize+i]
          //                + (float)VVEDGEG[j*vsize+i] *
          //                (float)VVEDGEG[j*vsize+i])
          //            +((float)HHEDGEB[j*vsize+i]*(float)HHEDGEB[j*vsize+i] +
          //            (float)VVEDGEB[j*vsize+i] *
          //            (float)VVEDGEB[j*vsize+i])));
          if (sqrt(((float)HHEDGER[j * vsize + i] *
                        (float)HHEDGER[j * vsize + i] +
                    (float)VVEDGER[j * vsize + i] *
                        (float)VVEDGER[j * vsize + i]) +
                   ((float)HHEDGEG[j * vsize + i] *
                        (float)HHEDGEG[j * vsize + i] +
                    (float)VVEDGEG[j * vsize + i] *
                        (float)VVEDGEG[j * vsize + i]) +
                   ((float)HHEDGEB[j * vsize + i] *
                        (float)HHEDGEB[j * vsize + i] +
                    (float)VVEDGEB[j * vsize + i] *
                        (float)VVEDGEB[j * vsize + i])) >=
              et * 4.0 * sqrt(6.0)) {
            *(index3 + j * vsize + i) = -1;
            EDGE_GASOSUU++;
          } else {
            *(index3 + j * vsize + i) = 1;
          }
        }
      }
    }
  }

  double *EDGERASISAY, REDGE, GEDGE, BEDGE;
  if (vvv == 4) {
    EDGERASISAY = (double *)malloc(sizeof(double) * hsize * vsize);
    // for(i=0,k=0;i<vsize;i++){
    // for(j=0;j<hsize;j++){
    // printf("%f\n",sqrt((float)HHEDGER[j*vsize+i]*(float)HHEDGER[j*vsize+i] +
    // (float)VVEDGER[j*vsize+i] * (float)VVEDGER[j*vsize+i]));
    // printf("%f\n",sqrt((float)HHEDGEG[j*vsize+i]*(float)HHEDGEG[j*vsize+i] +
    // (float)VVEDGEG[j*vsize+i] * (float)VVEDGEG[j*vsize+i]));
    // printf("%f\n",sqrt((float)HHEDGEB[j*vsize+i]*(float)HHEDGEB[j*vsize+i] +
    // (float)VVEDGEB[j*vsize+i] * (float)VVEDGEB[j*vsize+i])); if(
    // (sqrt((float)HHEDGER[j*vsize+i]*(float)HHEDGER[j*vsize+i] +
    // (float)VVEDGER[j*vsize+i] * (float)VVEDGER[j*vsize+i]) >= et*4.0) ||
    //(sqrt((float)HHEDGEG[j*vsize+i]*(float)HHEDGEG[j*vsize+i] +
    //(float)VVEDGEG[j*vsize+i] * (float)VVEDGEG[j*vsize+i]) >= et*4.0) ||
    //(sqrt((float)HHEDGEB[j*vsize+i]*(float)HHEDGEB[j*vsize+i] +
    //(float)VVEDGEB[j*vsize+i] * (float)VVEDGEB[j*vsize+i]) >= et*4.0)
    //){
    //*(index3+j*vsize+i)=-1;
    // k++;
    //} else {
    //*(index3+j*vsize+i)=1;
    //}
    //}
    //}

    //    EDGERASMAX = -9999999.9e64;
    for (int i = 0; i < vsize; i++) {
      for (int j = 0; j < hsize; j++) {
        *(index3 + j * vsize + i) = 1;
        REDGE = /*sqrt*/ (
            (float)HHEDGER[j * vsize + i] * (float)HHEDGER[j * vsize + i] +
            (float)VVEDGER[j * vsize + i] *
                (float)VVEDGER[j * vsize + i]) /*/ 4.0 / sqrt(2.0)*/;
        GEDGE = /*sqrt*/ (
            (float)HHEDGEG[j * vsize + i] * (float)HHEDGEG[j * vsize + i] +
            (float)VVEDGEG[j * vsize + i] *
                (float)VVEDGEG[j * vsize + i]) /*/ 4.0 / sqrt(2.0)*/;
        BEDGE = /*sqrt*/ (
            (float)HHEDGEB[j * vsize + i] * (float)HHEDGEB[j * vsize + i] +
            (float)VVEDGEB[j * vsize + i] *
                (float)VVEDGEB[j * vsize + i]) /*/ 4.0 / sqrt(2.0)*/;
        //            EDGERASISAY[j*vsize+i] =
        //            0.29891*REDGE+0.58661*GEDGE+0.11448*BEDGE;
        EDGERASISAY[j * vsize + i] =
            (sqrt(REDGE + GEDGE + BEDGE)) / 4.0 / sqrt(6.0);
        //            if(EDGERASISAY[j*vsize+i] > EDGERASMAX){
        //                EDGERASMAX = EDGERASISAY[j*vsize+i];
        //            }
      }
    }
    //    printf("EDGERASMAX=%f\n",EDGERASMAX);while(1);
  }

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
              if(*(index3+i) !=-1){
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
              if(*(index3+i) !=-1){
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
  U_R = 0.0;
  U_G = 0.0;
  U_B = 0.0;
  for (int i = 0; i < hsize * vsize; i++) {
    if (*(index3 + i) != -1) {
      U_R += *(YIN + i);
      U_G += *(UIN + i);
      U_B += *(VIN + i);
    }
  }
  U_R /= (double)(hsize * vsize - EDGE_GASOSUU);
  U_G /= (double)(hsize * vsize - EDGE_GASOSUU);
  U_B /= (double)(hsize * vsize - EDGE_GASOSUU);
  // U_R = (U_R>>RDIV);
  // U_G = (U_G>>GDIV);
  // U_B = (U_B>>BDIV);
  // U_R,U_G,U_Bは平均
  ct = 1000;
  eps = 1.0e-10;
  n = 3;
  A = new double *[n];
  A1 = new double *[n];
  A2 = new double *[n];
  X1 = new double *[n];
  X2 = new double *[n];
  for (i1 = 0; i1 < n; i1++) {
    A[i1] = new double[n];
    A1[i1] = new double[n];
    A2[i1] = new double[n];
    X1[i1] = new double[n];
    X2[i1] = new double[n];
  }
  // A[0][0],A[0][1],A[0][2],A[1][0],A[1][1],A[1][2],A[2][0],A[2][1],A[2][2]に分散共分散を入れる
  TMP_RR = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(RIN+i))>>RDIV)-U_R);
  TMP_GG = 0.0; //(((*(GIN+i))>>GDIV)-U_G)*(((*(GIN+i))>>GDIV)-U_G);
  TMP_BB = 0.0; //(((*(BIN+i))>>BDIV)-U_B)*(((*(BIN+i))>>BDIV)-U_B);
  TMP_RG = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(GIN+i))>>GDIV)-U_G);
  TMP_RB = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(BIN+i))>>BDIV)-U_B);
  TMP_GB = 0.0; //(((*(GIN+i))>>GDIV)-U_G)*(((*(BIN+i))>>BDIV)-U_B);
  for (int i = 0; i < hsize * vsize; i++) {
    if (*(index3 + i) != -1) {
      TMP_RR += (double)((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R));
      TMP_GG += (double)((((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G));
      TMP_BB += (double)((((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B));
      TMP_RG += (double)((((*(YIN + i))) - U_R) * (((*(UIN + i))) - U_G));
      TMP_RB += (double)((((*(YIN + i))) - U_R) * (((*(VIN + i))) - U_B));
      TMP_GB += (double)((((*(UIN + i))) - U_G) * (((*(VIN + i))) - U_B));
    }
  }
  // debug start
  // fprintf(stderr,"BUNSAN=%f %f\n",TMP_RR,TMP_RB);
  // while(1);
  // debug end
  TMP_RR /= (double)(hsize * vsize - EDGE_GASOSUU);
  TMP_GG /= (double)(hsize * vsize - EDGE_GASOSUU);
  TMP_BB /= (double)(hsize * vsize - EDGE_GASOSUU);
  TMP_RG /= (double)(hsize * vsize - EDGE_GASOSUU);
  TMP_RB /= (double)(hsize * vsize - EDGE_GASOSUU);
  TMP_GB /= (double)(hsize * vsize - EDGE_GASOSUU);
  // debug start
  // fprintf(stderr,"BUNSAN=%f %f\n",TMP_RR,TMP_RB);
  // while(1);
  // debug end
  A[0][0] = TMP_RR;
  A[0][1] = TMP_RG;
  A[1][0] = TMP_RG;
  A[1][1] = TMP_GG;
  A[2][2] = TMP_BB;
  A[0][2] = TMP_RB;
  A[2][0] = TMP_RB;
  A[1][2] = TMP_GB;
  A[2][1] = TMP_GB;

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
  double *RR;
  double *GG;
  double *BB;
  double *RRR;
  double *GGG;
  double *BBB;
  double V[3];
  double *XXX;
  double *YYY;
  double *ZZZ;
  double MAXX, MINN;
  int THRESH;
  double DTHRESH;
  RR = (double *)malloc(sizeof(double) * hsize * vsize);
  GG = (double *)malloc(sizeof(double) * hsize * vsize);
  BB = (double *)malloc(sizeof(double) * hsize * vsize);
  XXX = (double *)malloc(sizeof(double));
  YYY = (double *)malloc(sizeof(double));
  ZZZ = (double *)malloc(sizeof(double));
  RRR = (double *)malloc(sizeof(double) * hsize * vsize);
  GGG = (double *)malloc(sizeof(double) * hsize * vsize);
  BBB = (double *)malloc(sizeof(double) * hsize * vsize);
  for (int i = 0; i < hsize * vsize; i++) {
    if (*(index3 + i) != -1) {
      *(RR + i) = (double)(((*(YIN + i))) - U_R);
      *(GG + i) = (double)(((*(UIN + i))) - U_G);
      *(BB + i) = (double)(((*(VIN + i))) - U_B);
    }
  }
  V[0] = X1[0][Y];
  V[1] = X1[1][Y];
  V[2] = X1[2][Y];
  int *RRRR;
  RRRR = (int *)malloc(sizeof(int) * (hsize * vsize));
  for (int i = 0; i < hsize * vsize; i++) {
    kaiten(V, *(RR + i), *(GG + i), *(BB + i), RRR + i, GGG + i, BBB + i);
  }
  for (int i = 0, l = 0; i < hsize * vsize; i++) {
    if (*(index3 + i) != -1) {
      *(RRRR + l) = (int)(*(RRR + i) + 0.5);
      l++;
    }
  }

  double *RRRRRR, *GGGGGG, *BBBBBB;
  int GASOSUU = hsize * vsize - EDGE_GASOSUU;
  if (omh == 3 || omh == 4) {

    RRRRRR = (double *)malloc(sizeof(double) * hsize * vsize);
    GGGGGG = (double *)malloc(sizeof(double) * hsize * vsize);
    BBBBBB = (double *)malloc(sizeof(double) * hsize * vsize);
    for (int i = 0, l = 0; i < hsize * vsize; i++) {
      if (*(index3 + i) != -1) {
        *(RRRRRR + l) = *(RRR + i);
        *(GGGGGG + l) = *(GGG + i);
        *(BBBBBB + l) = *(BBB + i);
        l++;
      }
    }
    THRESH = ohtsu2(GASOSUU, RRRRRR, GGGGGG, BBBBBB, omh);
    printf("THRESH OHTSU2=%d\n", THRESH);
  }
  if (omh == 0) {
    THRESH = ohtsu(GASOSUU, RRRR);
    printf("THRESH OHTSU=%d\n", THRESH);
  } else if (omh == 1) {
    THRESH = media(GASOSUU, RRRR);
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
    *(INDEX + i) = -1;
  }
  int k2 = 0, l2 = 0;
  for (int i = 0; i < hsize * vsize; i++) {
    if (*(index3 + i) != -1) {
      // ZZ = X1[0][Y]*((double)(((*(YIN+i)))) - ((double)(U_R)+(*XXX))) +
      // X1[1][Y]*((double)(((*(UIN+i)))) - ((double)(U_G)+(*YYY))) +
      // X1[2][Y]*((double)(((*(VIN+i)))) - ((double)(U_B)+(*ZZZ))) ;
      // if(ZZ>=0.0){
      if ((double)(*(RRR + i)) >= DTHRESH) {
        *(INDEX + i) = PT[MEN][0].INDEXNO;
        k2++;
      } else {
        *(INDEX + i) = PT[MEN][1].INDEXNO;
        l2++;
      }
    }
  }
  PT[MEN][0].INDEXNUM = k2;
  PT[MEN][1].INDEXNUM = l2;
  // maxdistance tuika start
  U_R = 0.0;
  U_G = 0.0;
  U_B = 0.0;
  // 0側のmaxdistanceを求める
  for (int i = 0; i < hsize * vsize; i++) {
    if ((*(INDEX + i)) == PT[MEN][0].INDEXNO) {
      U_R += *(YIN + i);
      U_G += *(UIN + i);
      U_B += *(VIN + i);
    }
  }

  double MAXD;
  double TEMP;
  if (PT[MEN][0].INDEXNUM != 0) {
    U_R /= (double)PT[MEN][0].INDEXNUM;
    U_G /= (double)PT[MEN][0].INDEXNUM;
    U_B /= (double)PT[MEN][0].INDEXNUM;
    // U_R = (U_R>>RDIV);
    // U_G = (U_G>>GDIV);
    // U_B = (U_B>>BDIV);
    if (bun == 0) {
      MAXD = 0.0;
      for (int i = 0; i < hsize * vsize; i++) {
        if ((*(INDEX + i)) == PT[MEN][0].INDEXNO) {
          TEMP = (((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                 (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                 (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B);
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
            if ((*(INDEX + i)) == PT[MEN][0].INDEXNO) {
              TEMP += (((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                      (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                      (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B);
              // if( MAXD < TEMP){
              // MAXD = TEMP;
              // }
            }
          }
        }
        if (div == 1) {
          TEMP = 0.0;
          for (int i = 0; i < hsize * vsize; i++) {
            if ((*(INDEX + i)) == PT[MEN][0].INDEXNO) {
              TEMP += pow(((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                           (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                           (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B)),
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
            if ((*(INDEX + i)) == PT[MEN][0].INDEXNO) {
              TEMP += ((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                       (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                       (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B)) *
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
            if ((*(INDEX + i)) == PT[MEN][0].INDEXNO) {
              EDGERASISAYD = EDGERASISAY[i];
              EDGERASISAYD *= nbai;
              if (EDGERASISAYD > 255.0) {
                EDGERASISAYD = 255.0;
              }
              TEMP += pow(((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                           (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                           (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B)),
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

      // while(1);

      TMP_RR = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(RIN+i))>>RDIV)-U_R);
      TMP_GG = 0.0; //(((*(GIN+i))>>GDIV)-U_G)*(((*(GIN+i))>>GDIV)-U_G);
      TMP_BB = 0.0; //(((*(BIN+i))>>BDIV)-U_B)*(((*(BIN+i))>>BDIV)-U_B);
      TMP_RG = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(GIN+i))>>GDIV)-U_G);
      TMP_RB = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(BIN+i))>>BDIV)-U_B);
      TMP_GB = 0.0; //(((*(GIN+i))>>GDIV)-U_G)*(((*(BIN+i))>>BDIV)-U_B);
      for (int i = 0; i < hsize * vsize; i++) {
        if (*(INDEX + i) == PT[MEN][0].INDEXNO) {
          TMP_RR += (double)((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R));
          TMP_GG += (double)((((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G));
          TMP_BB += (double)((((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B));
          TMP_RG += (double)((((*(YIN + i))) - U_R) * (((*(UIN + i))) - U_G));
          TMP_RB += (double)((((*(YIN + i))) - U_R) * (((*(VIN + i))) - U_B));
          TMP_GB += (double)((((*(UIN + i))) - U_G) * (((*(VIN + i))) - U_B));
        }
      }
      // debug start
      // fprintf(stderr,"BUNSAN=%f %f\n",TMP_RR,TMP_RB);
      // while(1);
      // debug end
      TMP_RR /= (double)(PT[MEN][0].INDEXNUM);
      TMP_GG /= (double)(PT[MEN][0].INDEXNUM);
      TMP_BB /= (double)(PT[MEN][0].INDEXNUM);
      TMP_RG /= (double)(PT[MEN][0].INDEXNUM);
      TMP_RB /= (double)(PT[MEN][0].INDEXNUM);
      TMP_GB /= (double)(PT[MEN][0].INDEXNUM);
      // debug start
      // fprintf(stderr,"BUNSAN=%f %f\n",TMP_RR,TMP_RB);
      // while(1);
      // debug end
      A[0][0] = TMP_RR;
      A[0][1] = TMP_RG;
      A[1][0] = TMP_RG;
      A[1][1] = TMP_GG;
      A[2][2] = TMP_BB;
      A[0][2] = TMP_RB;
      A[2][0] = TMP_RB;
      A[1][2] = TMP_GB;
      A[2][1] = TMP_GB;

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
          if (*(INDEX + i) == PT[MEN][0].INDEXNO) {
            *(RR + i) = (double)(((*(YIN + i))) - U_R);
            *(GG + i) = (double)(((*(UIN + i))) - U_G);
            *(BB + i) = (double)(((*(VIN + i))) - U_B);
          }
        }
        V[0] = X1[0][Y];
        V[1] = X1[1][Y];
        V[2] = X1[2][Y];
        for (int i = 0; i < hsize * vsize; i++) {
          kaiten(V, *(RR + i), *(GG + i), *(BB + i), RRR + i, GGG + i, BBB + i);
        }
        // for(i=0,l=0;i<hsize*vsize;i++){
        // if(*(index3+i)!=-1){
        //*(RRRR+l) = (int)(*(RRR+i)+0.5);
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
          if (*(INDEX + i) == PT[MEN][0].INDEXNO) {
            *(RRRRRR + l2) = *(RRR + i);
            *(GGGGGG + l2) = *(GGG + i);
            *(BBBBBB + l2) = *(BBB + i);
            l2++;
          }
        }
        MAXX = -99.9e64;
        MINN = 99.9e64;
        for (int i = 0; i < l2; i++) {
          if (MAXX < *(RRRRRR + i)) {
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

  U_R = 0.0;
  U_G = 0.0;
  U_B = 0.0;
  // 1側のmaxdistanceを求める
  for (int i = 0; i < hsize * vsize; i++) {
    if ((*(INDEX + i)) == PT[MEN][1].INDEXNO) {
      U_R += *(YIN + i);
      U_G += *(UIN + i);
      U_B += *(VIN + i);
    }
  }

  if (PT[MEN][1].INDEXNUM != 0) {
    U_R /= (double)PT[MEN][1].INDEXNUM;
    U_G /= (double)PT[MEN][1].INDEXNUM;
    U_B /= (double)PT[MEN][1].INDEXNUM;
    // U_R = (U_R>>RDIV);
    // U_G = (U_G>>GDIV);
    // U_B = (U_B>>BDIV);
    if (bun == 0) {
      MAXD = 0.0;
      for (int i = 0; i < hsize * vsize; i++) {
        if ((*(INDEX + i)) == PT[MEN][1].INDEXNO) {
          TEMP = ((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                  (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                  (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B));
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
            if ((*(INDEX + i)) == PT[MEN][1].INDEXNO) {
              TEMP += ((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                       (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                       (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B));
              // if( MAXD < TEMP){
              // MAXD = TEMP;
              // }
            }
          }
        }
        if (div == 1) {
          TEMP = 0.0;
          for (int i = 0; i < hsize * vsize; i++) {
            if ((*(INDEX + i)) == PT[MEN][1].INDEXNO) {
              TEMP += pow(((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                           (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                           (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B)),
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
            if ((*(INDEX + i)) == PT[MEN][1].INDEXNO) {
              TEMP += ((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                       (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                       (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B)) *
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
            if ((*(INDEX + i)) == PT[MEN][1].INDEXNO) {
              EDGERASISAYD = EDGERASISAY[i];
              EDGERASISAYD *= nbai;
              if (EDGERASISAYD > 255.0) {
                EDGERASISAYD = 255.0;
              }

              TEMP += pow(((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                           (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                           (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B)),
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

      TMP_RR = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(RIN+i))>>RDIV)-U_R);
      TMP_GG = 0.0; //(((*(GIN+i))>>GDIV)-U_G)*(((*(GIN+i))>>GDIV)-U_G);
      TMP_BB = 0.0; //(((*(BIN+i))>>BDIV)-U_B)*(((*(BIN+i))>>BDIV)-U_B);
      TMP_RG = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(GIN+i))>>GDIV)-U_G);
      TMP_RB = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(BIN+i))>>BDIV)-U_B);
      TMP_GB = 0.0; //(((*(GIN+i))>>GDIV)-U_G)*(((*(BIN+i))>>BDIV)-U_B);
      for (int i = 0; i < hsize * vsize; i++) {
        if (*(INDEX + i) == PT[MEN][1].INDEXNO) {
          TMP_RR += (double)((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R));
          TMP_GG += (double)((((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G));
          TMP_BB += (double)((((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B));
          TMP_RG += (double)((((*(YIN + i))) - U_R) * (((*(UIN + i))) - U_G));
          TMP_RB += (double)((((*(YIN + i))) - U_R) * (((*(VIN + i))) - U_B));
          TMP_GB += (double)((((*(UIN + i))) - U_G) * (((*(VIN + i))) - U_B));
        }
      }
      // debug start
      // fprintf(stderr,"BUNSAN=%f %f\n",TMP_RR,TMP_RB);
      // while(1);
      // debug end
      TMP_RR /= (double)(PT[MEN][1].INDEXNUM);
      TMP_GG /= (double)(PT[MEN][1].INDEXNUM);
      TMP_BB /= (double)(PT[MEN][1].INDEXNUM);
      TMP_RG /= (double)(PT[MEN][1].INDEXNUM);
      TMP_RB /= (double)(PT[MEN][1].INDEXNUM);
      TMP_GB /= (double)(PT[MEN][1].INDEXNUM);
      // debug start
      // fprintf(stderr,"BUNSAN=%f %f\n",TMP_RR,TMP_RB);
      // while(1);
      // debug end
      A[0][0] = TMP_RR;
      A[0][1] = TMP_RG;
      A[1][0] = TMP_RG;
      A[1][1] = TMP_GG;
      A[2][2] = TMP_BB;
      A[0][2] = TMP_RB;
      A[2][0] = TMP_RB;
      A[1][2] = TMP_GB;
      A[2][1] = TMP_GB;

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
          if (*(INDEX + i) == PT[MEN][1].INDEXNO) {
            *(RR + i) = (double)(((*(YIN + i))) - U_R);
            *(GG + i) = (double)(((*(UIN + i))) - U_G);
            *(BB + i) = (double)(((*(VIN + i))) - U_B);
          }
        }
        V[0] = X1[0][Y];
        V[1] = X1[1][Y];
        V[2] = X1[2][Y];
        for (int i = 0; i < hsize * vsize; i++) {
          kaiten(V, *(RR + i), *(GG + i), *(BB + i), RRR + i, GGG + i, BBB + i);
        }
        // for(i=0,l=0;i<hsize*vsize;i++){
        // if(*(index3+i)!=-1){
        //*(RRRR+l) = (int)(*(RRR+i)+0.5);
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
          if (*(INDEX + i) == PT[MEN][1].INDEXNO) {
            *(RRRRRR + l2) = *(RRR + i);
            *(GGGGGG + l2) = *(GGG + i);
            *(BBBBBB + l2) = *(BBB + i);
            l2++;
          }
        }
        MAXX = -99.9e64;
        MINN = 99.9e64;
        for (int i = 0; i < l2; i++) {
          if (MAXX < *(RRRRRR + i)) {
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

  U_R = 0.0;
  U_G = 0.0;
  U_B = 0.0;
  for (int i = 0; i < hsize * vsize; i++) {
    if ((*(INDEX + i)) == j2) {
      U_R += *(YIN + i);
      U_G += *(UIN + i);
      U_B += *(VIN + i);
    }
  }

  U_R /= (double)PT[MEN][j2].INDEXNUM;
  U_G /= (double)PT[MEN][j2].INDEXNUM;
  U_B /= (double)PT[MEN][j2].INDEXNUM;
  // U_R = (U_R>>RDIV);
  // U_G = (U_G>>GDIV);
  // U_B = (U_B>>BDIV);
  // debug start
  // fprintf(stderr,"2kaimeheikin %d %d %d\n",U_R,U_G,U_B);
  // debug end
  TMP_RR = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(RIN+i))>>RDIV)-U_R);
  TMP_GG = 0.0; //(((*(GIN+i))>>GDIV)-U_G)*(((*(GIN+i))>>GDIV)-U_G);
  TMP_BB = 0.0; //(((*(BIN+i))>>BDIV)-U_B)*(((*(BIN+i))>>BDIV)-U_B);
  TMP_RG = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(GIN+i))>>GDIV)-U_G);
  TMP_RB = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(BIN+i))>>BDIV)-U_B);
  TMP_GB = 0.0; //(((*(GIN+i))>>GDIV)-U_G)*(((*(BIN+i))>>BDIV)-U_B);
  for (int i = 0; i < hsize * vsize; i++) {
    if ((*(INDEX + i)) == j2) {
      TMP_RR += (double)((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R));
      TMP_GG += (double)((((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G));
      TMP_BB += (double)((((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B));
      TMP_RG += (double)((((*(YIN + i))) - U_R) * (((*(UIN + i))) - U_G));
      TMP_RB += (double)((((*(YIN + i))) - U_R) * (((*(VIN + i))) - U_B));
      TMP_GB += (double)((((*(UIN + i))) - U_G) * (((*(VIN + i))) - U_B));
    }
  }
  TMP_RR /= (double)(PT[MEN][j2].INDEXNUM);
  TMP_GG /= (double)(PT[MEN][j2].INDEXNUM);
  TMP_BB /= (double)(PT[MEN][j2].INDEXNUM);
  TMP_RG /= (double)(PT[MEN][j2].INDEXNUM);
  TMP_RB /= (double)(PT[MEN][j2].INDEXNUM);
  TMP_GB /= (double)(PT[MEN][j2].INDEXNUM);

  // 2番目の分轄の軸を求める（固有ベクトルの計算）
  A[0][0] = TMP_RR;
  A[0][1] = TMP_RG;
  A[1][0] = TMP_RG;
  A[1][1] = TMP_GG;
  A[2][2] = TMP_BB;
  A[0][2] = TMP_RB;
  A[2][0] = TMP_RB;
  A[1][2] = TMP_GB;
  A[2][1] = TMP_GB;

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
    if ((*(INDEX + i)) == j2) {
      *(RR + i) = (double)(((*(YIN + i))) - U_R);
      *(GG + i) = (double)(((*(UIN + i))) - U_G);
      *(BB + i) = (double)(((*(VIN + i))) - U_B);
    }
  }
  V[0] = X1[0][Y];
  V[1] = X1[1][Y];
  V[2] = X1[2][Y];
  for (int i = 0; i < hsize * vsize; i++) {
    if ((*(INDEX + i)) == j2) {
      kaiten(V, *(RR + i), *(GG + i), *(BB + i), RRR + i, GGG + i, BBB + i);
    }
  }
  GASOSUU = PT[MEN][j2].INDEXNUM;
  for (int i = 0; i < hsize * vsize; i++) {
    *(RRR33 + i) = (double)(*(RRR + i));
  }
  for (int i = 0, k = 0; i < hsize * vsize; i++) {
    if ((*(INDEX + i)) == j2) {
      *(RR + k) = (double)(*(RRR + i));
      *(GG + k) = (double)(*(GGG + i));
      *(BB + k) = (double)(*(BBB + i));
      k++;
    }
  }
  for (int i = 0; i < hsize * vsize; i++) {
    *(RRRR + i) = (int)(*(RR + i) + 0.5);
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
    THRESH = ohtsu2(GASOSUU, RR, GG, BB, omh);
  }
  if (omh == 0) {
    THRESH = ohtsu(GASOSUU, RRRR);
  } else if (omh == 1) {
    THRESH = media(GASOSUU, RRRR);
  } else if (omh == 2) {
    THRESH = 0;
  }
  DTHRESH = (double)(THRESH);
  modoshi(V, DTHRESH, 0.0, 0.0, XXX, YYY, ZZZ);

  int k3 = 0, l3 = 0;
  for (int i = 0; i < hsize * vsize; i++) {
    if ((*(INDEX + i)) == j2) {
      // ZZ = X1[0][Y]*((double)(((*(YIN+i)))) - ((double)(U_R)+(*XXX))) +
      // X1[1][Y]*((double)(((*(UIN+i)))) - ((double)(U_G)+(*YYY))) +
      // X1[2][Y]*((double)(((*(VIN+i)))) - ((double)(U_B)+(*ZZZ))) ;
      // if(ZZ>=0.0){
      if (*(RRR33 + i) >= DTHRESH) {
        *(INDEX + i) = DIVIDENUM - 1;
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
  U_R = 0.0;
  U_G = 0.0;
  U_B = 0.0;
  // 0側のmaxdistanceを求める
  for (int i = 0; i < hsize * vsize; i++) {
    if ((*(INDEX + i)) == PT[MEN][0].INDEXNO) {
      U_R += *(YIN + i);
      U_G += *(UIN + i);
      U_B += *(VIN + i);
    }
  }
  if (PT[MEN][0].INDEXNUM != 0) {
    U_R /= (double)PT[MEN][0].INDEXNUM;
    U_G /= (double)PT[MEN][0].INDEXNUM;
    U_B /= (double)PT[MEN][0].INDEXNUM;
    // U_R = (U_R>>RDIV);
    // U_G = (U_G>>GDIV);
    // U_B = (U_B>>BDIV);
    if (bun == 0) {
      MAXD = 0.0;
      for (int i = 0; i < hsize * vsize; i++) {
        if ((*(INDEX + i)) == PT[MEN][0].INDEXNO) {
          TEMP = (((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                 (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                 (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B);
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
            if ((*(INDEX + i)) == PT[MEN][0].INDEXNO) {
              TEMP += (((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                      (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                      (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B);
              // if( MAXD < TEMP){
              // MAXD = TEMP;
              // }
            }
          }
        }
        if (div == 1) {
          TEMP = 0.0;
          for (int i = 0; i < hsize * vsize; i++) {
            if ((*(INDEX + i)) == PT[MEN][0].INDEXNO) {
              TEMP += pow(((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                           (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                           (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B)),
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
            if ((*(INDEX + i)) == PT[MEN][0].INDEXNO) {
              TEMP += ((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                       (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                       (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B)) *
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
            if ((*(INDEX + i)) == PT[MEN][0].INDEXNO) {
              EDGERASISAYD = EDGERASISAY[i];
              EDGERASISAYD *= nbai;
              if (EDGERASISAYD > 255.0) {
                EDGERASISAYD = 255.0;
              }

              TEMP += pow(((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                           (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                           (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B)),
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

      TMP_RR = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(RIN+i))>>RDIV)-U_R);
      TMP_GG = 0.0; //(((*(GIN+i))>>GDIV)-U_G)*(((*(GIN+i))>>GDIV)-U_G);
      TMP_BB = 0.0; //(((*(BIN+i))>>BDIV)-U_B)*(((*(BIN+i))>>BDIV)-U_B);
      TMP_RG = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(GIN+i))>>GDIV)-U_G);
      TMP_RB = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(BIN+i))>>BDIV)-U_B);
      TMP_GB = 0.0; //(((*(GIN+i))>>GDIV)-U_G)*(((*(BIN+i))>>BDIV)-U_B);
      for (int i = 0; i < hsize * vsize; i++) {
        if ((*(INDEX + i)) == PT[MEN][0].INDEXNO) {
          TMP_RR += (double)((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R));
          TMP_GG += (double)((((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G));
          TMP_BB += (double)((((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B));
          TMP_RG += (double)((((*(YIN + i))) - U_R) * (((*(UIN + i))) - U_G));
          TMP_RB += (double)((((*(YIN + i))) - U_R) * (((*(VIN + i))) - U_B));
          TMP_GB += (double)((((*(UIN + i))) - U_G) * (((*(VIN + i))) - U_B));
        }
      }
      TMP_RR /= (double)(PT[MEN][0].INDEXNUM);
      TMP_GG /= (double)(PT[MEN][0].INDEXNUM);
      TMP_BB /= (double)(PT[MEN][0].INDEXNUM);
      TMP_RG /= (double)(PT[MEN][0].INDEXNUM);
      TMP_RB /= (double)(PT[MEN][0].INDEXNUM);
      TMP_GB /= (double)(PT[MEN][0].INDEXNUM);

      // 2番目の分轄の軸を求める（固有ベクトルの計算）
      A[0][0] = TMP_RR;
      A[0][1] = TMP_RG;
      A[1][0] = TMP_RG;
      A[1][1] = TMP_GG;
      A[2][2] = TMP_BB;
      A[0][2] = TMP_RB;
      A[2][0] = TMP_RB;
      A[1][2] = TMP_GB;
      A[2][1] = TMP_GB;

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
          if (*(INDEX + i) == PT[MEN][0].INDEXNO) {
            *(RR + i) = (double)(((*(YIN + i))) - U_R);
            *(GG + i) = (double)(((*(UIN + i))) - U_G);
            *(BB + i) = (double)(((*(VIN + i))) - U_B);
          }
        }
        V[0] = X1[0][Y];
        V[1] = X1[1][Y];
        V[2] = X1[2][Y];
        for (int i = 0; i < hsize * vsize; i++) {
          kaiten(V, *(RR + i), *(GG + i), *(BB + i), RRR + i, GGG + i, BBB + i);
        }
        // for(i=0,l=0;i<hsize*vsize;i++){
        // if(*(index3+i)!=-1){
        //*(RRRR+l) = (int)(*(RRR+i)+0.5);
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
          if (*(INDEX + i) == PT[MEN][0].INDEXNO) {
            *(RRRRRR + l2) = *(RRR + i);
            *(GGGGGG + l2) = *(GGG + i);
            *(BBBBBB + l2) = *(BBB + i);
            l2++;
          }
        }
        MAXX = -99.9e64;
        MINN = 99.9e64;
        for (int i = 0; i < l2; i++) {
          if (MAXX < *(RRRRRR + i)) {
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
  U_R = 0.0;
  U_G = 0.0;
  U_B = 0.0;

  // 1側のmaxdistanceを求める
  for (int i = 0; i < hsize * vsize; i++) {
    if ((*(INDEX + i)) == PT[MEN][1].INDEXNO) {
      U_R += *(YIN + i);
      U_G += *(UIN + i);
      U_B += *(VIN + i);
    }
  }

  if (PT[MEN][1].INDEXNUM != 0) {
    U_R /= (double)PT[MEN][1].INDEXNUM;
    U_G /= (double)PT[MEN][1].INDEXNUM;
    U_B /= (double)PT[MEN][1].INDEXNUM;
    // U_R = (U_R>>RDIV);
    // U_G = (U_G>>GDIV);
    // U_B = (U_B>>BDIV);
    if (bun == 0) {
      MAXD = 0.0;
      for (int i = 0; i < hsize * vsize; i++) {
        if ((*(INDEX + i)) == PT[MEN][1].INDEXNO) {
          TEMP = (((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                 (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                 (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B);
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
            if ((*(INDEX + i)) == PT[MEN][1].INDEXNO) {
              TEMP += (((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                      (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                      (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B);
              // if( MAXD < TEMP){
              // MAXD = TEMP;
              // }
            }
          }
        }
        if (div == 1) {
          TEMP = 0.0;
          for (int i = 0; i < hsize * vsize; i++) {
            if ((*(INDEX + i)) == PT[MEN][1].INDEXNO) {
              TEMP += pow(((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                           (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                           (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B)),
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
            if ((*(INDEX + i)) == PT[MEN][1].INDEXNO) {
              TEMP += ((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                       (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                       (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B)) *
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
            if ((*(INDEX + i)) == PT[MEN][1].INDEXNO) {
              EDGERASISAYD = EDGERASISAY[i];
              EDGERASISAYD *= nbai;
              if (EDGERASISAYD > 255.0) {
                EDGERASISAYD = 255.0;
              }

              TEMP += pow(((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                           (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                           (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B)),
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

      TMP_RR = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(RIN+i))>>RDIV)-U_R);
      TMP_GG = 0.0; //(((*(GIN+i))>>GDIV)-U_G)*(((*(GIN+i))>>GDIV)-U_G);
      TMP_BB = 0.0; //(((*(BIN+i))>>BDIV)-U_B)*(((*(BIN+i))>>BDIV)-U_B);
      TMP_RG = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(GIN+i))>>GDIV)-U_G);
      TMP_RB = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(BIN+i))>>BDIV)-U_B);
      TMP_GB = 0.0; //(((*(GIN+i))>>GDIV)-U_G)*(((*(BIN+i))>>BDIV)-U_B);
      for (int i = 0; i < hsize * vsize; i++) {
        if ((*(INDEX + i)) == PT[MEN][1].INDEXNO) {
          TMP_RR += (double)((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R));
          TMP_GG += (double)((((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G));
          TMP_BB += (double)((((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B));
          TMP_RG += (double)((((*(YIN + i))) - U_R) * (((*(UIN + i))) - U_G));
          TMP_RB += (double)((((*(YIN + i))) - U_R) * (((*(VIN + i))) - U_B));
          TMP_GB += (double)((((*(UIN + i))) - U_G) * (((*(VIN + i))) - U_B));
        }
      }
      TMP_RR /= (double)(PT[MEN][1].INDEXNUM);
      TMP_GG /= (double)(PT[MEN][1].INDEXNUM);
      TMP_BB /= (double)(PT[MEN][1].INDEXNUM);
      TMP_RG /= (double)(PT[MEN][1].INDEXNUM);
      TMP_RB /= (double)(PT[MEN][1].INDEXNUM);
      TMP_GB /= (double)(PT[MEN][1].INDEXNUM);

      // 2番目の分轄の軸を求める（固有ベクトルの計算）
      A[0][0] = TMP_RR;
      A[0][1] = TMP_RG;
      A[1][0] = TMP_RG;
      A[1][1] = TMP_GG;
      A[2][2] = TMP_BB;
      A[0][2] = TMP_RB;
      A[2][0] = TMP_RB;
      A[1][2] = TMP_GB;
      A[2][1] = TMP_GB;

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
          if (*(INDEX + i) == PT[MEN][1].INDEXNO) {
            *(RR + i) = (double)(((*(YIN + i))) - U_R);
            *(GG + i) = (double)(((*(UIN + i))) - U_G);
            *(BB + i) = (double)(((*(VIN + i))) - U_B);
          }
        }
        V[0] = X1[0][Y];
        V[1] = X1[1][Y];
        V[2] = X1[2][Y];
        for (int i = 0; i < hsize * vsize; i++) {
          kaiten(V, *(RR + i), *(GG + i), *(BB + i), RRR + i, GGG + i, BBB + i);
        }
        // for(i=0,l=0;i<hsize*vsize;i++){
        // if(*(index3+i)!=-1){
        //*(RRRR+l) = (int)(*(RRR+i)+0.5);
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
          if (*(INDEX + i) == PT[MEN][1].INDEXNO) {
            *(RRRRRR + l2) = *(RRR + i);
            *(GGGGGG + l2) = *(GGG + i);
            *(BBBBBB + l2) = *(BBB + i);
            l2++;
          }
        }
        MAXX = -99.9e64;
        MINN = 99.9e64;
        for (int i = 0; i < l2; i++) {
          if (MAXX < *(RRRRRR + i)) {
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
    U_R = 0.0;
    U_G = 0.0;
    U_B = 0.0;
    for (int i = 0; i < hsize * vsize; i++) {
      if ((*(INDEX + i)) == PT[MEN][NUM].INDEXNO) {
        U_R += *(YIN + i);
        U_G += *(UIN + i);
        U_B += *(VIN + i);
      }
    }
    U_R /= (double)PT[MEN][NUM].INDEXNUM;
    U_G /= (double)PT[MEN][NUM].INDEXNUM;
    U_B /= (double)PT[MEN][NUM].INDEXNUM;
    // U_R = (U_R>>RDIV);
    // U_G = (U_G>>GDIV);
    // U_B = (U_B>>BDIV);
    // debug start
    // if(DIVIDENUM == 41){
    // fprintf(stderr,"AVE= %d %d %d\n",U_R,U_G,U_B);
    // fprintf(stderr,"NUM= %d\n",NUM);
    // }
    // debug end
    TMP_RR = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(RIN+i))>>RDIV)-U_R);
    TMP_GG = 0.0; //(((*(GIN+i))>>GDIV)-U_G)*(((*(GIN+i))>>GDIV)-U_G);
    TMP_BB = 0.0; //(((*(BIN+i))>>BDIV)-U_B)*(((*(BIN+i))>>BDIV)-U_B);
    TMP_RG = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(GIN+i))>>GDIV)-U_G);
    TMP_RB = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(BIN+i))>>BDIV)-U_B);
    TMP_GB = 0.0; //(((*(GIN+i))>>GDIV)-U_G)*(((*(BIN+i))>>BDIV)-U_B);
    for (int i = 0; i < hsize * vsize; i++) {
      if ((*(INDEX + i)) == PT[MEN][NUM].INDEXNO) {
        TMP_RR += (double)((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R));
        TMP_GG += (double)((((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G));
        TMP_BB += (double)((((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B));
        TMP_RG += (double)((((*(YIN + i))) - U_R) * (((*(UIN + i))) - U_G));
        TMP_RB += (double)((((*(YIN + i))) - U_R) * (((*(VIN + i))) - U_B));
        TMP_GB += (double)((((*(UIN + i))) - U_G) * (((*(VIN + i))) - U_B));
      }
    }
    TMP_RR /= (double)(PT[MEN][NUM].INDEXNUM);
    TMP_GG /= (double)(PT[MEN][NUM].INDEXNUM);
    TMP_BB /= (double)(PT[MEN][NUM].INDEXNUM);
    TMP_RG /= (double)(PT[MEN][NUM].INDEXNUM);
    TMP_RB /= (double)(PT[MEN][NUM].INDEXNUM);
    TMP_GB /= (double)(PT[MEN][NUM].INDEXNUM);

    // DIVIDENUM-1番目の分轄の軸を求める（固有ベクトルの計算）
    A[0][0] = TMP_RR;
    A[0][1] = TMP_RG;
    A[1][0] = TMP_RG;
    A[1][1] = TMP_GG;
    A[2][2] = TMP_BB;
    A[0][2] = TMP_RB;
    A[2][0] = TMP_RB;
    A[1][2] = TMP_GB;
    A[2][1] = TMP_GB;
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
      if ((*(INDEX + i)) == PT[MEN][NUM].INDEXNO) {
        *(RR + i) = (double)(((*(YIN + i))) - U_R);
        *(GG + i) = (double)(((*(UIN + i))) - U_G);
        *(BB + i) = (double)(((*(VIN + i))) - U_B);
      }
    }
    V[0] = X1[0][Y];
    V[1] = X1[1][Y];
    V[2] = X1[2][Y];
    for (int i = 0; i < hsize * vsize; i++) {
      if ((*(INDEX + i)) == PT[MEN][NUM].INDEXNO) {
        kaiten(V, *(RR + i), *(GG + i), *(BB + i), RRR + i, GGG + i, BBB + i);
      }
    }
    GASOSUU = PT[MEN][NUM].INDEXNUM;
    for (int i = 0; i < hsize * vsize; i++) {
      *(RRR33 + i) = (double)(*(RRR + i));
    }
    for (int i = 0, k = 0; i < hsize * vsize; i++) {
      if ((*(INDEX + i)) == PT[MEN][NUM].INDEXNO) {
        *(RR + k) = (double)(*(RRR + i));
        *(GG + k) = (double)(*(GGG + i));
        *(BB + k) = (double)(*(BBB + i));
        k++;
      }
    }
    for (int i = 0; i < hsize * vsize; i++) {
      *(RRRR + i) = (int)(*(RR + i) + 0.5);
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
    if (omh == 3 || omh == 4) {
      THRESH = ohtsu2(GASOSUU, RR, GG, BB, omh);
    }
    if (omh == 0) {
      THRESH = ohtsu(GASOSUU, RRRR);
    } else if (omh == 1) {
      THRESH = media(GASOSUU, RRRR);
    } else if (omh == 2) {
      THRESH = 0;
    }
    DTHRESH = (double)(THRESH);
    modoshi(V, DTHRESH, 0.0, 0.0, XXX, YYY, ZZZ);

    int k2 = 0, l2 = 0;
    for (int i = 0; i < hsize * vsize; i++) {
      if ((*(INDEX + i)) == PT[MEN][NUM].INDEXNO) {
        // ZZ = X1[0][Y]*((double)(((*(YIN+i)))) - ((double)(U_R)+(*XXX))) +
        // X1[1][Y]*((double)(((*(UIN+i)))) - ((double)(U_G)+(*YYY))) +
        // X1[2][Y]*((double)(((*(VIN+i)))) - ((double)(U_B)+(*ZZZ))) ;
        // if(ZZ>=0.0){
        if (*(RRR33 + i) >= DTHRESH) {
          *(INDEX + i) = DIVIDENUM - 1;
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
    U_R = 0.0;
    U_G = 0.0;
    U_B = 0.0;
    // 0側のmaxdistanceを求める
    for (int i = 0; i < hsize * vsize; i++) {
      if ((*(INDEX + i)) == PT[MEN][0].INDEXNO) {
        U_R += *(YIN + i);
        U_G += *(UIN + i);
        U_B += *(VIN + i);
      }
    }

    if (PT[MEN][0].INDEXNUM != 0) {
      U_R /= PT[MEN][0].INDEXNUM;
      U_G /= PT[MEN][0].INDEXNUM;
      U_B /= PT[MEN][0].INDEXNUM;
      // U_R = (U_R>>RDIV);
      // U_G = (U_G>>GDIV);
      // U_B = (U_B>>BDIV);
      if (bun == 0) {
        MAXD = 0.0;
        for (int i = 0; i < hsize * vsize; i++) {
          if ((*(INDEX + i)) == PT[MEN][0].INDEXNO) {
            TEMP = (((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                   (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                   (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B);
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
              if ((*(INDEX + i)) == PT[MEN][0].INDEXNO) {
                TEMP += (((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                        (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                        (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B);
                // if( MAXD < TEMP){
                // MAXD = TEMP;
                // }
              }
            }
          }
          if (div == 1) {
            TEMP = 0.0;
            for (int i = 0; i < hsize * vsize; i++) {
              if ((*(INDEX + i)) == PT[MEN][0].INDEXNO) {
                TEMP += pow(((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                             (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                             (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B)),
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
              if ((*(INDEX + i)) == PT[MEN][0].INDEXNO) {
                TEMP += ((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                         (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                         (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B)) *
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
              if ((*(INDEX + i)) == PT[MEN][0].INDEXNO) {
                EDGERASISAYD = EDGERASISAY[i];
                EDGERASISAYD *= nbai;
                if (EDGERASISAYD > 255.0) {
                  EDGERASISAYD = 255.0;
                }

                TEMP += pow(((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                             (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                             (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B)),
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

        TMP_RR = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(RIN+i))>>RDIV)-U_R);
        TMP_GG = 0.0; //(((*(GIN+i))>>GDIV)-U_G)*(((*(GIN+i))>>GDIV)-U_G);
        TMP_BB = 0.0; //(((*(BIN+i))>>BDIV)-U_B)*(((*(BIN+i))>>BDIV)-U_B);
        TMP_RG = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(GIN+i))>>GDIV)-U_G);
        TMP_RB = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(BIN+i))>>BDIV)-U_B);
        TMP_GB = 0.0; //(((*(GIN+i))>>GDIV)-U_G)*(((*(BIN+i))>>BDIV)-U_B);
        for (int i = 0; i < hsize * vsize; i++) {
          if ((*(INDEX + i)) == PT[MEN][0].INDEXNO) {
            TMP_RR += (double)((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R));
            TMP_GG += (double)((((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G));
            TMP_BB += (double)((((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B));
            TMP_RG += (double)((((*(YIN + i))) - U_R) * (((*(UIN + i))) - U_G));
            TMP_RB += (double)((((*(YIN + i))) - U_R) * (((*(VIN + i))) - U_B));
            TMP_GB += (double)((((*(UIN + i))) - U_G) * (((*(VIN + i))) - U_B));
          }
        }
        TMP_RR /= (double)(PT[MEN][0].INDEXNUM);
        TMP_GG /= (double)(PT[MEN][0].INDEXNUM);
        TMP_BB /= (double)(PT[MEN][0].INDEXNUM);
        TMP_RG /= (double)(PT[MEN][0].INDEXNUM);
        TMP_RB /= (double)(PT[MEN][0].INDEXNUM);
        TMP_GB /= (double)(PT[MEN][0].INDEXNUM);

        // DIVIDENUM-1番目の分轄の軸を求める（固有ベクトルの計算）
        A[0][0] = TMP_RR;
        A[0][1] = TMP_RG;
        A[1][0] = TMP_RG;
        A[1][1] = TMP_GG;
        A[2][2] = TMP_BB;
        A[0][2] = TMP_RB;
        A[2][0] = TMP_RB;
        A[1][2] = TMP_GB;
        A[2][1] = TMP_GB;
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
            if (*(INDEX + i) == PT[MEN][0].INDEXNO) {
              *(RR + i) = (double)(((*(YIN + i))) - U_R);
              *(GG + i) = (double)(((*(UIN + i))) - U_G);
              *(BB + i) = (double)(((*(VIN + i))) - U_B);
            }
          }
          V[0] = X1[0][Y];
          V[1] = X1[1][Y];
          V[2] = X1[2][Y];
          for (int i = 0; i < hsize * vsize; i++) {
            kaiten(V, *(RR + i), *(GG + i), *(BB + i), RRR + i, GGG + i,
                   BBB + i);
          }
          // for(i=0,l=0;i<hsize*vsize;i++){
          // if(*(index3+i)!=-1){
          //*(RRRR+l) = (int)(*(RRR+i)+0.5);
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
            if (*(INDEX + i) == PT[MEN][0].INDEXNO) {
              *(RRRRRR + l2) = *(RRR + i);
              *(GGGGGG + l2) = *(GGG + i);
              *(BBBBBB + l2) = *(BBB + i);
              l2++;
            }
          }
          MAXX = -99.9e64;
          MINN = 99.9e64;
          for (int i = 0; i < l2; i++) {
            if (MAXX < *(RRRRRR + i)) {
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
    } else {
      PT[MEN][0].MAXDISTANCE = 0.0;
    }

    U_R = 0.0;
    U_G = 0.0;
    U_B = 0.0;
    // 1側のmaxdistanceを求める
    for (int i = 0; i < hsize * vsize; i++) {
      if ((*(INDEX + i)) == PT[MEN][1].INDEXNO) {
        U_R += *(YIN + i);
        U_G += *(UIN + i);
        U_B += *(VIN + i);
      }
    }
    if (PT[MEN][1].INDEXNUM != 0) {
      U_R /= PT[MEN][1].INDEXNUM;
      U_G /= PT[MEN][1].INDEXNUM;
      U_B /= PT[MEN][1].INDEXNUM;
      // U_R = (U_R>>RDIV);
      // U_G = (U_G>>GDIV);
      // U_B = (U_B>>BDIV);
      if (bun == 0) {
        MAXD = 0.0;
        for (int i = 0; i < hsize * vsize; i++) {
          if ((*(INDEX + i)) == PT[MEN][1].INDEXNO) {
            TEMP = (((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                   (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                   (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B);
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
              if ((*(INDEX + i)) == PT[MEN][1].INDEXNO) {
                TEMP += (((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                        (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                        (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B);
                // if( MAXD < TEMP){
                // MAXD = TEMP;
                // }
              }
            }
          }
          if (div == 1) {
            TEMP = 0.0;
            for (int i = 0; i < hsize * vsize; i++) {
              if ((*(INDEX + i)) == PT[MEN][1].INDEXNO) {
                TEMP += pow(((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                             (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                             (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B)),
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
              if ((*(INDEX + i)) == PT[MEN][1].INDEXNO) {
                TEMP += ((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                         (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                         (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B)) *
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
              if ((*(INDEX + i)) == PT[MEN][1].INDEXNO) {
                EDGERASISAYD = EDGERASISAY[i];
                EDGERASISAYD *= nbai;
                if (EDGERASISAYD > 255.0) {
                  EDGERASISAYD = 255.0;
                }

                TEMP += pow(((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R) +
                             (((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G) +
                             (((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B)),
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

        TMP_RR = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(RIN+i))>>RDIV)-U_R);
        TMP_GG = 0.0; //(((*(GIN+i))>>GDIV)-U_G)*(((*(GIN+i))>>GDIV)-U_G);
        TMP_BB = 0.0; //(((*(BIN+i))>>BDIV)-U_B)*(((*(BIN+i))>>BDIV)-U_B);
        TMP_RG = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(GIN+i))>>GDIV)-U_G);
        TMP_RB = 0.0; //(((*(RIN+i))>>RDIV)-U_R)*(((*(BIN+i))>>BDIV)-U_B);
        TMP_GB = 0.0; //(((*(GIN+i))>>GDIV)-U_G)*(((*(BIN+i))>>BDIV)-U_B);
        for (int i = 0; i < hsize * vsize; i++) {
          if ((*(INDEX + i)) == PT[MEN][1].INDEXNO) {
            TMP_RR += (double)((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R));
            TMP_GG += (double)((((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G));
            TMP_BB += (double)((((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B));
            TMP_RG += (double)((((*(YIN + i))) - U_R) * (((*(UIN + i))) - U_G));
            TMP_RB += (double)((((*(YIN + i))) - U_R) * (((*(VIN + i))) - U_B));
            TMP_GB += (double)((((*(UIN + i))) - U_G) * (((*(VIN + i))) - U_B));
          }
        }
        TMP_RR /= (double)(PT[MEN][1].INDEXNUM);
        TMP_GG /= (double)(PT[MEN][1].INDEXNUM);
        TMP_BB /= (double)(PT[MEN][1].INDEXNUM);
        TMP_RG /= (double)(PT[MEN][1].INDEXNUM);
        TMP_RB /= (double)(PT[MEN][1].INDEXNUM);
        TMP_GB /= (double)(PT[MEN][1].INDEXNUM);

        // DIVIDENUM-1番目の分轄の軸を求める（固有ベクトルの計算）
        A[0][0] = TMP_RR;
        A[0][1] = TMP_RG;
        A[1][0] = TMP_RG;
        A[1][1] = TMP_GG;
        A[2][2] = TMP_BB;
        A[0][2] = TMP_RB;
        A[2][0] = TMP_RB;
        A[1][2] = TMP_GB;
        A[2][1] = TMP_GB;
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
            if (*(INDEX + i) == PT[MEN][1].INDEXNO) {
              *(RR + i) = (double)(((*(YIN + i))) - U_R);
              *(GG + i) = (double)(((*(UIN + i))) - U_G);
              *(BB + i) = (double)(((*(VIN + i))) - U_B);
            }
          }
          V[0] = X1[0][Y];
          V[1] = X1[1][Y];
          V[2] = X1[2][Y];
          for (int i = 0; i < hsize * vsize; i++) {
            kaiten(V, *(RR + i), *(GG + i), *(BB + i), RRR + i, GGG + i,
                   BBB + i);
          }
          // for(i=0,l=0;i<hsize*vsize;i++){
          // if(*(index3+i)!=-1){
          //*(RRRR+l) = (int)(*(RRR+i)+0.5);
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
            if (*(INDEX + i) == PT[MEN][1].INDEXNO) {
              *(RRRRRR + l2) = *(RRR + i);
              *(GGGGGG + l2) = *(GGG + i);
              *(BBBBBB + l2) = *(BBB + i);
              l2++;
            }
          }
          MAXX = -99.9e64;
          MINN = 99.9e64;
          for (int i = 0; i < l2; i++) {
            if (MAXX < *(RRRRRR + i)) {
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
        if ((*(INDEX + i)) == j) {
          Y_JYUSHIN[j] += *(YIN + i);
          U_JYUSHIN[j] += *(UIN + i);
          V_JYUSHIN[j] += *(VIN + i);
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
          MINZ = 60e60;
          for (int j = 0; j < IROSUU; j++) {
            Z[j] = ((((*(YIN + i)))) - Y_JYUSHIN[j]) *
                       ((((*(YIN + i)))) - Y_JYUSHIN[j]) +
                   ((((*(UIN + i)))) - U_JYUSHIN[j]) *
                       ((((*(UIN + i)))) - U_JYUSHIN[j]) +
                   ((((*(VIN + i)))) - V_JYUSHIN[j]) *
                       ((((*(VIN + i)))) - V_JYUSHIN[j]);
          }
          for (int j = 0; j < IROSUU; j++) {
            if (Z[j] < MINZ) {
              MINZ = Z[j];
              X = j;
            }
          }
          *(PALETGAZOU + i) = X;
        }
      }
      // debug start
      int KARI[IROSUU];
      if (1) {
        for (int i = 0; i < IROSUU; i++) {
          KARI[i] = 0;
        }
        for (int i = 0; i < hsize * vsize; i++) {
          KARI[*(PALETGAZOU + i)]++;
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
            if ((*(PALETGAZOU + i)) == j) {
              SUM_R[j] += (*(YIN + i));
              SUM_G[j] += (*(UIN + i));
              SUM_B[j] += (*(VIN + i));
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
  for (int i = 0; i < hsize * vsize; i++) {
    // if(index3[i]!=-1){
    MINZ = 60e60;
    for (int j = 0; j < IROSUU; j++) {
      Z[j] =
          ((((*(YIN + i)))) - Y_JYUSHIN[j]) *
              ((((*(YIN + i)))) - Y_JYUSHIN[j]) +
          ((((*(UIN + i)))) - U_JYUSHIN[j]) *
              ((((*(UIN + i)))) - U_JYUSHIN[j]) +
          ((((*(VIN + i)))) - V_JYUSHIN[j]) * ((((*(VIN + i)))) - V_JYUSHIN[j]);
    }
    for (int j = 0; j < IROSUU; j++) {
      if (Z[j] < MINZ) {
        MINZ = Z[j];
        X = j;
      }
    }
    *(PALETGAZOU + i) = X;
    // }
  }
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
  // *(ROUT+i) = REDUCE_R[*(PALETGAZOU+i)];
  // *(GOUT+i) = REDUCE_G[*(PALETGAZOU+i)];
  // *(BOUT+i) = REDUCE_B[*(PALETGAZOU+i)];
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
      *(IIRIN + i) = (*(RIN + i));
      *(IIGIN + i) = (*(GIN + i));
      *(IIBIN + i) = (*(BIN + i));
    }
    int mx;
    if (dither == 0 || dither == 1) {
      mx = hsize + 2;
    } else if (dither == 2 || dither == 3) {
      mx = hsize + 4;
    }
    double *errorR;
    double *errorG;
    double *errorB;
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
    double ee = 0.0;
    double est = 9999999999.9e33;
    int bst = 0;
    int adr = 0;
    if (dither == 0 || dither == 1) {
      for (int i = 0; i < mx * 2; i++) {
        *(errorR + i) = 0.0;
        *(errorG + i) = 0.0;
        *(errorB + i) = 0.0;
      }
    } else if (dither == 2 || dither == 3) {
      for (int i = 0; i < mx * 3; i++) {
        *(errorR + i) = 0.0;
        *(errorG + i) = 0.0;
        *(errorB + i) = 0.0;
      }
    }
    for (int y = 0; y < vsize; y++) {
      for (int x = 0; x < hsize; x++) {
        if (dither == 0 || dither == 1) {
          adr = x + 1;
        } else if (dither == 2 || dither == 3) {
          adr = x + 2;
        }
        // r= (*(IIRIN+x*vsize+y))+((*(errorR+adr))/16);
        // g= (*(IIGIN+x*vsize+y))+((*(errorG+adr))/16);
        // b= (*(IIBIN+x*vsize+y))+((*(errorB+adr))/16);
        if (dither == 0) {
          r = (double)(*(RIN + x * vsize + y)) +
              ((*(errorR + adr)) / (16.0 / (per + 0.0000001)) /*32.0*/);
          g = (double)(*(GIN + x * vsize + y)) +
              ((*(errorG + adr)) / (16.0 / (per + 0.0000001)) /*32.0*/);
          b = (double)(*(BIN + x * vsize + y)) +
              ((*(errorB + adr)) / (16.0 / (per + 0.0000001)) /*32.0*/);
        } else if (dither == 1) {
          r = (double)(*(RIN + x * vsize + y)) +
              ((*(errorR + adr)) / (4.0 / (per + 0.0000001)));
          g = (double)(*(GIN + x * vsize + y)) +
              ((*(errorG + adr)) / (4.0 / (per + 0.0000001)));
          b = (double)(*(BIN + x * vsize + y)) +
              ((*(errorB + adr)) / (4.0 / (per + 0.0000001)));
        } else if (dither == 2) {
          r = (double)(*(RIN + x * vsize + y)) +
              ((*(errorR + adr)) / (42.0 / (per + 0.0000001)));
          g = (double)(*(GIN + x * vsize + y)) +
              ((*(errorG + adr)) / (42.0 / (per + 0.0000001)));
          b = (double)(*(BIN + x * vsize + y)) +
              ((*(errorB + adr)) / (42.0 / (per + 0.0000001)));
        } else if (dither == 3) {
          r = (double)(*(RIN + x * vsize + y)) +
              ((*(errorR + adr)) / (48.0 / (per + 0.0000001)));
          g = (double)(*(GIN + x * vsize + y)) +
              ((*(errorG + adr)) / (48.0 / (per + 0.0000001)));
          b = (double)(*(BIN + x * vsize + y)) +
              ((*(errorB + adr)) / (48.0 / (per + 0.0000001)));
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
        // r= (*(IIRIN+x*vsize+y));//+((*(errorR+adr))/16);
        // g= (*(IIGIN+x*vsize+y));//+((*(errorG+adr))/16);
        // b= (*(IIBIN+x*vsize+y));//+((*(errorB+adr))/16);
        bst = 0;
        est = 999999.9e33;
        for (int i = 0; i < IROSUU; i++) {
          ee = (r - (double)(REDUCE_R[i])) * (r - (double)(REDUCE_R[i])) +
               (g - (double)(REDUCE_G[i])) * (g - (double)(REDUCE_G[i])) +
               (b - (double)(REDUCE_B[i])) * (b - (double)(REDUCE_B[i]));
          if (est > ee) {
            bst = i;
            est = ee;
          }
        }
        // debug start
        // printf("gosa=%d\n",(int)est);
        // debug end
        re = r - (double)(REDUCE_R[bst]);
        ge = g - (double)(REDUCE_G[bst]);
        be = b - (double)(REDUCE_B[bst]);
        // (*(errorR+adr+1)) += re*7;
        // (*(errorG+adr+1)) += ge*7;
        // (*(errorB+adr+1)) += be*7;

        // (*(errorR+adr+mx-1)) += re*3;
        // (*(errorG+adr+mx-1)) += ge*3;
        // (*(errorB+adr+mx-1)) += be*3;

        // (*(errorR+adr+mx)) += re*5;
        // (*(errorG+adr+mx)) += ge*5;
        // (*(errorB+adr+mx)) += be*5;

        // (*(errorR+adr+mx+1)) += re;
        // (*(errorG+adr+mx+1)) += ge;
        // (*(errorB+adr+mx+1)) += be;
        if (dither == 0) {
          (*(errorR + adr + 1)) += re * 7.0;
          (*(errorG + adr + 1)) += ge * 7.0;
          (*(errorB + adr + 1)) += be * 7.0;

          (*(errorR + adr + mx - 1)) += re * 3.0;
          (*(errorG + adr + mx - 1)) += ge * 3.0;
          (*(errorB + adr + mx - 1)) += be * 3.0;

          (*(errorR + adr + mx)) += re * 5.0;
          (*(errorG + adr + mx)) += ge * 5.0;
          (*(errorB + adr + mx)) += be * 5.0;

          (*(errorR + adr + mx + 1)) += re;
          (*(errorG + adr + mx + 1)) += ge;
          (*(errorB + adr + mx + 1)) += be;
          // (*(errorR+adr+mx)) += re*1;
          // (*(errorG+adr+mx)) += ge*1;
          // (*(errorB+adr+mx)) += be*1;
        } else if (dither == 1) {
          (*(errorR + adr + 1)) += re * 2.0;
          (*(errorG + adr + 1)) += ge * 2.0;
          (*(errorB + adr + 1)) += be * 2.0;

          (*(errorR + adr + mx - 1)) += re;
          (*(errorG + adr + mx - 1)) += ge;
          (*(errorB + adr + mx - 1)) += be;

          (*(errorR + adr + mx)) += re;
          (*(errorG + adr + mx)) += ge;
          (*(errorB + adr + mx)) += be;

          // (*(errorR+adr+mx+1)) += re;
          // (*(errorG+adr+mx+1)) += ge;
          // (*(errorB+adr+mx+1)) += be;

        } else if (dither == 2) {
          (*(errorR + adr + 1)) += re * 8.0;
          (*(errorG + adr + 1)) += ge * 8.0;
          (*(errorB + adr + 1)) += be * 8.0;
          (*(errorR + adr + 2)) += re * 4.0;
          (*(errorG + adr + 2)) += ge * 4.0;
          (*(errorB + adr + 2)) += be * 4.0;
          (*(errorR + adr + mx - 2)) += re * 2.0;
          (*(errorG + adr + mx - 2)) += ge * 2.0;
          (*(errorB + adr + mx - 2)) += be * 2.0;
          (*(errorR + adr + mx - 1)) += re * 4.0;
          (*(errorG + adr + mx - 1)) += ge * 4.0;
          (*(errorB + adr + mx - 1)) += be * 4.0;
          (*(errorR + adr + mx)) += re * 8.0;
          (*(errorG + adr + mx)) += ge * 8.0;
          (*(errorB + adr + mx)) += be * 8.0;
          (*(errorR + adr + mx + 1)) += re * 4.0;
          (*(errorG + adr + mx + 1)) += ge * 4.0;
          (*(errorB + adr + mx + 1)) += be * 4.0;
          (*(errorR + adr + mx + 2)) += re * 2.0;
          (*(errorG + adr + mx + 2)) += ge * 2.0;
          (*(errorB + adr + mx + 2)) += be * 2.0;

          (*(errorR + adr + mx * 2 - 2)) += re;
          (*(errorG + adr + mx * 2 - 2)) += ge;
          (*(errorB + adr + mx * 2 - 2)) += be;
          (*(errorR + adr + mx * 2 - 1)) += re * 2.0;
          (*(errorG + adr + mx * 2 - 1)) += ge * 2.0;
          (*(errorB + adr + mx * 2 - 1)) += be * 2.0;
          (*(errorR + adr + mx * 2)) += re * 4.0;
          (*(errorG + adr + mx * 2)) += ge * 4.0;
          (*(errorB + adr + mx * 2)) += be * 4.0;
          (*(errorR + adr + mx * 2 + 1)) += re * 2.0;
          (*(errorG + adr + mx * 2 + 1)) += ge * 2.0;
          (*(errorB + adr + mx * 2 + 1)) += be * 2.0;
          (*(errorR + adr + mx * 2 + 2)) += re;
          (*(errorG + adr + mx * 2 + 2)) += ge;
          (*(errorB + adr + mx * 2 + 2)) += be;
        } else if (dither == 3) {
          (*(errorR + adr + 1)) += re * 7.0;
          (*(errorG + adr + 1)) += ge * 7.0;
          (*(errorB + adr + 1)) += be * 7.0;
          (*(errorR + adr + 2)) += re * 5.0;
          (*(errorG + adr + 2)) += ge * 5.0;
          (*(errorB + adr + 2)) += be * 5.0;
          (*(errorR + adr + mx - 2)) += re * 3.0;
          (*(errorG + adr + mx - 2)) += ge * 3.0;
          (*(errorB + adr + mx - 2)) += be * 3.0;
          (*(errorR + adr + mx - 1)) += re * 4.0;
          (*(errorG + adr + mx - 1)) += ge * 4.0;
          (*(errorB + adr + mx - 1)) += be * 4.0;
          (*(errorR + adr + mx)) += re * 7.0;
          (*(errorG + adr + mx)) += ge * 7.0;
          (*(errorB + adr + mx)) += be * 7.0;
          (*(errorR + adr + mx + 1)) += re * 5.0;
          (*(errorG + adr + mx + 1)) += ge * 5.0;
          (*(errorB + adr + mx + 1)) += be * 5.0;
          (*(errorR + adr + mx + 2)) += re * 3.0;
          (*(errorG + adr + mx + 2)) += ge * 3.0;
          (*(errorB + adr + mx + 2)) += be * 3.0;

          (*(errorR + adr + mx * 2 - 2)) += re;
          (*(errorG + adr + mx * 2 - 2)) += ge;
          (*(errorB + adr + mx * 2 - 2)) += be;
          (*(errorR + adr + mx * 2 - 1)) += re * 3.0;
          (*(errorG + adr + mx * 2 - 1)) += ge * 3.0;
          (*(errorB + adr + mx * 2 - 1)) += be * 3.0;
          (*(errorR + adr + mx * 2)) += re * 5.0;
          (*(errorG + adr + mx * 2)) += ge * 5.0;
          (*(errorB + adr + mx * 2)) += be * 5.0;
          (*(errorR + adr + mx * 2 + 1)) += re * 3.0;
          (*(errorG + adr + mx * 2 + 1)) += ge * 3.0;
          (*(errorB + adr + mx * 2 + 1)) += be * 3.0;
          (*(errorR + adr + mx * 2 + 2)) += re;
          (*(errorG + adr + mx * 2 + 2)) += ge;
          (*(errorB + adr + mx * 2 + 2)) += be;
        }

        (*(IIRIN + x * vsize + y)) = (REDUCE_R[bst]);
        (*(IIGIN + x * vsize + y)) = (REDUCE_G[bst]);
        (*(IIBIN + x * vsize + y)) = (REDUCE_B[bst]);
        // (*(IIRIN+x*vsize+y)) = (int)(REDUCE_R[*(PALETGAZOU+x*vsize+y)]);
        // (*(IIGIN+x*vsize+y)) = (int)(REDUCE_G[*(PALETGAZOU+x*vsize+y)]);
        // (*(IIBIN+x*vsize+y)) = (int)(REDUCE_B[*(PALETGAZOU+x*vsize+y)]);
      }
      if (dither == 0 || dither == 1) {
        for (int j = 0; j < mx; j++) {
          *(errorR + j) = *(errorR + j + mx);
          *(errorG + j) = *(errorG + j + mx);
          *(errorB + j) = *(errorB + j + mx);
          *(errorR + j + mx) = 0.0;
          *(errorG + j + mx) = 0.0;
          *(errorB + j + mx) = 0.0;
        }
      } else if (dither == 2 || dither == 3) {
        for (int j = 0; j < mx; j++) {
          *(errorR + j) = *(errorR + j + mx);
          *(errorG + j) = *(errorG + j + mx);
          *(errorB + j) = *(errorB + j + mx);
          *(errorR + j + mx) = *(errorR + j + 2 * mx);
          *(errorG + j + mx) = *(errorG + j + 2 * mx);
          *(errorB + j + mx) = *(errorB + j + 2 * mx);
          *(errorR + j + 2 * mx) = 0.0;
          *(errorG + j + 2 * mx) = 0.0;
          *(errorB + j + 2 * mx) = 0.0;
        }
      }
    }

    // *(PALETGAZOU+i)
    // IIRINに画像データがはいっている

    for (int i = 0; i < hsize * vsize; i++) {
      for (int j = 0; j < IROSUU; j++) {
        if (((*(IIRIN + i)) == (REDUCE_R[j])) &&
            ((*(IIGIN + i)) == (REDUCE_G[j])) &&
            ((*(IIBIN + i)) == (REDUCE_B[j]))) {
          *(PALETGAZOU + i) = j;
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
      *(IIRIN + i) = (*(VIN + i));
      *(IIGIN + i) = (*(YIN + i));
      *(IIBIN + i) = (*(UIN + i));
    }
    int mx;
    if (dither == 0 || dither == 1) {
      mx = hsize + 2;
    } else if (dither == 2 || dither == 3) {
      mx = hsize + 4;
    }
    double *errorR;
    double *errorG;
    double *errorB;
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
    double ee = 0.0;
    double est = 9999999999.9e33;
    int bst = 0;
    int adr = 0;
    if (dither == 0 || dither == 1) {
      for (int i = 0; i < mx * 2; i++) {
        *(errorR + i) = 0.0;
        *(errorG + i) = 0.0;
        *(errorB + i) = 0.0;
      }
    } else if (dither == 2 || dither == 3) {
      for (int i = 0; i < mx * 3; i++) {
        *(errorR + i) = 0.0;
        *(errorG + i) = 0.0;
        *(errorB + i) = 0.0;
      }
    }
    for (int y = 0; y < vsize; y++) {
      for (int x = 0; x < hsize; x++) {
        if (dither == 0 || dither == 1) {
          adr = x + 1;
        } else if (dither == 2 || dither == 3) {
          adr = x + 2;
        }
        // r= (*(IIRIN+x*vsize+y))+((*(errorR+adr))/16);
        // g= (*(IIGIN+x*vsize+y))+((*(errorG+adr))/16);
        // b= (*(IIBIN+x*vsize+y))+((*(errorB+adr))/16);
        if (dither == 0) {
          r = (double)(*(VIN + x * vsize + y)) +
              ((*(errorR + adr)) / (16.0 / (per + 0.0000001)) /*32.0*/);
          g = (double)(*(YIN + x * vsize + y)) +
              ((*(errorG + adr)) / (16.0 / (per + 0.0000001)) /*32.0*/);
          b = (double)(*(UIN + x * vsize + y)) +
              ((*(errorB + adr)) / (16.0 / (per + 0.0000001)) /*32.0*/);
        } else if (dither == 1) {
          r = (double)(*(VIN + x * vsize + y)) +
              ((*(errorR + adr)) / (4.0 / (per + 0.0000001)));
          g = (double)(*(YIN + x * vsize + y)) +
              ((*(errorG + adr)) / (4.0 / (per + 0.0000001)));
          b = (double)(*(UIN + x * vsize + y)) +
              ((*(errorB + adr)) / (4.0 / (per + 0.0000001)));
        } else if (dither == 2) {
          r = (double)(*(VIN + x * vsize + y)) +
              ((*(errorR + adr)) / (42.0 / (per + 0.0000001)));
          g = (double)(*(YIN + x * vsize + y)) +
              ((*(errorG + adr)) / (42.0 / (per + 0.0000001)));
          b = (double)(*(UIN + x * vsize + y)) +
              ((*(errorB + adr)) / (42.0 / (per + 0.0000001)));
        } else if (dither == 3) {
          r = (double)(*(VIN + x * vsize + y)) +
              ((*(errorR + adr)) / (48.0 / (per + 0.0000001)));
          g = (double)(*(YIN + x * vsize + y)) +
              ((*(errorG + adr)) / (48.0 / (per + 0.0000001)));
          b = (double)(*(UIN + x * vsize + y)) +
              ((*(errorB + adr)) / (48.0 / (per + 0.0000001)));
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
        // r= (*(IIRIN+x*vsize+y));//+((*(errorR+adr))/16);
        // g= (*(IIGIN+x*vsize+y));//+((*(errorG+adr))/16);
        // b= (*(IIBIN+x*vsize+y));//+((*(errorB+adr))/16);
        bst = 0;
        est = 999999.9e33;
        for (int i = 0; i < IROSUU; i++) {
          ee = (r - (double)(V_JYUSHIN1[i])) * (r - (double)(V_JYUSHIN1[i])) +
               (g - (double)(Y_JYUSHIN1[i])) * (g - (double)(Y_JYUSHIN1[i])) +
               (b - (double)(U_JYUSHIN1[i])) * (b - (double)(U_JYUSHIN1[i]));
          if (est > ee) {
            bst = i;
            est = ee;
          }
        }
        // debug start
        // printf("gosa=%d\n",(int)est);
        // debug end
        re = r - (double)(V_JYUSHIN1[bst]);
        ge = g - (double)(Y_JYUSHIN1[bst]);
        be = b - (double)(U_JYUSHIN1[bst]);
        // (*(errorR+adr+1)) += re*7;
        // (*(errorG+adr+1)) += ge*7;
        // (*(errorB+adr+1)) += be*7;

        // (*(errorR+adr+mx-1)) += re*3;
        // (*(errorG+adr+mx-1)) += ge*3;
        // (*(errorB+adr+mx-1)) += be*3;

        // (*(errorR+adr+mx)) += re*5;
        // (*(errorG+adr+mx)) += ge*5;
        // (*(errorB+adr+mx)) += be*5;

        // (*(errorR+adr+mx+1)) += re;
        // (*(errorG+adr+mx+1)) += ge;
        // (*(errorB+adr+mx+1)) += be;
        if (dither == 0) {
          (*(errorR + adr + 1)) += re * 7.0;
          (*(errorG + adr + 1)) += ge * 7.0;
          (*(errorB + adr + 1)) += be * 7.0;

          (*(errorR + adr + mx - 1)) += re * 3.0;
          (*(errorG + adr + mx - 1)) += ge * 3.0;
          (*(errorB + adr + mx - 1)) += be * 3.0;

          (*(errorR + adr + mx)) += re * 5.0;
          (*(errorG + adr + mx)) += ge * 5.0;
          (*(errorB + adr + mx)) += be * 5.0;

          (*(errorR + adr + mx + 1)) += re;
          (*(errorG + adr + mx + 1)) += ge;
          (*(errorB + adr + mx + 1)) += be;
          // (*(errorR+adr+mx)) += re*1;
          // (*(errorG+adr+mx)) += ge*1;
          // (*(errorB+adr+mx)) += be*1;
        } else if (dither == 1) {
          (*(errorR + adr + 1)) += re * 2.0;
          (*(errorG + adr + 1)) += ge * 2.0;
          (*(errorB + adr + 1)) += be * 2.0;

          (*(errorR + adr + mx - 1)) += re;
          (*(errorG + adr + mx - 1)) += ge;
          (*(errorB + adr + mx - 1)) += be;

          (*(errorR + adr + mx)) += re;
          (*(errorG + adr + mx)) += ge;
          (*(errorB + adr + mx)) += be;

          // (*(errorR+adr+mx+1)) += re;
          // (*(errorG+adr+mx+1)) += ge;
          // (*(errorB+adr+mx+1)) += be;

        } else if (dither == 2) {
          (*(errorR + adr + 1)) += re * 8.0;
          (*(errorG + adr + 1)) += ge * 8.0;
          (*(errorB + adr + 1)) += be * 8.0;
          (*(errorR + adr + 2)) += re * 4.0;
          (*(errorG + adr + 2)) += ge * 4.0;
          (*(errorB + adr + 2)) += be * 4.0;
          (*(errorR + adr + mx - 2)) += re * 2.0;
          (*(errorG + adr + mx - 2)) += ge * 2.0;
          (*(errorB + adr + mx - 2)) += be * 2.0;
          (*(errorR + adr + mx - 1)) += re * 4.0;
          (*(errorG + adr + mx - 1)) += ge * 4.0;
          (*(errorB + adr + mx - 1)) += be * 4.0;
          (*(errorR + adr + mx)) += re * 8.0;
          (*(errorG + adr + mx)) += ge * 8.0;
          (*(errorB + adr + mx)) += be * 8.0;
          (*(errorR + adr + mx + 1)) += re * 4.0;
          (*(errorG + adr + mx + 1)) += ge * 4.0;
          (*(errorB + adr + mx + 1)) += be * 4.0;
          (*(errorR + adr + mx + 2)) += re * 2.0;
          (*(errorG + adr + mx + 2)) += ge * 2.0;
          (*(errorB + adr + mx + 2)) += be * 2.0;

          (*(errorR + adr + mx * 2 - 2)) += re;
          (*(errorG + adr + mx * 2 - 2)) += ge;
          (*(errorB + adr + mx * 2 - 2)) += be;
          (*(errorR + adr + mx * 2 - 1)) += re * 2.0;
          (*(errorG + adr + mx * 2 - 1)) += ge * 2.0;
          (*(errorB + adr + mx * 2 - 1)) += be * 2.0;
          (*(errorR + adr + mx * 2)) += re * 4.0;
          (*(errorG + adr + mx * 2)) += ge * 4.0;
          (*(errorB + adr + mx * 2)) += be * 4.0;
          (*(errorR + adr + mx * 2 + 1)) += re * 2.0;
          (*(errorG + adr + mx * 2 + 1)) += ge * 2.0;
          (*(errorB + adr + mx * 2 + 1)) += be * 2.0;
          (*(errorR + adr + mx * 2 + 2)) += re;
          (*(errorG + adr + mx * 2 + 2)) += ge;
          (*(errorB + adr + mx * 2 + 2)) += be;
        } else if (dither == 3) {
          (*(errorR + adr + 1)) += re * 7.0;
          (*(errorG + adr + 1)) += ge * 7.0;
          (*(errorB + adr + 1)) += be * 7.0;
          (*(errorR + adr + 2)) += re * 5.0;
          (*(errorG + adr + 2)) += ge * 5.0;
          (*(errorB + adr + 2)) += be * 5.0;
          (*(errorR + adr + mx - 2)) += re * 3.0;
          (*(errorG + adr + mx - 2)) += ge * 3.0;
          (*(errorB + adr + mx - 2)) += be * 3.0;
          (*(errorR + adr + mx - 1)) += re * 4.0;
          (*(errorG + adr + mx - 1)) += ge * 4.0;
          (*(errorB + adr + mx - 1)) += be * 4.0;
          (*(errorR + adr + mx)) += re * 7.0;
          (*(errorG + adr + mx)) += ge * 7.0;
          (*(errorB + adr + mx)) += be * 7.0;
          (*(errorR + adr + mx + 1)) += re * 5.0;
          (*(errorG + adr + mx + 1)) += ge * 5.0;
          (*(errorB + adr + mx + 1)) += be * 5.0;
          (*(errorR + adr + mx + 2)) += re * 3.0;
          (*(errorG + adr + mx + 2)) += ge * 3.0;
          (*(errorB + adr + mx + 2)) += be * 3.0;

          (*(errorR + adr + mx * 2 - 2)) += re;
          (*(errorG + adr + mx * 2 - 2)) += ge;
          (*(errorB + adr + mx * 2 - 2)) += be;
          (*(errorR + adr + mx * 2 - 1)) += re * 3.0;
          (*(errorG + adr + mx * 2 - 1)) += ge * 3.0;
          (*(errorB + adr + mx * 2 - 1)) += be * 3.0;
          (*(errorR + adr + mx * 2)) += re * 5.0;
          (*(errorG + adr + mx * 2)) += ge * 5.0;
          (*(errorB + adr + mx * 2)) += be * 5.0;
          (*(errorR + adr + mx * 2 + 1)) += re * 3.0;
          (*(errorG + adr + mx * 2 + 1)) += ge * 3.0;
          (*(errorB + adr + mx * 2 + 1)) += be * 3.0;
          (*(errorR + adr + mx * 2 + 2)) += re;
          (*(errorG + adr + mx * 2 + 2)) += ge;
          (*(errorB + adr + mx * 2 + 2)) += be;
        }

        (*(IIRIN + x * vsize + y)) = (V_JYUSHIN1[bst]);
        (*(IIGIN + x * vsize + y)) = (Y_JYUSHIN1[bst]);
        (*(IIBIN + x * vsize + y)) = (U_JYUSHIN1[bst]);
        // (*(IIRIN+x*vsize+y)) = (int)(REDUCE_R[*(PALETGAZOU+x*vsize+y)]);
        // (*(IIGIN+x*vsize+y)) = (int)(REDUCE_G[*(PALETGAZOU+x*vsize+y)]);
        // (*(IIBIN+x*vsize+y)) = (int)(REDUCE_B[*(PALETGAZOU+x*vsize+y)]);
      }
      if (dither == 0 || dither == 1) {
        for (int j = 0; j < mx; j++) {
          *(errorR + j) = *(errorR + j + mx);
          *(errorG + j) = *(errorG + j + mx);
          *(errorB + j) = *(errorB + j + mx);
          *(errorR + j + mx) = 0.0;
          *(errorG + j + mx) = 0.0;
          *(errorB + j + mx) = 0.0;
        }
      } else if (dither == 2 || dither == 3) {
        for (int j = 0; j < mx; j++) {
          *(errorR + j) = *(errorR + j + mx);
          *(errorG + j) = *(errorG + j + mx);
          *(errorB + j) = *(errorB + j + mx);
          *(errorR + j + mx) = *(errorR + j + 2 * mx);
          *(errorG + j + mx) = *(errorG + j + 2 * mx);
          *(errorB + j + mx) = *(errorB + j + 2 * mx);
          *(errorR + j + 2 * mx) = 0.0;
          *(errorG + j + 2 * mx) = 0.0;
          *(errorB + j + 2 * mx) = 0.0;
        }
      }
    }

    // *(PALETGAZOU+i)
    // IIRINに画像データがはいっている

    for (int i = 0; i < hsize * vsize; i++) {
      for (int j = 0; j < IROSUU; j++) {
        if (((int)((*(IIRIN + i)) /* *65384.0*/) ==
             (int)((V_JYUSHIN1[j]) /* *65384.0*/)) &&
            ((int)((*(IIGIN + i)) /* *65384.0*/) ==
             (int)((Y_JYUSHIN1[j]) /**65384.0*/)) &&
            ((int)((*(IIBIN + i)) /* *65384.0*/) ==
             (int)((U_JYUSHIN1[j]) /* *65384.0*/))) {
          *(PALETGAZOU + i) = j;
        }
      }
    }
  } // edon if end

} // median cut end
