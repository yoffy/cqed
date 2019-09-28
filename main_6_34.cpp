// ./*** -i 入力24bitBMPファイル名 -o 出力bmpファイル名（.bmpはいらない）

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct control_data {
  unsigned int edge;
  unsigned int sizey;
  int et;
  double a;
  int edon;
  int dither;
  int lpfon;
  double pw;
  double pw2;
  int bun;
  int kmon;
  double lpfk;
  double per;
  double bai;
  int omh;
  int cs;
  double t;
  int div;
  int hasyori;
  double r;
  double g;
  double b;
  double gamma;
  double bgamma;
  double rgamma;
  double ygamma;
  char InputFileName[128];
  char OutputFileName[128];
} CNTL;
int hsize;
int vsize;
char inputfile[128];
char outputfile[128];

//#include "mediancut8.c"
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
void print_help(CNTL *ptr) {
  fprintf(stderr, "Usage: bmp_color_reducer[options]\n");
  fprintf(stderr, "\toptions: \n");
  fprintf(stderr, "\t -help : \n");
  fprintf(stderr, "\t -i : input file name(BMP 24bit true color)\n");
  fprintf(stderr, "\t -o : outputfilename(BMP 8bit pallet color [kakutyousi "
                  ".bmp ga jidoude tukerareru]\n");
  fprintf(stderr,
          "\t -edge : edge kensyutu no houhou sentaku 0:rinsetu gaso kan sabun "
          "sqrt(suihei^2+suityoku^2) 1,2=rinsetu gasokansabun dokuritu 3=sobel "
          "filter 4=edgerasisa 5=sobel filter2 6=rinsetu gaso kan sabun sqrt2, "
          " 4 wo sentaku sitatiki kmean ha off surukoto\n");
  // fprintf( stderr, "\t -et : edge threshold -edge 0,1,2,3,5 notoki
  // (0<=et<=255) et ijyou wo edge to minasu -et 512 tositeokeba edge
  // kensyutusinai gennsyoku ga dekiru\n" );
  fprintf(stderr,
          "\t -et : edge threshold (0<=et<=255) et ijyou wo edge to minasu -et "
          "512 tositeokeba edge kensyutusinai gennsyoku ga dekiru\n");
  fprintf(stderr, "\t -bai : edgerasisa -edge 4 no toki yuukou edge no ookisa "
                  "wo bai bai suru tuneni 1.0 suisyou 0.0 de edgerasisa off 1 "
                  "yori ookihodo yowai edge made mushisuru  \n");
  fprintf(stderr, "\t -kmon : kmean onoff switch 1:on 0:off\n");
  fprintf(stderr, "\t -edon : error diffusion onoff switch 1:on ,RGB space de "
                  "dither 0:off 2:on ,henkan iro kuukan de dither \n");
  fprintf(stderr, "\t -r : R (v,B) zahyoujiku kakudairitu R3 V4 v1 B1 G4\n");
  fprintf(stderr, "\t -g : G (L,Y) zahyoujiku kakudairitu G6 Y4 L1 L1 C4\n");
  fprintf(stderr, "\t -b : B (u,A) zahyoujiku kakudairitu B1 U3 u1 A1 S4\n");
  fprintf(stderr,
          "\t -cs : color space 0:RGB 1:Lab 2:Luv 3:YUV 4.GCS(iZYINSgamma) "
          "5.YUV(S ji gamma -a option yuukou) 6.RGB(S ji gamma -a option "
          "yuukou) 7.RGB(iZyins gamma) 8.YUV(iZyins gamma) 9.YUV(izyins "
          "irokuukan gamma ha nashi)\n");
  fprintf(stderr,
          "\t -gamma :  RGB no green gamma.  1yori ookiito akaruibubunn ga "
          "komakaku iro wo kubetsushi kuraibubunn ga oozappa ni naru\n");
  fprintf(stderr,
          "\t -rgamma : RGB no red gamma.  1yori ookiito akaruibubunn ga "
          "komakaku iro wo kubetsushi kuraibubunn ga oozappa ni naru\n");
  fprintf(stderr,
          "\t -bgamma : RGB no blue gamma. 1yori ookiito akaruibubunn ga "
          "komakaku iro wo kubetsushi kuraibubunn ga oozappa ni naru\n");
  fprintf(stderr, "\t -ygamma : YUV only Y no gamma. LAB.Luv deha 1.0 ni "
                  "surukoto.1yori ookiito akaruibubunn ga komakaku iro wo "
                  "kubetsushi kuraibubunn ga oozappa ni naru\n");
  fprintf(stderr, "\t -omh : cut method 0:ohtsu 1:median 2:heikin "
                  "3.syousuu_ohtsu 4.syousuu_ohtsu2\n");
  fprintf(stderr, "\t -hasyori : RGB mode only. dark color threshold 20 kurai "
                  "ga ii. RGB igaiha 0 ni surukoto.\n");
  fprintf(stderr, "\t -dither : dither mode 0:floyd steinburg 1:Sierra Lite 2: "
                  "Stucki 3.Jarvis,Judice and Nink\n");
  fprintf(stderr,
          "\t -per : dither kyoudo 0 kara 1 made , 0.57 recommended. \n");
  // fprintf( stderr, "\t -div : 1. heikin karano kyori no 2jyouwa ga ookii mono
  // wo bunkatu.\n" ); fprintf( stderr, "\t      : 0. heikin karano kyori no wa
  // ga ookii mono wo bunkatu.\n" );
  fprintf(stderr, "\t -pw  : heikin kara gasoti no kyori wo nan jyou suruka "
                  "1.0 de kyori 2.0 de kyori no 2jyou \n");
  fprintf(stderr,
          "\t -pw2  : -edge 4 edgerasisa no toki yuukou 1.0 yori ookiito "
          "ookina edge ga iro sentaku no sai musi sareru. 1.0yori tiisai toki "
          "ookina edge dakedenaku tiisana edge mo musisareyasui \n");
  fprintf(stderr, "\t -bun  : bunkatsu ryouiki sentaku hou 0:maxdistance*num "
                  "max 1.distance^pw sum max 2.max bunsan*num max3.max bunsan "
                  "4.max bunsan no houkou no MAX-MIN ga max 5. max bunsan no "
                  "houkou no MAX-MIN *num ga max \n");
  fprintf(stderr, "\t -lpfon  : LPF SW 1:on 0:off \n");
  fprintf(
      stderr,
      "\t -lpfk   : LPF kyoudo 0.5(tuyoi) kara 1.0(yowai) 1.0 de lpfoff \n");
  fprintf(stderr, "\t -t  : color space kido jiku wo tyuusin ni kaiten sono "
                  "kakudo DEGREE 0~90\n");
  fprintf(stderr, "\t -a  : S ji tone curve gamma keisuu\n");
  fprintf(stderr,
          "\t -edge 0 notoki rinsetugaso kan sabun tate yoko square root\n");
  fprintf(stderr,
          "\t -edge 1,2 notoki rinsetugaso kan sabun tate yoko dokuritu\n");
  fprintf(stderr, "\t -edge 3 notoki sobel tate yoko square root\n");

  fprintf(stderr,
          "DEFAULT PARAM. "
          "edge=%d,et=%d,bai=%f,kmon=%d,edon=%d,r=%f,g=%f,b=%f,cs=%d,gamma=%f,"
          "rgamma=%f,bgamma=%f,ygamma=%f,omh=%d,hasyori=%d,dither=%d,per=%f,pw="
          "%f,pw2=%f,bun=%d,lpfon=%d,t=%f,a=%f,lpfk=%f\n",
          ptr->edge, ptr->et, ptr->bai, ptr->kmon, ptr->edon, ptr->r, ptr->g,
          ptr->b, ptr->cs, ptr->gamma, ptr->rgamma, ptr->bgamma, ptr->ygamma,
          ptr->omh, ptr->hasyori, ptr->dither, ptr->per, ptr->pw, ptr->pw2,
          ptr->bun, ptr->lpfon, ptr->t, ptr->a, ptr->lpfk);
}

void option_set(int argc, char *argv[], CNTL *ptr) {
  int i;
  for (i = 1; i < argc; i++) {
    if (argv[i][0] != '-') {
      continue;
    }
    if (!strcmp(&argv[i][1], "h")) {
      print_help(ptr);
      exit(0);
    } else if (!strcmp(&argv[i][1], "i")) {
      strcpy(ptr->InputFileName, argv[++i]);
    } else if (!strcmp(&argv[i][1], "o")) {
      strcpy(ptr->OutputFileName, argv[++i]);
    } else if (!strcmp(&argv[i][1], "edge")) {
      ptr->edge = atoi(argv[++i]);
    } else if (!strcmp(&argv[i][1], "edon")) {
      ptr->edon = atoi(argv[++i]);
    } else if (!strcmp(&argv[i][1], "kmon")) {
      ptr->kmon = atoi(argv[++i]);
    } else if (!strcmp(&argv[i][1], "cs")) {
      ptr->cs = atoi(argv[++i]);
    } else if (!strcmp(&argv[i][1], "et")) {
      ptr->et = atoi(argv[++i]);
    } else if (!strcmp(&argv[i][1], "omh")) {
      ptr->omh = atoi(argv[++i]);
    } else if (!strcmp(&argv[i][1], "hasyori")) {
      ptr->hasyori = atoi(argv[++i]);
    } else if (!strcmp(&argv[i][1], "y")) {
      ptr->sizey = atoi(argv[++i]);
    } else if (!strcmp(&argv[i][1], "dither")) {
      ptr->dither = atoi(argv[++i]);
    } else if (!strcmp(&argv[i][1], "div")) {
      ptr->div = atoi(argv[++i]);
    } else if (!strcmp(&argv[i][1], "bun")) {
      ptr->bun = atoi(argv[++i]);
    } else if (!strcmp(&argv[i][1], "lpfon")) {
      ptr->lpfon = atoi(argv[++i]);
    } else if (!strcmp(&argv[i][1], "r")) {
      ptr->r = atof(argv[++i]);
    } else if (!strcmp(&argv[i][1], "g")) {
      ptr->g = atof(argv[++i]);
    } else if (!strcmp(&argv[i][1], "b")) {
      ptr->b = atof(argv[++i]);
    } else if (!strcmp(&argv[i][1], "gamma")) {
      ptr->gamma = atof(argv[++i]);
    } else if (!strcmp(&argv[i][1], "rgamma")) {
      ptr->rgamma = atof(argv[++i]);
    } else if (!strcmp(&argv[i][1], "bgamma")) {
      ptr->bgamma = atof(argv[++i]);
    } else if (!strcmp(&argv[i][1], "ygamma")) {
      ptr->ygamma = atof(argv[++i]);
    } else if (!strcmp(&argv[i][1], "per")) {
      ptr->per = atof(argv[++i]);
    } else if (!strcmp(&argv[i][1], "pw")) {
      ptr->pw = atof(argv[++i]);
    } else if (!strcmp(&argv[i][1], "pw2")) {
      ptr->pw2 = atof(argv[++i]);
    } else if (!strcmp(&argv[i][1], "t")) {
      ptr->t = atof(argv[++i]);
    } else if (!strcmp(&argv[i][1], "a")) {
      ptr->a = atof(argv[++i]);
    } else if (!strcmp(&argv[i][1], "bai")) {
      ptr->bai = atof(argv[++i]);
    } else if (!strcmp(&argv[i][1], "lpfk")) {
      ptr->lpfk = atof(argv[++i]);
    } else {
      exit(0);
    }
  }
}
struct IRO {
  int INDEX;
  int CNT;
  int RGB;
  unsigned char R;
  unsigned char G;
  unsigned char B;
  double CONVV;
  double CONVY;
  double CONVU;
};
// struct IRO *iro;
// iro =(struct IRO*)malloc(sizeof(struct IRO)*hsize*vsize);

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
struct palet *PTP;

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

  double X, Y, Z, a, b, c, d, ud, vd, u0, v0, TEMP, L1;
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

  double X, Y, Z, fX, fY, fZ, tempX, tempY, tempZ;
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

int ohtsu(int NUM, int *X) {
  int i, j, k, l, m, n, o, p;
  // int SUM1;
  // int X[NUM] = {-10,-10,-10,-10,-10,
  // -10,-10,-10,-10,-10,
  // -3,-3,-3,-3,-3,
  // -1,-1,0,1,1,
  // 3,3,3,3,3,
  // 10,10,10,10,10,
  // 10,10,10,10,10};
  int MAX = -10000;
  int MIN = 99999;
  for (i = 0; i < NUM; i++) {
    if (MIN > *(X + i)) {
      MIN = *(X + i);
    }
    if (MAX < *(X + i)) {
      MAX = *(X + i);
    }
  }
  int NODEHANI;
  NODEHANI = MAX - MIN + 1;
  int *HIST;
  HIST = (int *)malloc(sizeof(int) * NODEHANI);
  // ND
  // 0- -11
  // 1- -10

  // 11- 0

  // 22- 11
  for (i = 0; i < NODEHANI; i++) {
    *(HIST + i) = 0;
  }
  int TMP;
  for (i = 0; i < NUM; i++) {
    TMP = *(X + i) - MIN;
    (*(HIST + TMP))++;
  }
  // debug start
  // for(i=0;i<NODEHANI;i++){
  // fprintf(stderr,"NUM ATAI HIST %d %d %d\n",i,i+MIN,*(HIST+i));
  //}
  // while(1);
  // debug end

  float MAX2 = 0.0;
  float AVE1, AVE2;
  int TOTAL1, TOTAL2;
  float HANTEI;
  int THRESH;
  for (i = 1; i < NODEHANI; i++) {
    AVE1 = 0.0;
    TOTAL1 = 0;
    for (j = 0; j < i; j++) {
      AVE1 += (*(HIST + j)) * (j + MIN);
      TOTAL1 += (*(HIST + j));
    }
    AVE1 /= (float)(TOTAL1);
    AVE2 = 0.0;
    TOTAL2 = 0;
    for (j = i; j < NODEHANI; j++) {
      AVE2 += (*(HIST + j)) * (j + MIN);
      TOTAL2 += (*(HIST + j));
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
  return (THRESH);
  // return (0);
} // main kansuu end

int ohtsu2(int NUM, double *X, double *Y, double *Z, int omh) {
  int i, j, k, l, m, n, o, p;
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
  double MAX = -10000.0e32;
  double MIN = 99999.0e32;
  for (i = 0; i < NUM; i++) {
    if (MIN > *(X + i)) {
      MIN = *(X + i);
    }
    if (MAX < *(X + i)) {
      MAX = *(X + i);
    }
  }
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

  int NODEHANI;
  NODEHANI = CMAX - CMIN + 1;
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
  // return (0);
} // main kansuu end

int media(int NUM, int *X) {
  int i, j, k, l, m, n, o, p;
  // int SUM1;
  // int X[NUM] = {-10,-10,-10,-10,-10,
  // -10,-10,-10,-10,-10,
  // -3,-3,-3,-3,-3,
  // -1,-1,0,1,1,
  // 3,3,3,3,3,
  // 10,10,10,10,10,
  // 10,10,10,10,10};
  int MAX = -10000;
  int MIN = 99999;
  for (i = 0; i < NUM; i++) {
    if (MIN > *(X + i)) {
      MIN = *(X + i);
    }
    if (MAX < *(X + i)) {
      MAX = *(X + i);
    }
  }
  int NODEHANI;
  NODEHANI = MAX - MIN + 1;
  int *HIST;
  HIST = (int *)malloc(sizeof(int) * NODEHANI);
  // ND
  // 0- -11
  // 1- -10

  // 11- 0

  // 22- 11
  for (i = 0; i < NODEHANI; i++) {
    *(HIST + i) = 0;
  }
  int TMP;
  for (i = 0; i < NUM; i++) {
    TMP = *(X + i) - MIN;
    (*(HIST + TMP))++;
  }
  // debug start
  // for(i=0;i<NODEHANI;i++){
  // fprintf(stderr,"NUM ATAI HIST %d %d %d\n",i,i+MIN,*(HIST+i));
  //}
  // while(1);
  // debug end
  int GOUKEI, MED;
  GOUKEI = 0;
  for (i = 0; i < NODEHANI; i++) {
    GOUKEI += *(HIST + i);
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
      AVE1 += (*(HIST + j)) * (j + MIN);
      TOTAL1 += (*(HIST + j));
    }
    AVE1 /= (float)(TOTAL1);
    AVE2 = 0.0;
    TOTAL2 = 0;
    for (j = i; j < NODEHANI; j++) {
      AVE2 += (*(HIST + j)) * (j + MIN);
      TOTAL2 += (*(HIST + j));
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
  double THIETA, FAI;
  double X1, Y1, Z1, KARI;
  double PI = atan(1.0) * 4.0;
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
  double R1, G1, B1, KARI;
  double PI = atan(1.0) * 4.0;
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
  int RDIV;
  int GDIV;
  int BDIV;
  int i, j, k, l, m, n, o, p;
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
  int i1, i2, ind, ct;
  double TMP_RR;
  double TMP_GG;
  double TMP_BB;
  double TMP_RG;
  double TMP_RB;
  double TMP_GB;
  double IGENMAX;
  int Y;
  double ZZ;
  int *INDEX;
  double U_R;
  double U_G;
  double U_B;
  double MAXINDEXNUM;
  // int KARI;
  int KARI1;
  int KARI2;
  int KARI3;
  int KARI4;
  int KARI5;
  int KARI6;
  int FLAG;
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
    for (i = 0; i < vsize; i++) {
      for (j = 0; j < hsize; j++) {
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
    for (i = 0; i < hsize; i++) {
      for (j = 0; j < vsize; j++) {
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
  for (i = 0; i < hsize * vsize; i++) {
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
  for (i = 0; i < hsize * vsize; i++) {
    RrIN[i] = 255.0 * pow((double)RIN[i] / 255.0, rgamma);
    GgIN[i] = 255.0 * pow((double)GIN[i] / 255.0, gamma);
    BbIN[i] = 255.0 * pow((double)BIN[i] / 255.0, bgamma);
  }

  IRO_THIETA *= PI / 180;
  if (cs == 0) {
    for (i = 0; i < hsize * vsize; i++) {
      *(YIN + i) = *(GgIN + i); //(0.29891*(float)(*(RIN+i))+0.58661*(float)(*(GIN+i))+0.11448*(float)(*(BIN+i))
                                ///*+ 0.5*/);
      *(UIN + i) = *(BbIN + i); //((-0.16874*(float)(*(RIN+i))-0.33126*(float)(*(GIN+i))+0.50000*(float)(*(BIN+i)))
                                //+ 128.0);
      *(VIN + i) = *(RrIN + i); //(0.50000*(float)(*(RIN+i))-0.41869*(float)(*(GIN+i))-0.08131*(float)(*(BIN+i))
                                //+ 128.0);
    }
  } else if (cs == 1) {
    for (i = 0; i < hsize * vsize; i++) {
      RGBtoLAB(RrIN[i], GgIN[i], BbIN[i], VIN + i, YIN + i, UIN + i); //(0.29891*(float)(*(RIN+i))+0.58661*(float)(*(GIN+i))+0.11448*(float)(*(BIN+i))
                                                                      //+ 0.5);
      //*(UIN+i) =
      //*(BIN+i);//((-0.16874*(float)(*(RIN+i))-0.33126*(float)(*(GIN+i))+0.50000*(float)(*(BIN+i)))
      //+ 128.5);
      //*(VIN+i) =
      //*(RIN+i);//(0.50000*(float)(*(RIN+i))-0.41869*(float)(*(GIN+i))-0.08131*(float)(*(BIN+i))
      //+ 128.5);
    }
  } else if (cs == 2) {
    for (i = 0; i < hsize * vsize; i++) {
      RGBtoLuv(RrIN[i], GgIN[i], BbIN[i], VIN + i, YIN + i, UIN + i); //(0.29891*(float)(*(RIN+i))+0.58661*(float)(*(GIN+i))+0.11448*(float)(*(BIN+i))
                                                                      //+ 0.5);
      //*(UIN+i) =
      //*(BIN+i);//((-0.16874*(float)(*(RIN+i))-0.33126*(float)(*(GIN+i))+0.50000*(float)(*(BIN+i)))
      //+ 128.5);
      //*(VIN+i) =
      //*(RIN+i);//(0.50000*(float)(*(RIN+i))-0.41869*(float)(*(GIN+i))-0.08131*(float)(*(BIN+i))
      //+ 128.5);
    }
  } else if (cs == 3 || cs == 9) {
    for (i = 0; i < hsize * vsize; i++) {
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
    for (i = 0; i < hsize * vsize; i++) {
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
    for (i = 0; i < hsize * vsize; i++) {
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
    for (i = 0; i < hsize * vsize; i++) {
      RrIN[i] = (255.0 / (1.0 + exp(-aaaa * ((double)RrIN[i] - 127.5))) -
                 255.0 / (1.0 + exp(-aaaa * (-127.5)))) *
                255.0 / (255.0 - 2.0 * 255.0 / (1.0 + exp(-aaaa * (-127.5))));
      GgIN[i] = (255.0 / (1.0 + exp(-aaaa * ((double)GgIN[i] - 127.5))) -
                 255.0 / (1.0 + exp(-aaaa * (-127.5)))) *
                255.0 / (255.0 - 2.0 * 255.0 / (1.0 + exp(-aaaa * (-127.5))));
      BbIN[i] = (255.0 / (1.0 + exp(-aaaa * ((double)BbIN[i] - 127.5))) -
                 255.0 / (1.0 + exp(-aaaa * (-127.5)))) *
                255.0 / (255.0 - 2.0 * 255.0 / (1.0 + exp(-aaaa * (-127.5))));

      *(YIN + i) = *(GgIN + i); //(0.29891*(float)(*(RIN+i))+0.58661*(float)(*(GIN+i))+0.11448*(float)(*(BIN+i))
                                ///*+ 0.5*/);
      *(UIN + i) = *(BbIN + i); //((-0.16874*(float)(*(RIN+i))-0.33126*(float)(*(GIN+i))+0.50000*(float)(*(BIN+i)))
                                //+ 128.0);
      *(VIN + i) = *(RrIN + i); //(0.50000*(float)(*(RIN+i))-0.41869*(float)(*(GIN+i))-0.08131*(float)(*(BIN+i))
                                //+ 128.0);
    }
  } else if (cs == 7) {
    for (i = 0; i < hsize * vsize; i++) {
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
    for (i = 0; i < hsize * vsize; i++) {
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
      //BbIN[i];//((-0.1350*(float)(*(RrIN+i))-0.2650*(float)(*(GgIN+i))+0.40000*(float)(*(BbIN+i))));
      //*(VIN+i) =
      //RrIN[i];//(0.4000*(float)(*(RrIN+i))-0.3346*(float)(*(GgIN+i))-0.0653*(float)(*(BbIN+i)));
    }
  }
  double *KARIU, *KARIV;
  KARIU = (double *)malloc(sizeof(double) * hsize * vsize);
  KARIV = (double *)malloc(sizeof(double) * hsize * vsize);
  if (cs == 4 || cs == 9) {
    for (i = 0; i < hsize * vsize; i++) {
      KARIU[i] = UIN[i] * /*cos(IRO_THIETA)*/ 0.941996495 +
                 VIN[i] * /*sin(IRO_THIETA)*/ 0.13632709;
      KARIV[i] = UIN[i] * /*sin(IRO_THIETA)*/ 0.118202579 +
                 VIN[i] * /*cos(IRO_THIETA)*/ 1.201426141;
      UIN[i] = KARIU[i];
      VIN[i] = KARIV[i];
    }
  } else {
    for (i = 0; i < hsize * vsize; i++) {
      KARIU[i] = UIN[i] * cos(IRO_THIETA) - VIN[i] * sin(IRO_THIETA);
      KARIV[i] = UIN[i] * sin(IRO_THIETA) + VIN[i] * cos(IRO_THIETA);
      UIN[i] = KARIU[i];
      VIN[i] = KARIV[i];
    }
  }
  for (i = 0; i < hsize * vsize; i++) {
    *(YIN + i) = 255.0 * pow(*(YIN + i) / 255.0, ygamma);
  }

  for (i = 0; i < hsize * vsize; i++) {
    *(YIN + i) *= gmult;
  }
  for (i = 0; i < hsize * vsize; i++) {
    *(VIN + i) *= rmult;
  }
  // for(i=0;i<hsize*vsize;i++){
  //*(VIN+i) = 255.0*rmult*pow(*(VIN+i)/(rmult*255.0),rgamma);
  //}
  for (i = 0; i < hsize * vsize; i++) {
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
  for (i = 0; i < hsize * vsize; i++) {
    *(IRIN + i) = (int)(*(RIN + i));
    *(IGIN + i) = (int)(*(GIN + i));
    *(IBIN + i) = (int)(*(BIN + i));
  }
  int *index3;
  index3 = (int *)malloc(sizeof(int) * hsize * vsize);
  // double EDGERASMAX;
  for (i = 0; i < hsize * vsize; i++) {
    HHEDGER[i] = 0;
    HHEDGEG[i] = 0;
    HHEDGEB[i] = 0;
    VVEDGER[i] = 0;
    VVEDGEG[i] = 0;
    VVEDGEB[i] = 0;
  }

  int SOBEL[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
  int SOBEL2[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
  if (vvv == 3 || vvv == 4 || vvv == 5) {
    for (i = 0; i < vsize; i++) {   // suityoku
      for (j = 0; j < hsize; j++) { // suihei
        for (k = 0; k < 3; k++) {   // suityoku
          for (l = 0; l < 3; l++) { // suihei
            m = i + k - 1;          //-1~1suiityoku
            n = j + l - 1;          //-1~1suihei
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
    for (i = 0; i < vsize; i++) {
      for (j = 0; j < hsize; j++) {
        k = j - 1;
        l = j + 1;
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
    for (i = 0; i < hsize; i++) {
      for (j = 0; j < vsize; j++) {
        k = j - 1;
        l = j + 1;
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
    for (i = 0, k = 0; i < vsize; i++) {
      for (j = 0; j < hsize; j++) {
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
            k++;
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
            k++;
          } else {
            *(index3 + j * vsize + i) = 1;
          }
        }
        if (vvv == 1) {
          if (HHEDGER[j * vsize + i] >= et || VVEDGER[j * vsize + i] >= et ||
              HHEDGEG[j * vsize + i] >= et || VVEDGEG[j * vsize + i] >= et ||
              HHEDGEB[j * vsize + i] >= et || VVEDGEB[j * vsize + i] >= et) {
            *(index3 + j * vsize + i) = -1;
            k++;
          } else {
            *(index3 + j * vsize + i) = 1;
          }
        }
        if (vvv == 2) {
          if (EDGER[j * vsize + i] == 1 || EDGEG[j * vsize + i] == 1 ||
              EDGEB[j * vsize + i] == 1 || VEDGER[j * vsize + i] == 1 ||
              VEDGEG[j * vsize + i] == 1 || VEDGEB[j * vsize + i] == 1) {
            *(index3 + j * vsize + i) = -1;
            k++;
          } else {
            *(index3 + j * vsize + i) = 1;
          }
        }
      }
    }
  } // cnt if end

  if (vvv == 3 || vvv == 5) {
    for (i = 0, k = 0; i < vsize; i++) {
      for (j = 0; j < hsize; j++) {
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
            k++;
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
            k++;
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
    for (i = 0; i < vsize; i++) {
      for (j = 0; j < hsize; j++) {
        *(index3 + j * vsize + i) = 1;
        REDGE = /*sqrt*/ (
            (float)HHEDGER[j * vsize + i] * (float)HHEDGER[j * vsize + i] +
            (float)VVEDGER[j * vsize + i] *
                (float)VVEDGER[j * vsize + i]) /*/ 4.0 / sqrt(2.0)*/;
        GEDGE = /*sqrt*/ (
            (float)HHEDGEG[j * vsize + i] * (float)HHEDGEG[j * vsize + i] +
            (float)VVEDGEG[j * vsize + i] *
                (float)VVEDGEG[j * vsize + i]) /*/ 4.0 / sqrt(2.0)*/;
        REDGE = /*sqrt*/ (
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

    k = 0;
  }

  // debug start]
  // if(vvv == 4){
  //    printf("EDGE rasisa mode\n");
  //} else {
  printf("EDGE gasosuu=%d\n", k);
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
  RDIV = RD;
  GDIV = GD;
  BDIV = BD;
  if (MEN == 1) {
    NMEN = 0;
  } else {
    NMEN = 1;
  }
  // MEN = 0,NMEN = 1
  U_R = 0.0;
  U_G = 0.0;
  U_B = 0.0;
  for (i = 0; i < hsize * vsize; i++) {
    if (*(index3 + i) != -1) {
      U_R += *(YIN + i);
      U_G += *(UIN + i);
      U_B += *(VIN + i);
    }
  }
  U_R /= (double)(hsize * vsize - k);
  U_G /= (double)(hsize * vsize - k);
  U_B /= (double)(hsize * vsize - k);
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
  for (i = 0; i < hsize * vsize; i++) {
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
  TMP_RR /= (double)(hsize * vsize - k);
  TMP_GG /= (double)(hsize * vsize - k);
  TMP_BB /= (double)(hsize * vsize - k);
  TMP_RG /= (double)(hsize * vsize - k);
  TMP_RB /= (double)(hsize * vsize - k);
  TMP_GB /= (double)(hsize * vsize - k);
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
  for (i = 0; i < 3; i++) {
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
  for (i = 0; i < hsize * vsize; i++) {
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
  for (i = 0; i < hsize * vsize; i++) {
    kaiten(V, *(RR + i), *(GG + i), *(BB + i), RRR + i, GGG + i, BBB + i);
  }
  for (i = 0, l = 0; i < hsize * vsize; i++) {
    if (*(index3 + i) != -1) {
      *(RRRR + l) = (int)(*(RRR + i) + 0.5);
      l++;
    }
  }

  double *RRRRRR, *GGGGGG, *BBBBBB;
  int GASOSUU = hsize * vsize - k;
  if (omh == 3 || omh == 4) {

    RRRRRR = (double *)malloc(sizeof(double) * hsize * vsize);
    GGGGGG = (double *)malloc(sizeof(double) * hsize * vsize);
    BBBBBB = (double *)malloc(sizeof(double) * hsize * vsize);
    for (i = 0, l = 0; i < hsize * vsize; i++) {
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
  int kk;
  kk = k;
  for (i = 0; i < hsize * vsize; i++) {
    *(INDEX + i) = -1;
  }
  for (i = 0, k = 0, l = 0; i < hsize * vsize; i++) {
    if (*(index3 + i) != -1) {
      // ZZ = X1[0][Y]*((double)(((*(YIN+i)))) - ((double)(U_R)+(*XXX))) +
      // X1[1][Y]*((double)(((*(UIN+i)))) - ((double)(U_G)+(*YYY))) +
      // X1[2][Y]*((double)(((*(VIN+i)))) - ((double)(U_B)+(*ZZZ))) ;
      // if(ZZ>=0.0){
      if ((double)(*(RRR + i)) >= DTHRESH) {
        *(INDEX + i) = PT[MEN][0].INDEXNO;
        k++;
      } else {
        *(INDEX + i) = PT[MEN][1].INDEXNO;
        l++;
      }
    }
  }
  PT[MEN][0].INDEXNUM = k;
  PT[MEN][1].INDEXNUM = l;
  // maxdistance tuika start
  U_R = 0.0;
  U_G = 0.0;
  U_B = 0.0;
  // 0側のmaxdistanceを求める
  for (i = 0; i < hsize * vsize; i++) {
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
      for (i = 0; i < hsize * vsize; i++) {
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
          for (i = 0; i < hsize * vsize; i++) {
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
          for (i = 0; i < hsize * vsize; i++) {
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
          for (i = 0; i < hsize * vsize; i++) {
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
          for (i = 0; i < hsize * vsize; i++) {
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
      for (i = 0; i < hsize * vsize; i++) {
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
      for (i = 0; i < 3; i++) {
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

        for (i = 0; i < hsize * vsize; i++) {
          if (*(INDEX + i) == PT[MEN][0].INDEXNO) {
            *(RR + i) = (double)(((*(YIN + i))) - U_R);
            *(GG + i) = (double)(((*(UIN + i))) - U_G);
            *(BB + i) = (double)(((*(VIN + i))) - U_B);
          }
        }
        V[0] = X1[0][Y];
        V[1] = X1[1][Y];
        V[2] = X1[2][Y];
        for (i = 0; i < hsize * vsize; i++) {
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
        for (i = 0, l = 0; i < hsize * vsize; i++) {
          if (*(INDEX + i) == PT[MEN][0].INDEXNO) {
            *(RRRRRR + l) = *(RRR + i);
            *(GGGGGG + l) = *(GGG + i);
            *(BBBBBB + l) = *(BBB + i);
            l++;
          }
        }
        MAXX = -99.9e64;
        MINN = 99.9e64;
        for (i = 0; i < l; i++) {
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
  for (i = 0; i < hsize * vsize; i++) {
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
      for (i = 0; i < hsize * vsize; i++) {
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
          for (i = 0; i < hsize * vsize; i++) {
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
          for (i = 0; i < hsize * vsize; i++) {
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
          for (i = 0; i < hsize * vsize; i++) {
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
          for (i = 0; i < hsize * vsize; i++) {
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
      for (i = 0; i < hsize * vsize; i++) {
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
      for (i = 0; i < 3; i++) {
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

        for (i = 0; i < hsize * vsize; i++) {
          if (*(INDEX + i) == PT[MEN][1].INDEXNO) {
            *(RR + i) = (double)(((*(YIN + i))) - U_R);
            *(GG + i) = (double)(((*(UIN + i))) - U_G);
            *(BB + i) = (double)(((*(VIN + i))) - U_B);
          }
        }
        V[0] = X1[0][Y];
        V[1] = X1[1][Y];
        V[2] = X1[2][Y];
        for (i = 0; i < hsize * vsize; i++) {
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
        for (i = 0, l = 0; i < hsize * vsize; i++) {
          if (*(INDEX + i) == PT[MEN][1].INDEXNO) {
            *(RRRRRR + l) = *(RRR + i);
            *(GGGGGG + l) = *(GGG + i);
            *(BBBBBB + l) = *(BBB + i);
            l++;
          }
        }
        MAXX = -99.9e64;
        MINN = 99.9e64;
        for (i = 0; i < l; i++) {
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

  if ((bun == 1) || (bun == 2) || (bun == 3) || (bun == 4) || (bun == 5)) {
    if (PT[MEN][0].MAXDISTANCE /**PT[MEN][0].INDEXNUM*/ >=
        PT[MEN][1]
            .MAXDISTANCE /**PT[MEN][1].INDEXNUM*/) { // maxdistance div tuika
      j = 0;                                         // INDEXNO
    } else {
      j = 1;
    }
  } else if (bun == 0) {
    if (PT[MEN][0].MAXDISTANCE * PT[MEN][0].INDEXNUM >=
        PT[MEN][1].MAXDISTANCE * PT[MEN][1].INDEXNUM) { // maxdistance div tuika
      j = 0;                                            // INDEXNO
    } else {
      j = 1;
    }
  }

  U_R = 0.0;
  U_G = 0.0;
  U_B = 0.0;
  for (i = 0; i < hsize * vsize; i++) {
    if ((*(INDEX + i)) == j) {
      U_R += *(YIN + i);
      U_G += *(UIN + i);
      U_B += *(VIN + i);
    }
  }

  U_R /= (double)PT[MEN][j].INDEXNUM;
  U_G /= (double)PT[MEN][j].INDEXNUM;
  U_B /= (double)PT[MEN][j].INDEXNUM;
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
  for (i = 0; i < hsize * vsize; i++) {
    if ((*(INDEX + i)) == j) {
      TMP_RR += (double)((((*(YIN + i))) - U_R) * (((*(YIN + i))) - U_R));
      TMP_GG += (double)((((*(UIN + i))) - U_G) * (((*(UIN + i))) - U_G));
      TMP_BB += (double)((((*(VIN + i))) - U_B) * (((*(VIN + i))) - U_B));
      TMP_RG += (double)((((*(YIN + i))) - U_R) * (((*(UIN + i))) - U_G));
      TMP_RB += (double)((((*(YIN + i))) - U_R) * (((*(VIN + i))) - U_B));
      TMP_GB += (double)((((*(UIN + i))) - U_G) * (((*(VIN + i))) - U_B));
    }
  }
  TMP_RR /= (double)(PT[MEN][j].INDEXNUM);
  TMP_GG /= (double)(PT[MEN][j].INDEXNUM);
  TMP_BB /= (double)(PT[MEN][j].INDEXNUM);
  TMP_RG /= (double)(PT[MEN][j].INDEXNUM);
  TMP_RB /= (double)(PT[MEN][j].INDEXNUM);
  TMP_GB /= (double)(PT[MEN][j].INDEXNUM);

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
  for (i = 0; i < 3; i++) {
    if (A1[i][i] > IGENMAX) {
      IGENMAX = A1[i][i];
      Y = i;
    }
  }
  DIVIDENUM = 3;

  // ohtsu tuika start
  for (i = 0; i < hsize * vsize; i++) {
    if ((*(INDEX + i)) == j) {
      *(RR + i) = (double)(((*(YIN + i))) - U_R);
      *(GG + i) = (double)(((*(UIN + i))) - U_G);
      *(BB + i) = (double)(((*(VIN + i))) - U_B);
    }
  }
  V[0] = X1[0][Y];
  V[1] = X1[1][Y];
  V[2] = X1[2][Y];
  for (i = 0; i < hsize * vsize; i++) {
    if ((*(INDEX + i)) == j) {
      kaiten(V, *(RR + i), *(GG + i), *(BB + i), RRR + i, GGG + i, BBB + i);
    }
  }
  GASOSUU = PT[MEN][j].INDEXNUM;
  for (i = 0; i < hsize * vsize; i++) {
    *(RRR33 + i) = (double)(*(RRR + i));
  }
  for (i = 0, k = 0; i < hsize * vsize; i++) {
    if ((*(INDEX + i)) == j) {
      *(RR + k) = (double)(*(RRR + i));
      *(GG + k) = (double)(*(GGG + i));
      *(BB + k) = (double)(*(BBB + i));
      k++;
    }
  }
  for (i = 0; i < hsize * vsize; i++) {
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

  for (i = 0, k = 0, l = 0; i < hsize * vsize; i++) {
    if ((*(INDEX + i)) == j) {
      // ZZ = X1[0][Y]*((double)(((*(YIN+i)))) - ((double)(U_R)+(*XXX))) +
      // X1[1][Y]*((double)(((*(UIN+i)))) - ((double)(U_G)+(*YYY))) +
      // X1[2][Y]*((double)(((*(VIN+i)))) - ((double)(U_B)+(*ZZZ))) ;
      // if(ZZ>=0.0){
      if (*(RRR33 + i) >= DTHRESH) {
        *(INDEX + i) = DIVIDENUM - 1;
        k++;
      } else {
        ;
        l++;
      }
    }
  }
  PT[MEN][0].INDEXNO = j;
  PT[MEN][1].INDEXNO = DIVIDENUM - 1;
  PT[MEN][0].INDEXNUM = l;
  PT[MEN][1].INDEXNUM = k;

  // max distance div tuika

  // maxdistance tuika start
  U_R = 0.0;
  U_G = 0.0;
  U_B = 0.0;
  // 0側のmaxdistanceを求める
  for (i = 0; i < hsize * vsize; i++) {
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
      for (i = 0; i < hsize * vsize; i++) {
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
          for (i = 0; i < hsize * vsize; i++) {
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
          for (i = 0; i < hsize * vsize; i++) {
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
          for (i = 0; i < hsize * vsize; i++) {
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
          for (i = 0; i < hsize * vsize; i++) {
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
      for (i = 0; i < hsize * vsize; i++) {
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
      for (i = 0; i < 3; i++) {
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

        for (i = 0; i < hsize * vsize; i++) {
          if (*(INDEX + i) == PT[MEN][0].INDEXNO) {
            *(RR + i) = (double)(((*(YIN + i))) - U_R);
            *(GG + i) = (double)(((*(UIN + i))) - U_G);
            *(BB + i) = (double)(((*(VIN + i))) - U_B);
          }
        }
        V[0] = X1[0][Y];
        V[1] = X1[1][Y];
        V[2] = X1[2][Y];
        for (i = 0; i < hsize * vsize; i++) {
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
        for (i = 0, l = 0; i < hsize * vsize; i++) {
          if (*(INDEX + i) == PT[MEN][0].INDEXNO) {
            *(RRRRRR + l) = *(RRR + i);
            *(GGGGGG + l) = *(GGG + i);
            *(BBBBBB + l) = *(BBB + i);
            l++;
          }
        }
        MAXX = -99.9e64;
        MINN = 99.9e64;
        for (i = 0; i < l; i++) {
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
  for (i = 0; i < hsize * vsize; i++) {
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
      for (i = 0; i < hsize * vsize; i++) {
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
          for (i = 0; i < hsize * vsize; i++) {
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
          for (i = 0; i < hsize * vsize; i++) {
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
          for (i = 0; i < hsize * vsize; i++) {
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
          for (i = 0; i < hsize * vsize; i++) {
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
      for (i = 0; i < hsize * vsize; i++) {
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
      for (i = 0; i < 3; i++) {
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
        for (i = 0; i < hsize * vsize; i++) {
          if (*(INDEX + i) == PT[MEN][1].INDEXNO) {
            *(RR + i) = (double)(((*(YIN + i))) - U_R);
            *(GG + i) = (double)(((*(UIN + i))) - U_G);
            *(BB + i) = (double)(((*(VIN + i))) - U_B);
          }
        }
        V[0] = X1[0][Y];
        V[1] = X1[1][Y];
        V[2] = X1[2][Y];
        for (i = 0; i < hsize * vsize; i++) {
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
        for (i = 0, l = 0; i < hsize * vsize; i++) {
          if (*(INDEX + i) == PT[MEN][1].INDEXNO) {
            *(RRRRRR + l) = *(RRR + i);
            *(GGGGGG + l) = *(GGG + i);
            *(BBBBBB + l) = *(BBB + i);
            l++;
          }
        }
        MAXX = -99.9e64;
        MINN = 99.9e64;
        for (i = 0; i < l; i++) {
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

  if (j == 0) {
    for (i = 2; i <= DIVIDENUM - 1; i++) {
      PT[MEN][i].INDEXNO = PT[NMEN][i - 1].INDEXNO;
      PT[MEN][i].INDEXNUM = PT[NMEN][i - 1].INDEXNUM;
      PT[MEN][i].MAXDISTANCE =
          PT[NMEN][i - 1].MAXDISTANCE; // max distance div tuika
    }
  }
  if ((j >= 1) && (j <= (DIVIDENUM - 3))) {
    for (i = 2; i <= j + 1; i++) {
      //// i i-2
      PT[MEN][i].INDEXNO = PT[NMEN][i - 2].INDEXNO;
      PT[MEN][i].INDEXNUM = PT[NMEN][i - 2].INDEXNUM;
      PT[MEN][i].MAXDISTANCE =
          PT[NMEN][i - 2].MAXDISTANCE; // max distance div tuika
    }

    for (i = j + 2; i <= DIVIDENUM - 1; i++) {
      //// i i-1
      PT[MEN][i].INDEXNO = PT[NMEN][i - 1].INDEXNO;
      PT[MEN][i].INDEXNUM = PT[NMEN][i - 1].INDEXNUM;
      PT[MEN][i].MAXDISTANCE =
          PT[NMEN][i - 1].MAXDISTANCE; // max distance div tuika
    }
  }
  if (j == (DIVIDENUM - 2)) {
    for (i = 2; i <= j + 1; i++) {
      //// i i-2
      PT[MEN][i].INDEXNO = PT[NMEN][i - 2].INDEXNO;
      PT[MEN][i].INDEXNUM = PT[NMEN][i - 2].INDEXNUM;
      PT[MEN][i].MAXDISTANCE =
          PT[NMEN][i - 2].MAXDISTANCE; // max distance div tuika
    }
  }
  double MAXBUNSAN = 100.0;
  // 9/9 kokokara
  while (DIVIDENUM < 256) {
    //最大画素数のブロックをさがし、そのINDEXNOを調べる。
    MAXINDEXNUM = -1.0;

    if ((bun == 1) || (bun == 2) || (bun == 3) || (bun == 4) || (bun == 5)) {
      for (i = 0; i < DIVIDENUM; i++) {
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
      for (i = 0; i < DIVIDENUM; i++) {
        // if(PT[MEN][i].INDEXNUM > MAXINDEXNUM)
        // MAXINDEXNUM = PT[MEN][i].INDEXNUM;
        if (PT[MEN][i].MAXDISTANCE * PT[MEN][i].INDEXNUM > MAXINDEXNUM) {
          MAXINDEXNUM = PT[MEN][i].MAXDISTANCE *
                        PT[MEN][i].INDEXNUM; // max distance tuika
          NUM = i;
        }
      }
    }

    FLAG = 0;
  label:
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
    for (i = 0; i < hsize * vsize; i++) {
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
    for (i = 0; i < hsize * vsize; i++) {
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
      for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
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
    for (i = 0; i < 3; i++) {
      if (A1[i][i] > IGENMAX) {
        IGENMAX = A1[i][i];
        Y = i;
      }
    }
    DIVIDENUM++;

    // ohtsu tuika start
    for (i = 0; i < hsize * vsize; i++) {
      if ((*(INDEX + i)) == PT[MEN][NUM].INDEXNO) {
        *(RR + i) = (double)(((*(YIN + i))) - U_R);
        *(GG + i) = (double)(((*(UIN + i))) - U_G);
        *(BB + i) = (double)(((*(VIN + i))) - U_B);
      }
    }
    V[0] = X1[0][Y];
    V[1] = X1[1][Y];
    V[2] = X1[2][Y];
    for (i = 0; i < hsize * vsize; i++) {
      if ((*(INDEX + i)) == PT[MEN][NUM].INDEXNO) {
        kaiten(V, *(RR + i), *(GG + i), *(BB + i), RRR + i, GGG + i, BBB + i);
      }
    }
    GASOSUU = PT[MEN][NUM].INDEXNUM;
    for (i = 0; i < hsize * vsize; i++) {
      *(RRR33 + i) = (double)(*(RRR + i));
    }
    for (i = 0, k = 0; i < hsize * vsize; i++) {
      if ((*(INDEX + i)) == PT[MEN][NUM].INDEXNO) {
        *(RR + k) = (double)(*(RRR + i));
        *(GG + k) = (double)(*(GGG + i));
        *(BB + k) = (double)(*(BBB + i));
        k++;
      }
    }
    for (i = 0; i < hsize * vsize; i++) {
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

    for (i = 0, k = 0, l = 0; i < hsize * vsize; i++) {
      if ((*(INDEX + i)) == PT[MEN][NUM].INDEXNO) {
        // ZZ = X1[0][Y]*((double)(((*(YIN+i)))) - ((double)(U_R)+(*XXX))) +
        // X1[1][Y]*((double)(((*(UIN+i)))) - ((double)(U_G)+(*YYY))) +
        // X1[2][Y]*((double)(((*(VIN+i)))) - ((double)(U_B)+(*ZZZ))) ;
        // if(ZZ>=0.0){
        if (*(RRR33 + i) >= DTHRESH) {
          *(INDEX + i) = DIVIDENUM - 1;
          k++;
        } else {
          ;
          l++;
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
    PT[MEN][0].INDEXNUM = l;
    PT[MEN][1].INDEXNUM = k;
    // debug start
    if (DIVIDENUM == 5 || DIVIDENUM == 6) {
      fprintf(stderr, "INDEXNO%d %d\n", PT[MEN][0].INDEXNO, l);
      fprintf(stderr, "INDEXNO%d %d\n", PT[MEN][1].INDEXNO, k);
    }
    // while(1);
    // debug end

    // max distance divide tuika
    // maxdistance tuika start
    U_R = 0.0;
    U_G = 0.0;
    U_B = 0.0;
    // 0側のmaxdistanceを求める
    for (i = 0; i < hsize * vsize; i++) {
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
        for (i = 0; i < hsize * vsize; i++) {
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
            for (i = 0; i < hsize * vsize; i++) {
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
            for (i = 0; i < hsize * vsize; i++) {
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
            for (i = 0; i < hsize * vsize; i++) {
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
            for (i = 0; i < hsize * vsize; i++) {
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
        for (i = 0; i < hsize * vsize; i++) {
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
        for (i = 0; i < 3; i++) {
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
          for (i = 0; i < hsize * vsize; i++) {
            if (*(INDEX + i) == PT[MEN][0].INDEXNO) {
              *(RR + i) = (double)(((*(YIN + i))) - U_R);
              *(GG + i) = (double)(((*(UIN + i))) - U_G);
              *(BB + i) = (double)(((*(VIN + i))) - U_B);
            }
          }
          V[0] = X1[0][Y];
          V[1] = X1[1][Y];
          V[2] = X1[2][Y];
          for (i = 0; i < hsize * vsize; i++) {
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
          for (i = 0, l = 0; i < hsize * vsize; i++) {
            if (*(INDEX + i) == PT[MEN][0].INDEXNO) {
              *(RRRRRR + l) = *(RRR + i);
              *(GGGGGG + l) = *(GGG + i);
              *(BBBBBB + l) = *(BBB + i);
              l++;
            }
          }
          MAXX = -99.9e64;
          MINN = 99.9e64;
          for (i = 0; i < l; i++) {
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
    for (i = 0; i < hsize * vsize; i++) {
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
        for (i = 0; i < hsize * vsize; i++) {
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
            for (i = 0; i < hsize * vsize; i++) {
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
            for (i = 0; i < hsize * vsize; i++) {
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
            for (i = 0; i < hsize * vsize; i++) {
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
            for (i = 0; i < hsize * vsize; i++) {
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
        for (i = 0; i < hsize * vsize; i++) {
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
        for (i = 0; i < 3; i++) {
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
          for (i = 0; i < hsize * vsize; i++) {
            if (*(INDEX + i) == PT[MEN][1].INDEXNO) {
              *(RR + i) = (double)(((*(YIN + i))) - U_R);
              *(GG + i) = (double)(((*(UIN + i))) - U_G);
              *(BB + i) = (double)(((*(VIN + i))) - U_B);
            }
          }
          V[0] = X1[0][Y];
          V[1] = X1[1][Y];
          V[2] = X1[2][Y];
          for (i = 0; i < hsize * vsize; i++) {
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
          for (i = 0, l = 0; i < hsize * vsize; i++) {
            if (*(INDEX + i) == PT[MEN][1].INDEXNO) {
              *(RRRRRR + l) = *(RRR + i);
              *(GGGGGG + l) = *(GGG + i);
              *(BBBBBB + l) = *(BBB + i);
              l++;
            }
          }
          MAXX = -99.9e64;
          MINN = 99.9e64;
          for (i = 0; i < l; i++) {
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
      for (i = 2; i <= DIVIDENUM - 1; i++) {
        PT[MEN][i].INDEXNO = PT[NMEN][i - 1].INDEXNO;
        PT[MEN][i].INDEXNUM = PT[NMEN][i - 1].INDEXNUM;
        PT[MEN][i].MAXDISTANCE = PT[NMEN][i - 1].MAXDISTANCE;
      }
    }
    if ((NUM >= 1) && (NUM <= (DIVIDENUM - 3))) {
      for (i = 2; i <= NUM + 1; i++) {
        // i i-2
        PT[MEN][i].INDEXNO = PT[NMEN][i - 2].INDEXNO;
        PT[MEN][i].INDEXNUM = PT[NMEN][i - 2].INDEXNUM;
        PT[MEN][i].MAXDISTANCE = PT[NMEN][i - 2].MAXDISTANCE;
      }
      for (i = NUM + 2; i <= DIVIDENUM - 1; i++) {
        // i i-1
        PT[MEN][i].INDEXNO = PT[NMEN][i - 1].INDEXNO;
        PT[MEN][i].INDEXNUM = PT[NMEN][i - 1].INDEXNUM;
        PT[MEN][i].MAXDISTANCE = PT[NMEN][i - 1].MAXDISTANCE;
      }
    }
    if (NUM == (DIVIDENUM - 2)) {
      for (i = 2; i <= NUM + 1; i++) {
        // i i-2
        PT[MEN][i].INDEXNO = PT[NMEN][i - 2].INDEXNO;
        PT[MEN][i].INDEXNUM = PT[NMEN][i - 2].INDEXNUM;
        PT[MEN][i].MAXDISTANCE = PT[NMEN][i - 2].MAXDISTANCE;
      }
    }
    // debug start
    if (DIVIDENUM == 256) {
      for (i = 0; i < DIVIDENUM; i++) {
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
  for (i = 0; i < IROSUU; i++) {
    Y_JYUSHIN[i] = 0.0;
    U_JYUSHIN[i] = 0.0;
    V_JYUSHIN[i] = 0.0;
  }
  int w;
  w = 0;
  for (i = 0; i < hsize * vsize; i++) {

    if (index3[i] != -1) {
      for (j = 0; j < IROSUU; j++) {
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
  for (i = 0; i < IROSUU; i++) {
    y += PT[MEN][i].INDEXNUM;
  }
  printf("palet sou gasosuu=%d\n", y);

  for (m = 0; m < IROSUU; m++) {
    k = PT[MEN][m].INDEXNO;
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
  int s, q;
  p = 0;
  s = 0;
  if (kmon == 1) {
    while (1) {
      p++;
      fprintf(stderr, "%d iter\n", p);
      for (i = 0; i < hsize * vsize; i++) {
        if (index3[i] != -1) {
          MINZ = 60e60;
          for (j = 0; j < IROSUU; j++) {
            Z[j] = ((((*(YIN + i)))) - Y_JYUSHIN[j]) *
                       ((((*(YIN + i)))) - Y_JYUSHIN[j]) +
                   ((((*(UIN + i)))) - U_JYUSHIN[j]) *
                       ((((*(UIN + i)))) - U_JYUSHIN[j]) +
                   ((((*(VIN + i)))) - V_JYUSHIN[j]) *
                       ((((*(VIN + i)))) - V_JYUSHIN[j]);
          }
          for (j = 0; j < IROSUU; j++) {
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
        for (i = 0; i < IROSUU; i++) {
          KARI[i] = 0;
        }
        for (i = 0; i < hsize * vsize; i++) {
          KARI[*(PALETGAZOU + i)]++;
        }
        for (i = 0; i < IROSUU; i++) {
          fprintf(stderr, "i=%d %d\n", i, KARI[i]);
        }
      }
      // debug end

      // kmeans法
      l = 0;
      for (i = 0; i < IROSUU; i++) {
        SUM_R[i] = 0.0;
        SUM_G[i] = 0.0;
        SUM_B[i] = 0.0;
      }
      for (j = 0; j < IROSUU; j++) {
        for (i = 0, k = 0; i < hsize * vsize; i++) {
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
          l++;
          HEIKIN_R[j] = (double)((rand()) % ((int)(gmult * 256)));
          HEIKIN_G[j] = (double)((rand()) % ((int)(bmult * 256)));
          HEIKIN_B[j] = (double)((rand()) % ((int)(rmult * 256)));
          // debug start
          printf("RAND de dasita palet = %f %f %f\n", HEIKIN_R[j], HEIKIN_G[j],
                 HEIKIN_B[j]);
          // debug end
        }
      }
      if (l == 0) {
        if (s == 0) {
          q = p + 50;
          s = 1;
        }
      }
      if (((int)HEIKIN_R[0] == (int)Y_JYUSHIN[0]) &&
          ((int)HEIKIN_R[1] == (int)Y_JYUSHIN[1]) &&
          ((int)HEIKIN_R[2] == (int)Y_JYUSHIN[2]) &&
          ((int)HEIKIN_R[3] == (int)Y_JYUSHIN[3]) &&
          ((int)HEIKIN_R[4] == (int)Y_JYUSHIN[4]) &&
          ((int)HEIKIN_R[5] == (int)Y_JYUSHIN[5]) &&
          ((int)HEIKIN_R[6] == (int)Y_JYUSHIN[6]) &&
          ((int)HEIKIN_R[7] == (int)Y_JYUSHIN[7]) &&
          ((int)HEIKIN_R[8] == (int)Y_JYUSHIN[8]) &&
          ((int)HEIKIN_R[9] == (int)Y_JYUSHIN[9]) &&
          ((int)HEIKIN_R[10] == (int)Y_JYUSHIN[10]) &&
          ((int)HEIKIN_R[11] == (int)Y_JYUSHIN[11]) &&
          ((int)HEIKIN_R[12] == (int)Y_JYUSHIN[12]) &&
          ((int)HEIKIN_R[13] == (int)Y_JYUSHIN[13]) &&
          ((int)HEIKIN_R[14] == (int)Y_JYUSHIN[14]) &&
          ((int)HEIKIN_R[15] == (int)Y_JYUSHIN[15]) &&
          ((int)HEIKIN_R[16] == (int)Y_JYUSHIN[16]) &&
          ((int)HEIKIN_R[17] == (int)Y_JYUSHIN[17]) &&
          ((int)HEIKIN_R[18] == (int)Y_JYUSHIN[18]) &&
          ((int)HEIKIN_R[19] == (int)Y_JYUSHIN[19]) &&
          ((int)HEIKIN_R[20] == (int)Y_JYUSHIN[20]) &&
          ((int)HEIKIN_R[21] == (int)Y_JYUSHIN[21]) &&
          ((int)HEIKIN_R[22] == (int)Y_JYUSHIN[22]) &&
          ((int)HEIKIN_R[23] == (int)Y_JYUSHIN[23]) &&
          ((int)HEIKIN_R[24] == (int)Y_JYUSHIN[24]) &&
          ((int)HEIKIN_R[25] == (int)Y_JYUSHIN[25]) &&
          ((int)HEIKIN_R[26] == (int)Y_JYUSHIN[26]) &&
          ((int)HEIKIN_R[27] == (int)Y_JYUSHIN[27]) &&
          ((int)HEIKIN_R[28] == (int)Y_JYUSHIN[28]) &&
          ((int)HEIKIN_R[29] == (int)Y_JYUSHIN[29]) &&
          ((int)HEIKIN_R[30] == (int)Y_JYUSHIN[30]) &&
          ((int)HEIKIN_R[31] == (int)Y_JYUSHIN[31]) &&
          ((int)HEIKIN_R[32] == (int)Y_JYUSHIN[32]) &&
          ((int)HEIKIN_R[33] == (int)Y_JYUSHIN[33]) &&
          ((int)HEIKIN_R[34] == (int)Y_JYUSHIN[34]) &&
          ((int)HEIKIN_R[35] == (int)Y_JYUSHIN[35]) &&
          ((int)HEIKIN_R[36] == (int)Y_JYUSHIN[36]) &&
          ((int)HEIKIN_R[37] == (int)Y_JYUSHIN[37]) &&
          ((int)HEIKIN_R[38] == (int)Y_JYUSHIN[38]) &&
          ((int)HEIKIN_R[39] == (int)Y_JYUSHIN[39]) &&
          ((int)HEIKIN_R[40] == (int)Y_JYUSHIN[40]) &&
          ((int)HEIKIN_R[41] == (int)Y_JYUSHIN[41]) &&
          ((int)HEIKIN_R[42] == (int)Y_JYUSHIN[42]) &&
          ((int)HEIKIN_R[43] == (int)Y_JYUSHIN[43]) &&
          ((int)HEIKIN_R[44] == (int)Y_JYUSHIN[44]) &&
          ((int)HEIKIN_R[45] == (int)Y_JYUSHIN[45]) &&
          ((int)HEIKIN_R[46] == (int)Y_JYUSHIN[46]) &&
          ((int)HEIKIN_R[47] == (int)Y_JYUSHIN[47]) &&
          ((int)HEIKIN_R[48] == (int)Y_JYUSHIN[48]) &&
          ((int)HEIKIN_R[49] == (int)Y_JYUSHIN[49]) &&
          ((int)HEIKIN_R[50] == (int)Y_JYUSHIN[50]) &&
          ((int)HEIKIN_R[51] == (int)Y_JYUSHIN[51]) &&
          ((int)HEIKIN_R[52] == (int)Y_JYUSHIN[52]) &&
          ((int)HEIKIN_R[53] == (int)Y_JYUSHIN[53]) &&
          ((int)HEIKIN_R[54] == (int)Y_JYUSHIN[54]) &&
          ((int)HEIKIN_R[55] == (int)Y_JYUSHIN[55]) &&
          ((int)HEIKIN_R[56] == (int)Y_JYUSHIN[56]) &&
          ((int)HEIKIN_R[57] == (int)Y_JYUSHIN[57]) &&
          ((int)HEIKIN_R[58] == (int)Y_JYUSHIN[58]) &&
          ((int)HEIKIN_R[59] == (int)Y_JYUSHIN[59]) &&
          ((int)HEIKIN_R[60] == (int)Y_JYUSHIN[60]) &&
          ((int)HEIKIN_R[61] == (int)Y_JYUSHIN[61]) &&
          ((int)HEIKIN_R[62] == (int)Y_JYUSHIN[62]) &&
          ((int)HEIKIN_R[63] == (int)Y_JYUSHIN[63]) &&
          ((int)HEIKIN_R[64] == (int)Y_JYUSHIN[64]) &&
          ((int)HEIKIN_R[65] == (int)Y_JYUSHIN[65]) &&
          ((int)HEIKIN_R[66] == (int)Y_JYUSHIN[66]) &&
          ((int)HEIKIN_R[67] == (int)Y_JYUSHIN[67]) &&
          ((int)HEIKIN_R[68] == (int)Y_JYUSHIN[68]) &&
          ((int)HEIKIN_R[69] == (int)Y_JYUSHIN[69]) &&
          ((int)HEIKIN_R[70] == (int)Y_JYUSHIN[70]) &&
          ((int)HEIKIN_R[71] == (int)Y_JYUSHIN[71]) &&
          ((int)HEIKIN_R[72] == (int)Y_JYUSHIN[72]) &&
          ((int)HEIKIN_R[73] == (int)Y_JYUSHIN[73]) &&
          ((int)HEIKIN_R[74] == (int)Y_JYUSHIN[74]) &&
          ((int)HEIKIN_R[75] == (int)Y_JYUSHIN[75]) &&
          ((int)HEIKIN_R[76] == (int)Y_JYUSHIN[76]) &&
          ((int)HEIKIN_R[77] == (int)Y_JYUSHIN[77]) &&
          ((int)HEIKIN_R[78] == (int)Y_JYUSHIN[78]) &&
          ((int)HEIKIN_R[79] == (int)Y_JYUSHIN[79]) &&
          ((int)HEIKIN_R[80] == (int)Y_JYUSHIN[80]) &&
          ((int)HEIKIN_R[81] == (int)Y_JYUSHIN[81]) &&
          ((int)HEIKIN_R[82] == (int)Y_JYUSHIN[82]) &&
          ((int)HEIKIN_R[83] == (int)Y_JYUSHIN[83]) &&
          ((int)HEIKIN_R[84] == (int)Y_JYUSHIN[84]) &&
          ((int)HEIKIN_R[85] == (int)Y_JYUSHIN[85]) &&
          ((int)HEIKIN_R[86] == (int)Y_JYUSHIN[86]) &&
          ((int)HEIKIN_R[87] == (int)Y_JYUSHIN[87]) &&
          ((int)HEIKIN_R[88] == (int)Y_JYUSHIN[88]) &&
          ((int)HEIKIN_R[89] == (int)Y_JYUSHIN[89]) &&
          ((int)HEIKIN_R[90] == (int)Y_JYUSHIN[90]) &&
          ((int)HEIKIN_R[91] == (int)Y_JYUSHIN[91]) &&
          ((int)HEIKIN_R[92] == (int)Y_JYUSHIN[92]) &&
          ((int)HEIKIN_R[93] == (int)Y_JYUSHIN[93]) &&
          ((int)HEIKIN_R[94] == (int)Y_JYUSHIN[94]) &&
          ((int)HEIKIN_R[95] == (int)Y_JYUSHIN[95]) &&
          ((int)HEIKIN_R[96] == (int)Y_JYUSHIN[96]) &&
          ((int)HEIKIN_R[97] == (int)Y_JYUSHIN[97]) &&
          ((int)HEIKIN_R[98] == (int)Y_JYUSHIN[98]) &&
          ((int)HEIKIN_R[99] == (int)Y_JYUSHIN[99]) &&
          ((int)HEIKIN_R[100] == (int)Y_JYUSHIN[100]) &&
          ((int)HEIKIN_R[101] == (int)Y_JYUSHIN[101]) &&
          ((int)HEIKIN_R[102] == (int)Y_JYUSHIN[102]) &&
          ((int)HEIKIN_R[103] == (int)Y_JYUSHIN[103]) &&
          ((int)HEIKIN_R[104] == (int)Y_JYUSHIN[104]) &&
          ((int)HEIKIN_R[105] == (int)Y_JYUSHIN[105]) &&
          ((int)HEIKIN_R[106] == (int)Y_JYUSHIN[106]) &&
          ((int)HEIKIN_R[107] == (int)Y_JYUSHIN[107]) &&
          ((int)HEIKIN_R[108] == (int)Y_JYUSHIN[108]) &&
          ((int)HEIKIN_R[109] == (int)Y_JYUSHIN[109]) &&
          ((int)HEIKIN_R[110] == (int)Y_JYUSHIN[110]) &&
          ((int)HEIKIN_R[111] == (int)Y_JYUSHIN[111]) &&
          ((int)HEIKIN_R[112] == (int)Y_JYUSHIN[112]) &&
          ((int)HEIKIN_R[113] == (int)Y_JYUSHIN[113]) &&
          ((int)HEIKIN_R[114] == (int)Y_JYUSHIN[114]) &&
          ((int)HEIKIN_R[115] == (int)Y_JYUSHIN[115]) &&
          ((int)HEIKIN_R[116] == (int)Y_JYUSHIN[116]) &&
          ((int)HEIKIN_R[117] == (int)Y_JYUSHIN[117]) &&
          ((int)HEIKIN_R[118] == (int)Y_JYUSHIN[118]) &&
          ((int)HEIKIN_R[119] == (int)Y_JYUSHIN[119]) &&
          ((int)HEIKIN_R[120] == (int)Y_JYUSHIN[120]) &&
          ((int)HEIKIN_R[121] == (int)Y_JYUSHIN[121]) &&
          ((int)HEIKIN_R[122] == (int)Y_JYUSHIN[122]) &&
          ((int)HEIKIN_R[123] == (int)Y_JYUSHIN[123]) &&
          ((int)HEIKIN_R[124] == (int)Y_JYUSHIN[124]) &&
          ((int)HEIKIN_R[125] == (int)Y_JYUSHIN[125]) &&
          ((int)HEIKIN_R[126] == (int)Y_JYUSHIN[126]) &&
          ((int)HEIKIN_R[127] == (int)Y_JYUSHIN[127]) &&
          ((int)HEIKIN_R[128] == (int)Y_JYUSHIN[128]) &&
          ((int)HEIKIN_R[129] == (int)Y_JYUSHIN[129]) &&
          ((int)HEIKIN_R[130] == (int)Y_JYUSHIN[130]) &&
          ((int)HEIKIN_R[131] == (int)Y_JYUSHIN[131]) &&
          ((int)HEIKIN_R[132] == (int)Y_JYUSHIN[132]) &&
          ((int)HEIKIN_R[133] == (int)Y_JYUSHIN[133]) &&
          ((int)HEIKIN_R[134] == (int)Y_JYUSHIN[134]) &&
          ((int)HEIKIN_R[135] == (int)Y_JYUSHIN[135]) &&
          ((int)HEIKIN_R[136] == (int)Y_JYUSHIN[136]) &&
          ((int)HEIKIN_R[137] == (int)Y_JYUSHIN[137]) &&
          ((int)HEIKIN_R[138] == (int)Y_JYUSHIN[138]) &&
          ((int)HEIKIN_R[139] == (int)Y_JYUSHIN[139]) &&
          ((int)HEIKIN_R[140] == (int)Y_JYUSHIN[140]) &&
          ((int)HEIKIN_R[141] == (int)Y_JYUSHIN[141]) &&
          ((int)HEIKIN_R[142] == (int)Y_JYUSHIN[142]) &&
          ((int)HEIKIN_R[143] == (int)Y_JYUSHIN[143]) &&
          ((int)HEIKIN_R[144] == (int)Y_JYUSHIN[144]) &&
          ((int)HEIKIN_R[145] == (int)Y_JYUSHIN[145]) &&
          ((int)HEIKIN_R[146] == (int)Y_JYUSHIN[146]) &&
          ((int)HEIKIN_R[147] == (int)Y_JYUSHIN[147]) &&
          ((int)HEIKIN_R[148] == (int)Y_JYUSHIN[148]) &&
          ((int)HEIKIN_R[149] == (int)Y_JYUSHIN[149]) &&
          ((int)HEIKIN_R[150] == (int)Y_JYUSHIN[150]) &&
          ((int)HEIKIN_R[151] == (int)Y_JYUSHIN[151]) &&
          ((int)HEIKIN_R[152] == (int)Y_JYUSHIN[152]) &&
          ((int)HEIKIN_R[153] == (int)Y_JYUSHIN[153]) &&
          ((int)HEIKIN_R[154] == (int)Y_JYUSHIN[154]) &&
          ((int)HEIKIN_R[155] == (int)Y_JYUSHIN[155]) &&
          ((int)HEIKIN_R[156] == (int)Y_JYUSHIN[156]) &&
          ((int)HEIKIN_R[157] == (int)Y_JYUSHIN[157]) &&
          ((int)HEIKIN_R[158] == (int)Y_JYUSHIN[158]) &&
          ((int)HEIKIN_R[159] == (int)Y_JYUSHIN[159]) &&
          ((int)HEIKIN_R[160] == (int)Y_JYUSHIN[160]) &&
          ((int)HEIKIN_R[161] == (int)Y_JYUSHIN[161]) &&
          ((int)HEIKIN_R[162] == (int)Y_JYUSHIN[162]) &&
          ((int)HEIKIN_R[163] == (int)Y_JYUSHIN[163]) &&
          ((int)HEIKIN_R[164] == (int)Y_JYUSHIN[164]) &&
          ((int)HEIKIN_R[165] == (int)Y_JYUSHIN[165]) &&
          ((int)HEIKIN_R[166] == (int)Y_JYUSHIN[166]) &&
          ((int)HEIKIN_R[167] == (int)Y_JYUSHIN[167]) &&
          ((int)HEIKIN_R[168] == (int)Y_JYUSHIN[168]) &&
          ((int)HEIKIN_R[169] == (int)Y_JYUSHIN[169]) &&
          ((int)HEIKIN_R[170] == (int)Y_JYUSHIN[170]) &&
          ((int)HEIKIN_R[171] == (int)Y_JYUSHIN[171]) &&
          ((int)HEIKIN_R[172] == (int)Y_JYUSHIN[172]) &&
          ((int)HEIKIN_R[173] == (int)Y_JYUSHIN[173]) &&
          ((int)HEIKIN_R[174] == (int)Y_JYUSHIN[174]) &&
          ((int)HEIKIN_R[175] == (int)Y_JYUSHIN[175]) &&
          ((int)HEIKIN_R[176] == (int)Y_JYUSHIN[176]) &&
          ((int)HEIKIN_R[177] == (int)Y_JYUSHIN[177]) &&
          ((int)HEIKIN_R[178] == (int)Y_JYUSHIN[178]) &&
          ((int)HEIKIN_R[179] == (int)Y_JYUSHIN[179]) &&
          ((int)HEIKIN_R[180] == (int)Y_JYUSHIN[180]) &&
          ((int)HEIKIN_R[181] == (int)Y_JYUSHIN[181]) &&
          ((int)HEIKIN_R[182] == (int)Y_JYUSHIN[182]) &&
          ((int)HEIKIN_R[183] == (int)Y_JYUSHIN[183]) &&
          ((int)HEIKIN_R[184] == (int)Y_JYUSHIN[184]) &&
          ((int)HEIKIN_R[185] == (int)Y_JYUSHIN[185]) &&
          ((int)HEIKIN_R[186] == (int)Y_JYUSHIN[186]) &&
          ((int)HEIKIN_R[187] == (int)Y_JYUSHIN[187]) &&
          ((int)HEIKIN_R[188] == (int)Y_JYUSHIN[188]) &&
          ((int)HEIKIN_R[189] == (int)Y_JYUSHIN[189]) &&
          ((int)HEIKIN_R[190] == (int)Y_JYUSHIN[190]) &&
          ((int)HEIKIN_R[191] == (int)Y_JYUSHIN[191]) &&
          ((int)HEIKIN_R[192] == (int)Y_JYUSHIN[192]) &&
          ((int)HEIKIN_R[193] == (int)Y_JYUSHIN[193]) &&
          ((int)HEIKIN_R[194] == (int)Y_JYUSHIN[194]) &&
          ((int)HEIKIN_R[195] == (int)Y_JYUSHIN[195]) &&
          ((int)HEIKIN_R[196] == (int)Y_JYUSHIN[196]) &&
          ((int)HEIKIN_R[197] == (int)Y_JYUSHIN[197]) &&
          ((int)HEIKIN_R[198] == (int)Y_JYUSHIN[198]) &&
          ((int)HEIKIN_R[199] == (int)Y_JYUSHIN[199]) &&
          ((int)HEIKIN_R[200] == (int)Y_JYUSHIN[200]) &&
          ((int)HEIKIN_R[201] == (int)Y_JYUSHIN[201]) &&
          ((int)HEIKIN_R[202] == (int)Y_JYUSHIN[202]) &&
          ((int)HEIKIN_R[203] == (int)Y_JYUSHIN[203]) &&
          ((int)HEIKIN_R[204] == (int)Y_JYUSHIN[204]) &&
          ((int)HEIKIN_R[205] == (int)Y_JYUSHIN[205]) &&
          ((int)HEIKIN_R[206] == (int)Y_JYUSHIN[206]) &&
          ((int)HEIKIN_R[207] == (int)Y_JYUSHIN[207]) &&
          ((int)HEIKIN_R[208] == (int)Y_JYUSHIN[208]) &&
          ((int)HEIKIN_R[209] == (int)Y_JYUSHIN[209]) &&
          ((int)HEIKIN_R[210] == (int)Y_JYUSHIN[210]) &&
          ((int)HEIKIN_R[211] == (int)Y_JYUSHIN[211]) &&
          ((int)HEIKIN_R[212] == (int)Y_JYUSHIN[212]) &&
          ((int)HEIKIN_R[213] == (int)Y_JYUSHIN[213]) &&
          ((int)HEIKIN_R[214] == (int)Y_JYUSHIN[214]) &&
          ((int)HEIKIN_R[215] == (int)Y_JYUSHIN[215]) &&
          ((int)HEIKIN_R[216] == (int)Y_JYUSHIN[216]) &&
          ((int)HEIKIN_R[217] == (int)Y_JYUSHIN[217]) &&
          ((int)HEIKIN_R[218] == (int)Y_JYUSHIN[218]) &&
          ((int)HEIKIN_R[219] == (int)Y_JYUSHIN[219]) &&
          ((int)HEIKIN_R[220] == (int)Y_JYUSHIN[220]) &&
          ((int)HEIKIN_R[221] == (int)Y_JYUSHIN[221]) &&
          ((int)HEIKIN_R[222] == (int)Y_JYUSHIN[222]) &&
          ((int)HEIKIN_R[223] == (int)Y_JYUSHIN[223]) &&
          ((int)HEIKIN_R[224] == (int)Y_JYUSHIN[224]) &&
          ((int)HEIKIN_R[225] == (int)Y_JYUSHIN[225]) &&
          ((int)HEIKIN_R[226] == (int)Y_JYUSHIN[226]) &&
          ((int)HEIKIN_R[227] == (int)Y_JYUSHIN[227]) &&
          ((int)HEIKIN_R[228] == (int)Y_JYUSHIN[228]) &&
          ((int)HEIKIN_R[229] == (int)Y_JYUSHIN[229]) &&
          ((int)HEIKIN_R[230] == (int)Y_JYUSHIN[230]) &&
          ((int)HEIKIN_R[231] == (int)Y_JYUSHIN[231]) &&
          ((int)HEIKIN_R[232] == (int)Y_JYUSHIN[232]) &&
          ((int)HEIKIN_R[233] == (int)Y_JYUSHIN[233]) &&
          ((int)HEIKIN_R[234] == (int)Y_JYUSHIN[234]) &&
          ((int)HEIKIN_R[235] == (int)Y_JYUSHIN[235]) &&
          ((int)HEIKIN_R[236] == (int)Y_JYUSHIN[236]) &&
          ((int)HEIKIN_R[237] == (int)Y_JYUSHIN[237]) &&
          ((int)HEIKIN_R[238] == (int)Y_JYUSHIN[238]) &&
          ((int)HEIKIN_R[239] == (int)Y_JYUSHIN[239]) &&
          ((int)HEIKIN_R[240] == (int)Y_JYUSHIN[240]) &&
          ((int)HEIKIN_R[241] == (int)Y_JYUSHIN[241]) &&
          ((int)HEIKIN_R[242] == (int)Y_JYUSHIN[242]) &&
          ((int)HEIKIN_R[243] == (int)Y_JYUSHIN[243]) &&
          ((int)HEIKIN_R[244] == (int)Y_JYUSHIN[244]) &&
          ((int)HEIKIN_R[245] == (int)Y_JYUSHIN[245]) &&
          ((int)HEIKIN_R[246] == (int)Y_JYUSHIN[246]) &&
          ((int)HEIKIN_R[247] == (int)Y_JYUSHIN[247]) &&
          ((int)HEIKIN_R[248] == (int)Y_JYUSHIN[248]) &&
          ((int)HEIKIN_R[249] == (int)Y_JYUSHIN[249]) &&
          ((int)HEIKIN_R[250] == (int)Y_JYUSHIN[250]) &&
          ((int)HEIKIN_R[251] == (int)Y_JYUSHIN[251]) &&
          ((int)HEIKIN_R[252] == (int)Y_JYUSHIN[252]) &&
          ((int)HEIKIN_R[253] == (int)Y_JYUSHIN[253]) &&
          ((int)HEIKIN_R[254] == (int)Y_JYUSHIN[254]) &&
          ((int)HEIKIN_R[255] == (int)Y_JYUSHIN[255]) &&
          ((int)HEIKIN_G[0] == (int)U_JYUSHIN[0]) &&
          ((int)HEIKIN_G[1] == (int)U_JYUSHIN[1]) &&
          ((int)HEIKIN_G[2] == (int)U_JYUSHIN[2]) &&
          ((int)HEIKIN_G[3] == (int)U_JYUSHIN[3]) &&
          ((int)HEIKIN_G[4] == (int)U_JYUSHIN[4]) &&
          ((int)HEIKIN_G[5] == (int)U_JYUSHIN[5]) &&
          ((int)HEIKIN_G[6] == (int)U_JYUSHIN[6]) &&
          ((int)HEIKIN_G[7] == (int)U_JYUSHIN[7]) &&
          ((int)HEIKIN_G[8] == (int)U_JYUSHIN[8]) &&
          ((int)HEIKIN_G[9] == (int)U_JYUSHIN[9]) &&
          ((int)HEIKIN_G[10] == (int)U_JYUSHIN[10]) &&
          ((int)HEIKIN_G[11] == (int)U_JYUSHIN[11]) &&
          ((int)HEIKIN_G[12] == (int)U_JYUSHIN[12]) &&
          ((int)HEIKIN_G[13] == (int)U_JYUSHIN[13]) &&
          ((int)HEIKIN_G[14] == (int)U_JYUSHIN[14]) &&
          ((int)HEIKIN_G[15] == (int)U_JYUSHIN[15]) &&
          ((int)HEIKIN_G[16] == (int)U_JYUSHIN[16]) &&
          ((int)HEIKIN_G[17] == (int)U_JYUSHIN[17]) &&
          ((int)HEIKIN_G[18] == (int)U_JYUSHIN[18]) &&
          ((int)HEIKIN_G[19] == (int)U_JYUSHIN[19]) &&
          ((int)HEIKIN_G[20] == (int)U_JYUSHIN[20]) &&
          ((int)HEIKIN_G[21] == (int)U_JYUSHIN[21]) &&
          ((int)HEIKIN_G[22] == (int)U_JYUSHIN[22]) &&
          ((int)HEIKIN_G[23] == (int)U_JYUSHIN[23]) &&
          ((int)HEIKIN_G[24] == (int)U_JYUSHIN[24]) &&
          ((int)HEIKIN_G[25] == (int)U_JYUSHIN[25]) &&
          ((int)HEIKIN_G[26] == (int)U_JYUSHIN[26]) &&
          ((int)HEIKIN_G[27] == (int)U_JYUSHIN[27]) &&
          ((int)HEIKIN_G[28] == (int)U_JYUSHIN[28]) &&
          ((int)HEIKIN_G[29] == (int)U_JYUSHIN[29]) &&
          ((int)HEIKIN_G[30] == (int)U_JYUSHIN[30]) &&
          ((int)HEIKIN_G[31] == (int)U_JYUSHIN[31]) &&
          ((int)HEIKIN_G[32] == (int)U_JYUSHIN[32]) &&
          ((int)HEIKIN_G[33] == (int)U_JYUSHIN[33]) &&
          ((int)HEIKIN_G[34] == (int)U_JYUSHIN[34]) &&
          ((int)HEIKIN_G[35] == (int)U_JYUSHIN[35]) &&
          ((int)HEIKIN_G[36] == (int)U_JYUSHIN[36]) &&
          ((int)HEIKIN_G[37] == (int)U_JYUSHIN[37]) &&
          ((int)HEIKIN_G[38] == (int)U_JYUSHIN[38]) &&
          ((int)HEIKIN_G[39] == (int)U_JYUSHIN[39]) &&
          ((int)HEIKIN_G[40] == (int)U_JYUSHIN[40]) &&
          ((int)HEIKIN_G[41] == (int)U_JYUSHIN[41]) &&
          ((int)HEIKIN_G[42] == (int)U_JYUSHIN[42]) &&
          ((int)HEIKIN_G[43] == (int)U_JYUSHIN[43]) &&
          ((int)HEIKIN_G[44] == (int)U_JYUSHIN[44]) &&
          ((int)HEIKIN_G[45] == (int)U_JYUSHIN[45]) &&
          ((int)HEIKIN_G[46] == (int)U_JYUSHIN[46]) &&
          ((int)HEIKIN_G[47] == (int)U_JYUSHIN[47]) &&
          ((int)HEIKIN_G[48] == (int)U_JYUSHIN[48]) &&
          ((int)HEIKIN_G[49] == (int)U_JYUSHIN[49]) &&
          ((int)HEIKIN_G[50] == (int)U_JYUSHIN[50]) &&
          ((int)HEIKIN_G[51] == (int)U_JYUSHIN[51]) &&
          ((int)HEIKIN_G[52] == (int)U_JYUSHIN[52]) &&
          ((int)HEIKIN_G[53] == (int)U_JYUSHIN[53]) &&
          ((int)HEIKIN_G[54] == (int)U_JYUSHIN[54]) &&
          ((int)HEIKIN_G[55] == (int)U_JYUSHIN[55]) &&
          ((int)HEIKIN_G[56] == (int)U_JYUSHIN[56]) &&
          ((int)HEIKIN_G[57] == (int)U_JYUSHIN[57]) &&
          ((int)HEIKIN_G[58] == (int)U_JYUSHIN[58]) &&
          ((int)HEIKIN_G[59] == (int)U_JYUSHIN[59]) &&
          ((int)HEIKIN_G[60] == (int)U_JYUSHIN[60]) &&
          ((int)HEIKIN_G[61] == (int)U_JYUSHIN[61]) &&
          ((int)HEIKIN_G[62] == (int)U_JYUSHIN[62]) &&
          ((int)HEIKIN_G[63] == (int)U_JYUSHIN[63]) &&
          ((int)HEIKIN_G[64] == (int)U_JYUSHIN[64]) &&
          ((int)HEIKIN_G[65] == (int)U_JYUSHIN[65]) &&
          ((int)HEIKIN_G[66] == (int)U_JYUSHIN[66]) &&
          ((int)HEIKIN_G[67] == (int)U_JYUSHIN[67]) &&
          ((int)HEIKIN_G[68] == (int)U_JYUSHIN[68]) &&
          ((int)HEIKIN_G[69] == (int)U_JYUSHIN[69]) &&
          ((int)HEIKIN_G[70] == (int)U_JYUSHIN[70]) &&
          ((int)HEIKIN_G[71] == (int)U_JYUSHIN[71]) &&
          ((int)HEIKIN_G[72] == (int)U_JYUSHIN[72]) &&
          ((int)HEIKIN_G[73] == (int)U_JYUSHIN[73]) &&
          ((int)HEIKIN_G[74] == (int)U_JYUSHIN[74]) &&
          ((int)HEIKIN_G[75] == (int)U_JYUSHIN[75]) &&
          ((int)HEIKIN_G[76] == (int)U_JYUSHIN[76]) &&
          ((int)HEIKIN_G[77] == (int)U_JYUSHIN[77]) &&
          ((int)HEIKIN_G[78] == (int)U_JYUSHIN[78]) &&
          ((int)HEIKIN_G[79] == (int)U_JYUSHIN[79]) &&
          ((int)HEIKIN_G[80] == (int)U_JYUSHIN[80]) &&
          ((int)HEIKIN_G[81] == (int)U_JYUSHIN[81]) &&
          ((int)HEIKIN_G[82] == (int)U_JYUSHIN[82]) &&
          ((int)HEIKIN_G[83] == (int)U_JYUSHIN[83]) &&
          ((int)HEIKIN_G[84] == (int)U_JYUSHIN[84]) &&
          ((int)HEIKIN_G[85] == (int)U_JYUSHIN[85]) &&
          ((int)HEIKIN_G[86] == (int)U_JYUSHIN[86]) &&
          ((int)HEIKIN_G[87] == (int)U_JYUSHIN[87]) &&
          ((int)HEIKIN_G[88] == (int)U_JYUSHIN[88]) &&
          ((int)HEIKIN_G[89] == (int)U_JYUSHIN[89]) &&
          ((int)HEIKIN_G[90] == (int)U_JYUSHIN[90]) &&
          ((int)HEIKIN_G[91] == (int)U_JYUSHIN[91]) &&
          ((int)HEIKIN_G[92] == (int)U_JYUSHIN[92]) &&
          ((int)HEIKIN_G[93] == (int)U_JYUSHIN[93]) &&
          ((int)HEIKIN_G[94] == (int)U_JYUSHIN[94]) &&
          ((int)HEIKIN_G[95] == (int)U_JYUSHIN[95]) &&
          ((int)HEIKIN_G[96] == (int)U_JYUSHIN[96]) &&
          ((int)HEIKIN_G[97] == (int)U_JYUSHIN[97]) &&
          ((int)HEIKIN_G[98] == (int)U_JYUSHIN[98]) &&
          ((int)HEIKIN_G[99] == (int)U_JYUSHIN[99]) &&
          ((int)HEIKIN_G[100] == (int)U_JYUSHIN[100]) &&
          ((int)HEIKIN_G[101] == (int)U_JYUSHIN[101]) &&
          ((int)HEIKIN_G[102] == (int)U_JYUSHIN[102]) &&
          ((int)HEIKIN_G[103] == (int)U_JYUSHIN[103]) &&
          ((int)HEIKIN_G[104] == (int)U_JYUSHIN[104]) &&
          ((int)HEIKIN_G[105] == (int)U_JYUSHIN[105]) &&
          ((int)HEIKIN_G[106] == (int)U_JYUSHIN[106]) &&
          ((int)HEIKIN_G[107] == (int)U_JYUSHIN[107]) &&
          ((int)HEIKIN_G[108] == (int)U_JYUSHIN[108]) &&
          ((int)HEIKIN_G[109] == (int)U_JYUSHIN[109]) &&
          ((int)HEIKIN_G[110] == (int)U_JYUSHIN[110]) &&
          ((int)HEIKIN_G[111] == (int)U_JYUSHIN[111]) &&
          ((int)HEIKIN_G[112] == (int)U_JYUSHIN[112]) &&
          ((int)HEIKIN_G[113] == (int)U_JYUSHIN[113]) &&
          ((int)HEIKIN_G[114] == (int)U_JYUSHIN[114]) &&
          ((int)HEIKIN_G[115] == (int)U_JYUSHIN[115]) &&
          ((int)HEIKIN_G[116] == (int)U_JYUSHIN[116]) &&
          ((int)HEIKIN_G[117] == (int)U_JYUSHIN[117]) &&
          ((int)HEIKIN_G[118] == (int)U_JYUSHIN[118]) &&
          ((int)HEIKIN_G[119] == (int)U_JYUSHIN[119]) &&
          ((int)HEIKIN_G[120] == (int)U_JYUSHIN[120]) &&
          ((int)HEIKIN_G[121] == (int)U_JYUSHIN[121]) &&
          ((int)HEIKIN_G[122] == (int)U_JYUSHIN[122]) &&
          ((int)HEIKIN_G[123] == (int)U_JYUSHIN[123]) &&
          ((int)HEIKIN_G[124] == (int)U_JYUSHIN[124]) &&
          ((int)HEIKIN_G[125] == (int)U_JYUSHIN[125]) &&
          ((int)HEIKIN_G[126] == (int)U_JYUSHIN[126]) &&
          ((int)HEIKIN_G[127] == (int)U_JYUSHIN[127]) &&
          ((int)HEIKIN_G[128] == (int)U_JYUSHIN[128]) &&
          ((int)HEIKIN_G[129] == (int)U_JYUSHIN[129]) &&
          ((int)HEIKIN_G[130] == (int)U_JYUSHIN[130]) &&
          ((int)HEIKIN_G[131] == (int)U_JYUSHIN[131]) &&
          ((int)HEIKIN_G[132] == (int)U_JYUSHIN[132]) &&
          ((int)HEIKIN_G[133] == (int)U_JYUSHIN[133]) &&
          ((int)HEIKIN_G[134] == (int)U_JYUSHIN[134]) &&
          ((int)HEIKIN_G[135] == (int)U_JYUSHIN[135]) &&
          ((int)HEIKIN_G[136] == (int)U_JYUSHIN[136]) &&
          ((int)HEIKIN_G[137] == (int)U_JYUSHIN[137]) &&
          ((int)HEIKIN_G[138] == (int)U_JYUSHIN[138]) &&
          ((int)HEIKIN_G[139] == (int)U_JYUSHIN[139]) &&
          ((int)HEIKIN_G[140] == (int)U_JYUSHIN[140]) &&
          ((int)HEIKIN_G[141] == (int)U_JYUSHIN[141]) &&
          ((int)HEIKIN_G[142] == (int)U_JYUSHIN[142]) &&
          ((int)HEIKIN_G[143] == (int)U_JYUSHIN[143]) &&
          ((int)HEIKIN_G[144] == (int)U_JYUSHIN[144]) &&
          ((int)HEIKIN_G[145] == (int)U_JYUSHIN[145]) &&
          ((int)HEIKIN_G[146] == (int)U_JYUSHIN[146]) &&
          ((int)HEIKIN_G[147] == (int)U_JYUSHIN[147]) &&
          ((int)HEIKIN_G[148] == (int)U_JYUSHIN[148]) &&
          ((int)HEIKIN_G[149] == (int)U_JYUSHIN[149]) &&
          ((int)HEIKIN_G[150] == (int)U_JYUSHIN[150]) &&
          ((int)HEIKIN_G[151] == (int)U_JYUSHIN[151]) &&
          ((int)HEIKIN_G[152] == (int)U_JYUSHIN[152]) &&
          ((int)HEIKIN_G[153] == (int)U_JYUSHIN[153]) &&
          ((int)HEIKIN_G[154] == (int)U_JYUSHIN[154]) &&
          ((int)HEIKIN_G[155] == (int)U_JYUSHIN[155]) &&
          ((int)HEIKIN_G[156] == (int)U_JYUSHIN[156]) &&
          ((int)HEIKIN_G[157] == (int)U_JYUSHIN[157]) &&
          ((int)HEIKIN_G[158] == (int)U_JYUSHIN[158]) &&
          ((int)HEIKIN_G[159] == (int)U_JYUSHIN[159]) &&
          ((int)HEIKIN_G[160] == (int)U_JYUSHIN[160]) &&
          ((int)HEIKIN_G[161] == (int)U_JYUSHIN[161]) &&
          ((int)HEIKIN_G[162] == (int)U_JYUSHIN[162]) &&
          ((int)HEIKIN_G[163] == (int)U_JYUSHIN[163]) &&
          ((int)HEIKIN_G[164] == (int)U_JYUSHIN[164]) &&
          ((int)HEIKIN_G[165] == (int)U_JYUSHIN[165]) &&
          ((int)HEIKIN_G[166] == (int)U_JYUSHIN[166]) &&
          ((int)HEIKIN_G[167] == (int)U_JYUSHIN[167]) &&
          ((int)HEIKIN_G[168] == (int)U_JYUSHIN[168]) &&
          ((int)HEIKIN_G[169] == (int)U_JYUSHIN[169]) &&
          ((int)HEIKIN_G[170] == (int)U_JYUSHIN[170]) &&
          ((int)HEIKIN_G[171] == (int)U_JYUSHIN[171]) &&
          ((int)HEIKIN_G[172] == (int)U_JYUSHIN[172]) &&
          ((int)HEIKIN_G[173] == (int)U_JYUSHIN[173]) &&
          ((int)HEIKIN_G[174] == (int)U_JYUSHIN[174]) &&
          ((int)HEIKIN_G[175] == (int)U_JYUSHIN[175]) &&
          ((int)HEIKIN_G[176] == (int)U_JYUSHIN[176]) &&
          ((int)HEIKIN_G[177] == (int)U_JYUSHIN[177]) &&
          ((int)HEIKIN_G[178] == (int)U_JYUSHIN[178]) &&
          ((int)HEIKIN_G[179] == (int)U_JYUSHIN[179]) &&
          ((int)HEIKIN_G[180] == (int)U_JYUSHIN[180]) &&
          ((int)HEIKIN_G[181] == (int)U_JYUSHIN[181]) &&
          ((int)HEIKIN_G[182] == (int)U_JYUSHIN[182]) &&
          ((int)HEIKIN_G[183] == (int)U_JYUSHIN[183]) &&
          ((int)HEIKIN_G[184] == (int)U_JYUSHIN[184]) &&
          ((int)HEIKIN_G[185] == (int)U_JYUSHIN[185]) &&
          ((int)HEIKIN_G[186] == (int)U_JYUSHIN[186]) &&
          ((int)HEIKIN_G[187] == (int)U_JYUSHIN[187]) &&
          ((int)HEIKIN_G[188] == (int)U_JYUSHIN[188]) &&
          ((int)HEIKIN_G[189] == (int)U_JYUSHIN[189]) &&
          ((int)HEIKIN_G[190] == (int)U_JYUSHIN[190]) &&
          ((int)HEIKIN_G[191] == (int)U_JYUSHIN[191]) &&
          ((int)HEIKIN_G[192] == (int)U_JYUSHIN[192]) &&
          ((int)HEIKIN_G[193] == (int)U_JYUSHIN[193]) &&
          ((int)HEIKIN_G[194] == (int)U_JYUSHIN[194]) &&
          ((int)HEIKIN_G[195] == (int)U_JYUSHIN[195]) &&
          ((int)HEIKIN_G[196] == (int)U_JYUSHIN[196]) &&
          ((int)HEIKIN_G[197] == (int)U_JYUSHIN[197]) &&
          ((int)HEIKIN_G[198] == (int)U_JYUSHIN[198]) &&
          ((int)HEIKIN_G[199] == (int)U_JYUSHIN[199]) &&
          ((int)HEIKIN_G[200] == (int)U_JYUSHIN[200]) &&
          ((int)HEIKIN_G[201] == (int)U_JYUSHIN[201]) &&
          ((int)HEIKIN_G[202] == (int)U_JYUSHIN[202]) &&
          ((int)HEIKIN_G[203] == (int)U_JYUSHIN[203]) &&
          ((int)HEIKIN_G[204] == (int)U_JYUSHIN[204]) &&
          ((int)HEIKIN_G[205] == (int)U_JYUSHIN[205]) &&
          ((int)HEIKIN_G[206] == (int)U_JYUSHIN[206]) &&
          ((int)HEIKIN_G[207] == (int)U_JYUSHIN[207]) &&
          ((int)HEIKIN_G[208] == (int)U_JYUSHIN[208]) &&
          ((int)HEIKIN_G[209] == (int)U_JYUSHIN[209]) &&
          ((int)HEIKIN_G[210] == (int)U_JYUSHIN[210]) &&
          ((int)HEIKIN_G[211] == (int)U_JYUSHIN[211]) &&
          ((int)HEIKIN_G[212] == (int)U_JYUSHIN[212]) &&
          ((int)HEIKIN_G[213] == (int)U_JYUSHIN[213]) &&
          ((int)HEIKIN_G[214] == (int)U_JYUSHIN[214]) &&
          ((int)HEIKIN_G[215] == (int)U_JYUSHIN[215]) &&
          ((int)HEIKIN_G[216] == (int)U_JYUSHIN[216]) &&
          ((int)HEIKIN_G[217] == (int)U_JYUSHIN[217]) &&
          ((int)HEIKIN_G[218] == (int)U_JYUSHIN[218]) &&
          ((int)HEIKIN_G[219] == (int)U_JYUSHIN[219]) &&
          ((int)HEIKIN_G[220] == (int)U_JYUSHIN[220]) &&
          ((int)HEIKIN_G[221] == (int)U_JYUSHIN[221]) &&
          ((int)HEIKIN_G[222] == (int)U_JYUSHIN[222]) &&
          ((int)HEIKIN_G[223] == (int)U_JYUSHIN[223]) &&
          ((int)HEIKIN_G[224] == (int)U_JYUSHIN[224]) &&
          ((int)HEIKIN_G[225] == (int)U_JYUSHIN[225]) &&
          ((int)HEIKIN_G[226] == (int)U_JYUSHIN[226]) &&
          ((int)HEIKIN_G[227] == (int)U_JYUSHIN[227]) &&
          ((int)HEIKIN_G[228] == (int)U_JYUSHIN[228]) &&
          ((int)HEIKIN_G[229] == (int)U_JYUSHIN[229]) &&
          ((int)HEIKIN_G[230] == (int)U_JYUSHIN[230]) &&
          ((int)HEIKIN_G[231] == (int)U_JYUSHIN[231]) &&
          ((int)HEIKIN_G[232] == (int)U_JYUSHIN[232]) &&
          ((int)HEIKIN_G[233] == (int)U_JYUSHIN[233]) &&
          ((int)HEIKIN_G[234] == (int)U_JYUSHIN[234]) &&
          ((int)HEIKIN_G[235] == (int)U_JYUSHIN[235]) &&
          ((int)HEIKIN_G[236] == (int)U_JYUSHIN[236]) &&
          ((int)HEIKIN_G[237] == (int)U_JYUSHIN[237]) &&
          ((int)HEIKIN_G[238] == (int)U_JYUSHIN[238]) &&
          ((int)HEIKIN_G[239] == (int)U_JYUSHIN[239]) &&
          ((int)HEIKIN_G[240] == (int)U_JYUSHIN[240]) &&
          ((int)HEIKIN_G[241] == (int)U_JYUSHIN[241]) &&
          ((int)HEIKIN_G[242] == (int)U_JYUSHIN[242]) &&
          ((int)HEIKIN_G[243] == (int)U_JYUSHIN[243]) &&
          ((int)HEIKIN_G[244] == (int)U_JYUSHIN[244]) &&
          ((int)HEIKIN_G[245] == (int)U_JYUSHIN[245]) &&
          ((int)HEIKIN_G[246] == (int)U_JYUSHIN[246]) &&
          ((int)HEIKIN_G[247] == (int)U_JYUSHIN[247]) &&
          ((int)HEIKIN_G[248] == (int)U_JYUSHIN[248]) &&
          ((int)HEIKIN_G[249] == (int)U_JYUSHIN[249]) &&
          ((int)HEIKIN_G[250] == (int)U_JYUSHIN[250]) &&
          ((int)HEIKIN_G[251] == (int)U_JYUSHIN[251]) &&
          ((int)HEIKIN_G[252] == (int)U_JYUSHIN[252]) &&
          ((int)HEIKIN_G[253] == (int)U_JYUSHIN[253]) &&
          ((int)HEIKIN_G[254] == (int)U_JYUSHIN[254]) &&
          ((int)HEIKIN_G[255] == (int)U_JYUSHIN[255]) &&
          ((int)HEIKIN_B[0] == (int)V_JYUSHIN[0]) &&
          ((int)HEIKIN_B[1] == (int)V_JYUSHIN[1]) &&
          ((int)HEIKIN_B[2] == (int)V_JYUSHIN[2]) &&
          ((int)HEIKIN_B[3] == (int)V_JYUSHIN[3]) &&
          ((int)HEIKIN_B[4] == (int)V_JYUSHIN[4]) &&
          ((int)HEIKIN_B[5] == (int)V_JYUSHIN[5]) &&
          ((int)HEIKIN_B[6] == (int)V_JYUSHIN[6]) &&
          ((int)HEIKIN_B[7] == (int)V_JYUSHIN[7]) &&
          ((int)HEIKIN_B[8] == (int)V_JYUSHIN[8]) &&
          ((int)HEIKIN_B[9] == (int)V_JYUSHIN[9]) &&
          ((int)HEIKIN_B[10] == (int)V_JYUSHIN[10]) &&
          ((int)HEIKIN_B[11] == (int)V_JYUSHIN[11]) &&
          ((int)HEIKIN_B[12] == (int)V_JYUSHIN[12]) &&
          ((int)HEIKIN_B[13] == (int)V_JYUSHIN[13]) &&
          ((int)HEIKIN_B[14] == (int)V_JYUSHIN[14]) &&
          ((int)HEIKIN_B[15] == (int)V_JYUSHIN[15]) &&
          ((int)HEIKIN_B[16] == (int)V_JYUSHIN[16]) &&
          ((int)HEIKIN_B[17] == (int)V_JYUSHIN[17]) &&
          ((int)HEIKIN_B[18] == (int)V_JYUSHIN[18]) &&
          ((int)HEIKIN_B[19] == (int)V_JYUSHIN[19]) &&
          ((int)HEIKIN_B[20] == (int)V_JYUSHIN[20]) &&
          ((int)HEIKIN_B[21] == (int)V_JYUSHIN[21]) &&
          ((int)HEIKIN_B[22] == (int)V_JYUSHIN[22]) &&
          ((int)HEIKIN_B[23] == (int)V_JYUSHIN[23]) &&
          ((int)HEIKIN_B[24] == (int)V_JYUSHIN[24]) &&
          ((int)HEIKIN_B[25] == (int)V_JYUSHIN[25]) &&
          ((int)HEIKIN_B[26] == (int)V_JYUSHIN[26]) &&
          ((int)HEIKIN_B[27] == (int)V_JYUSHIN[27]) &&
          ((int)HEIKIN_B[28] == (int)V_JYUSHIN[28]) &&
          ((int)HEIKIN_B[29] == (int)V_JYUSHIN[29]) &&
          ((int)HEIKIN_B[30] == (int)V_JYUSHIN[30]) &&
          ((int)HEIKIN_B[31] == (int)V_JYUSHIN[31]) &&
          ((int)HEIKIN_B[32] == (int)V_JYUSHIN[32]) &&
          ((int)HEIKIN_B[33] == (int)V_JYUSHIN[33]) &&
          ((int)HEIKIN_B[34] == (int)V_JYUSHIN[34]) &&
          ((int)HEIKIN_B[35] == (int)V_JYUSHIN[35]) &&
          ((int)HEIKIN_B[36] == (int)V_JYUSHIN[36]) &&
          ((int)HEIKIN_B[37] == (int)V_JYUSHIN[37]) &&
          ((int)HEIKIN_B[38] == (int)V_JYUSHIN[38]) &&
          ((int)HEIKIN_B[39] == (int)V_JYUSHIN[39]) &&
          ((int)HEIKIN_B[40] == (int)V_JYUSHIN[40]) &&
          ((int)HEIKIN_B[41] == (int)V_JYUSHIN[41]) &&
          ((int)HEIKIN_B[42] == (int)V_JYUSHIN[42]) &&
          ((int)HEIKIN_B[43] == (int)V_JYUSHIN[43]) &&
          ((int)HEIKIN_B[44] == (int)V_JYUSHIN[44]) &&
          ((int)HEIKIN_B[45] == (int)V_JYUSHIN[45]) &&
          ((int)HEIKIN_B[46] == (int)V_JYUSHIN[46]) &&
          ((int)HEIKIN_B[47] == (int)V_JYUSHIN[47]) &&
          ((int)HEIKIN_B[48] == (int)V_JYUSHIN[48]) &&
          ((int)HEIKIN_B[49] == (int)V_JYUSHIN[49]) &&
          ((int)HEIKIN_B[50] == (int)V_JYUSHIN[50]) &&
          ((int)HEIKIN_B[51] == (int)V_JYUSHIN[51]) &&
          ((int)HEIKIN_B[52] == (int)V_JYUSHIN[52]) &&
          ((int)HEIKIN_B[53] == (int)V_JYUSHIN[53]) &&
          ((int)HEIKIN_B[54] == (int)V_JYUSHIN[54]) &&
          ((int)HEIKIN_B[55] == (int)V_JYUSHIN[55]) &&
          ((int)HEIKIN_B[56] == (int)V_JYUSHIN[56]) &&
          ((int)HEIKIN_B[57] == (int)V_JYUSHIN[57]) &&
          ((int)HEIKIN_B[58] == (int)V_JYUSHIN[58]) &&
          ((int)HEIKIN_B[59] == (int)V_JYUSHIN[59]) &&
          ((int)HEIKIN_B[60] == (int)V_JYUSHIN[60]) &&
          ((int)HEIKIN_B[61] == (int)V_JYUSHIN[61]) &&
          ((int)HEIKIN_B[62] == (int)V_JYUSHIN[62]) &&
          ((int)HEIKIN_B[63] == (int)V_JYUSHIN[63]) &&
          ((int)HEIKIN_B[64] == (int)V_JYUSHIN[64]) &&
          ((int)HEIKIN_B[65] == (int)V_JYUSHIN[65]) &&
          ((int)HEIKIN_B[66] == (int)V_JYUSHIN[66]) &&
          ((int)HEIKIN_B[67] == (int)V_JYUSHIN[67]) &&
          ((int)HEIKIN_B[68] == (int)V_JYUSHIN[68]) &&
          ((int)HEIKIN_B[69] == (int)V_JYUSHIN[69]) &&
          ((int)HEIKIN_B[70] == (int)V_JYUSHIN[70]) &&
          ((int)HEIKIN_B[71] == (int)V_JYUSHIN[71]) &&
          ((int)HEIKIN_B[72] == (int)V_JYUSHIN[72]) &&
          ((int)HEIKIN_B[73] == (int)V_JYUSHIN[73]) &&
          ((int)HEIKIN_B[74] == (int)V_JYUSHIN[74]) &&
          ((int)HEIKIN_B[75] == (int)V_JYUSHIN[75]) &&
          ((int)HEIKIN_B[76] == (int)V_JYUSHIN[76]) &&
          ((int)HEIKIN_B[77] == (int)V_JYUSHIN[77]) &&
          ((int)HEIKIN_B[78] == (int)V_JYUSHIN[78]) &&
          ((int)HEIKIN_B[79] == (int)V_JYUSHIN[79]) &&
          ((int)HEIKIN_B[80] == (int)V_JYUSHIN[80]) &&
          ((int)HEIKIN_B[81] == (int)V_JYUSHIN[81]) &&
          ((int)HEIKIN_B[82] == (int)V_JYUSHIN[82]) &&
          ((int)HEIKIN_B[83] == (int)V_JYUSHIN[83]) &&
          ((int)HEIKIN_B[84] == (int)V_JYUSHIN[84]) &&
          ((int)HEIKIN_B[85] == (int)V_JYUSHIN[85]) &&
          ((int)HEIKIN_B[86] == (int)V_JYUSHIN[86]) &&
          ((int)HEIKIN_B[87] == (int)V_JYUSHIN[87]) &&
          ((int)HEIKIN_B[88] == (int)V_JYUSHIN[88]) &&
          ((int)HEIKIN_B[89] == (int)V_JYUSHIN[89]) &&
          ((int)HEIKIN_B[90] == (int)V_JYUSHIN[90]) &&
          ((int)HEIKIN_B[91] == (int)V_JYUSHIN[91]) &&
          ((int)HEIKIN_B[92] == (int)V_JYUSHIN[92]) &&
          ((int)HEIKIN_B[93] == (int)V_JYUSHIN[93]) &&
          ((int)HEIKIN_B[94] == (int)V_JYUSHIN[94]) &&
          ((int)HEIKIN_B[95] == (int)V_JYUSHIN[95]) &&
          ((int)HEIKIN_B[96] == (int)V_JYUSHIN[96]) &&
          ((int)HEIKIN_B[97] == (int)V_JYUSHIN[97]) &&
          ((int)HEIKIN_B[98] == (int)V_JYUSHIN[98]) &&
          ((int)HEIKIN_B[99] == (int)V_JYUSHIN[99]) &&
          ((int)HEIKIN_B[100] == (int)V_JYUSHIN[100]) &&
          ((int)HEIKIN_B[101] == (int)V_JYUSHIN[101]) &&
          ((int)HEIKIN_B[102] == (int)V_JYUSHIN[102]) &&
          ((int)HEIKIN_B[103] == (int)V_JYUSHIN[103]) &&
          ((int)HEIKIN_B[104] == (int)V_JYUSHIN[104]) &&
          ((int)HEIKIN_B[105] == (int)V_JYUSHIN[105]) &&
          ((int)HEIKIN_B[106] == (int)V_JYUSHIN[106]) &&
          ((int)HEIKIN_B[107] == (int)V_JYUSHIN[107]) &&
          ((int)HEIKIN_B[108] == (int)V_JYUSHIN[108]) &&
          ((int)HEIKIN_B[109] == (int)V_JYUSHIN[109]) &&
          ((int)HEIKIN_B[110] == (int)V_JYUSHIN[110]) &&
          ((int)HEIKIN_B[111] == (int)V_JYUSHIN[111]) &&
          ((int)HEIKIN_B[112] == (int)V_JYUSHIN[112]) &&
          ((int)HEIKIN_B[113] == (int)V_JYUSHIN[113]) &&
          ((int)HEIKIN_B[114] == (int)V_JYUSHIN[114]) &&
          ((int)HEIKIN_B[115] == (int)V_JYUSHIN[115]) &&
          ((int)HEIKIN_B[116] == (int)V_JYUSHIN[116]) &&
          ((int)HEIKIN_B[117] == (int)V_JYUSHIN[117]) &&
          ((int)HEIKIN_B[118] == (int)V_JYUSHIN[118]) &&
          ((int)HEIKIN_B[119] == (int)V_JYUSHIN[119]) &&
          ((int)HEIKIN_B[120] == (int)V_JYUSHIN[120]) &&
          ((int)HEIKIN_B[121] == (int)V_JYUSHIN[121]) &&
          ((int)HEIKIN_B[122] == (int)V_JYUSHIN[122]) &&
          ((int)HEIKIN_B[123] == (int)V_JYUSHIN[123]) &&
          ((int)HEIKIN_B[124] == (int)V_JYUSHIN[124]) &&
          ((int)HEIKIN_B[125] == (int)V_JYUSHIN[125]) &&
          ((int)HEIKIN_B[126] == (int)V_JYUSHIN[126]) &&
          ((int)HEIKIN_B[127] == (int)V_JYUSHIN[127]) &&
          ((int)HEIKIN_B[128] == (int)V_JYUSHIN[128]) &&
          ((int)HEIKIN_B[129] == (int)V_JYUSHIN[129]) &&
          ((int)HEIKIN_B[130] == (int)V_JYUSHIN[130]) &&
          ((int)HEIKIN_B[131] == (int)V_JYUSHIN[131]) &&
          ((int)HEIKIN_B[132] == (int)V_JYUSHIN[132]) &&
          ((int)HEIKIN_B[133] == (int)V_JYUSHIN[133]) &&
          ((int)HEIKIN_B[134] == (int)V_JYUSHIN[134]) &&
          ((int)HEIKIN_B[135] == (int)V_JYUSHIN[135]) &&
          ((int)HEIKIN_B[136] == (int)V_JYUSHIN[136]) &&
          ((int)HEIKIN_B[137] == (int)V_JYUSHIN[137]) &&
          ((int)HEIKIN_B[138] == (int)V_JYUSHIN[138]) &&
          ((int)HEIKIN_B[139] == (int)V_JYUSHIN[139]) &&
          ((int)HEIKIN_B[140] == (int)V_JYUSHIN[140]) &&
          ((int)HEIKIN_B[141] == (int)V_JYUSHIN[141]) &&
          ((int)HEIKIN_B[142] == (int)V_JYUSHIN[142]) &&
          ((int)HEIKIN_B[143] == (int)V_JYUSHIN[143]) &&
          ((int)HEIKIN_B[144] == (int)V_JYUSHIN[144]) &&
          ((int)HEIKIN_B[145] == (int)V_JYUSHIN[145]) &&
          ((int)HEIKIN_B[146] == (int)V_JYUSHIN[146]) &&
          ((int)HEIKIN_B[147] == (int)V_JYUSHIN[147]) &&
          ((int)HEIKIN_B[148] == (int)V_JYUSHIN[148]) &&
          ((int)HEIKIN_B[149] == (int)V_JYUSHIN[149]) &&
          ((int)HEIKIN_B[150] == (int)V_JYUSHIN[150]) &&
          ((int)HEIKIN_B[151] == (int)V_JYUSHIN[151]) &&
          ((int)HEIKIN_B[152] == (int)V_JYUSHIN[152]) &&
          ((int)HEIKIN_B[153] == (int)V_JYUSHIN[153]) &&
          ((int)HEIKIN_B[154] == (int)V_JYUSHIN[154]) &&
          ((int)HEIKIN_B[155] == (int)V_JYUSHIN[155]) &&
          ((int)HEIKIN_B[156] == (int)V_JYUSHIN[156]) &&
          ((int)HEIKIN_B[157] == (int)V_JYUSHIN[157]) &&
          ((int)HEIKIN_B[158] == (int)V_JYUSHIN[158]) &&
          ((int)HEIKIN_B[159] == (int)V_JYUSHIN[159]) &&
          ((int)HEIKIN_B[160] == (int)V_JYUSHIN[160]) &&
          ((int)HEIKIN_B[161] == (int)V_JYUSHIN[161]) &&
          ((int)HEIKIN_B[162] == (int)V_JYUSHIN[162]) &&
          ((int)HEIKIN_B[163] == (int)V_JYUSHIN[163]) &&
          ((int)HEIKIN_B[164] == (int)V_JYUSHIN[164]) &&
          ((int)HEIKIN_B[165] == (int)V_JYUSHIN[165]) &&
          ((int)HEIKIN_B[166] == (int)V_JYUSHIN[166]) &&
          ((int)HEIKIN_B[167] == (int)V_JYUSHIN[167]) &&
          ((int)HEIKIN_B[168] == (int)V_JYUSHIN[168]) &&
          ((int)HEIKIN_B[169] == (int)V_JYUSHIN[169]) &&
          ((int)HEIKIN_B[170] == (int)V_JYUSHIN[170]) &&
          ((int)HEIKIN_B[171] == (int)V_JYUSHIN[171]) &&
          ((int)HEIKIN_B[172] == (int)V_JYUSHIN[172]) &&
          ((int)HEIKIN_B[173] == (int)V_JYUSHIN[173]) &&
          ((int)HEIKIN_B[174] == (int)V_JYUSHIN[174]) &&
          ((int)HEIKIN_B[175] == (int)V_JYUSHIN[175]) &&
          ((int)HEIKIN_B[176] == (int)V_JYUSHIN[176]) &&
          ((int)HEIKIN_B[177] == (int)V_JYUSHIN[177]) &&
          ((int)HEIKIN_B[178] == (int)V_JYUSHIN[178]) &&
          ((int)HEIKIN_B[179] == (int)V_JYUSHIN[179]) &&
          ((int)HEIKIN_B[180] == (int)V_JYUSHIN[180]) &&
          ((int)HEIKIN_B[181] == (int)V_JYUSHIN[181]) &&
          ((int)HEIKIN_B[182] == (int)V_JYUSHIN[182]) &&
          ((int)HEIKIN_B[183] == (int)V_JYUSHIN[183]) &&
          ((int)HEIKIN_B[184] == (int)V_JYUSHIN[184]) &&
          ((int)HEIKIN_B[185] == (int)V_JYUSHIN[185]) &&
          ((int)HEIKIN_B[186] == (int)V_JYUSHIN[186]) &&
          ((int)HEIKIN_B[187] == (int)V_JYUSHIN[187]) &&
          ((int)HEIKIN_B[188] == (int)V_JYUSHIN[188]) &&
          ((int)HEIKIN_B[189] == (int)V_JYUSHIN[189]) &&
          ((int)HEIKIN_B[190] == (int)V_JYUSHIN[190]) &&
          ((int)HEIKIN_B[191] == (int)V_JYUSHIN[191]) &&
          ((int)HEIKIN_B[192] == (int)V_JYUSHIN[192]) &&
          ((int)HEIKIN_B[193] == (int)V_JYUSHIN[193]) &&
          ((int)HEIKIN_B[194] == (int)V_JYUSHIN[194]) &&
          ((int)HEIKIN_B[195] == (int)V_JYUSHIN[195]) &&
          ((int)HEIKIN_B[196] == (int)V_JYUSHIN[196]) &&
          ((int)HEIKIN_B[197] == (int)V_JYUSHIN[197]) &&
          ((int)HEIKIN_B[198] == (int)V_JYUSHIN[198]) &&
          ((int)HEIKIN_B[199] == (int)V_JYUSHIN[199]) &&
          ((int)HEIKIN_B[200] == (int)V_JYUSHIN[200]) &&
          ((int)HEIKIN_B[201] == (int)V_JYUSHIN[201]) &&
          ((int)HEIKIN_B[202] == (int)V_JYUSHIN[202]) &&
          ((int)HEIKIN_B[203] == (int)V_JYUSHIN[203]) &&
          ((int)HEIKIN_B[204] == (int)V_JYUSHIN[204]) &&
          ((int)HEIKIN_B[205] == (int)V_JYUSHIN[205]) &&
          ((int)HEIKIN_B[206] == (int)V_JYUSHIN[206]) &&
          ((int)HEIKIN_B[207] == (int)V_JYUSHIN[207]) &&
          ((int)HEIKIN_B[208] == (int)V_JYUSHIN[208]) &&
          ((int)HEIKIN_B[209] == (int)V_JYUSHIN[209]) &&
          ((int)HEIKIN_B[210] == (int)V_JYUSHIN[210]) &&
          ((int)HEIKIN_B[211] == (int)V_JYUSHIN[211]) &&
          ((int)HEIKIN_B[212] == (int)V_JYUSHIN[212]) &&
          ((int)HEIKIN_B[213] == (int)V_JYUSHIN[213]) &&
          ((int)HEIKIN_B[214] == (int)V_JYUSHIN[214]) &&
          ((int)HEIKIN_B[215] == (int)V_JYUSHIN[215]) &&
          ((int)HEIKIN_B[216] == (int)V_JYUSHIN[216]) &&
          ((int)HEIKIN_B[217] == (int)V_JYUSHIN[217]) &&
          ((int)HEIKIN_B[218] == (int)V_JYUSHIN[218]) &&
          ((int)HEIKIN_B[219] == (int)V_JYUSHIN[219]) &&
          ((int)HEIKIN_B[220] == (int)V_JYUSHIN[220]) &&
          ((int)HEIKIN_B[221] == (int)V_JYUSHIN[221]) &&
          ((int)HEIKIN_B[222] == (int)V_JYUSHIN[222]) &&
          ((int)HEIKIN_B[223] == (int)V_JYUSHIN[223]) &&
          ((int)HEIKIN_B[224] == (int)V_JYUSHIN[224]) &&
          ((int)HEIKIN_B[225] == (int)V_JYUSHIN[225]) &&
          ((int)HEIKIN_B[226] == (int)V_JYUSHIN[226]) &&
          ((int)HEIKIN_B[227] == (int)V_JYUSHIN[227]) &&
          ((int)HEIKIN_B[228] == (int)V_JYUSHIN[228]) &&
          ((int)HEIKIN_B[229] == (int)V_JYUSHIN[229]) &&
          ((int)HEIKIN_B[230] == (int)V_JYUSHIN[230]) &&
          ((int)HEIKIN_B[231] == (int)V_JYUSHIN[231]) &&
          ((int)HEIKIN_B[232] == (int)V_JYUSHIN[232]) &&
          ((int)HEIKIN_B[233] == (int)V_JYUSHIN[233]) &&
          ((int)HEIKIN_B[234] == (int)V_JYUSHIN[234]) &&
          ((int)HEIKIN_B[235] == (int)V_JYUSHIN[235]) &&
          ((int)HEIKIN_B[236] == (int)V_JYUSHIN[236]) &&
          ((int)HEIKIN_B[237] == (int)V_JYUSHIN[237]) &&
          ((int)HEIKIN_B[238] == (int)V_JYUSHIN[238]) &&
          ((int)HEIKIN_B[239] == (int)V_JYUSHIN[239]) &&
          ((int)HEIKIN_B[240] == (int)V_JYUSHIN[240]) &&
          ((int)HEIKIN_B[241] == (int)V_JYUSHIN[241]) &&
          ((int)HEIKIN_B[242] == (int)V_JYUSHIN[242]) &&
          ((int)HEIKIN_B[243] == (int)V_JYUSHIN[243]) &&
          ((int)HEIKIN_B[244] == (int)V_JYUSHIN[244]) &&
          ((int)HEIKIN_B[245] == (int)V_JYUSHIN[245]) &&
          ((int)HEIKIN_B[246] == (int)V_JYUSHIN[246]) &&
          ((int)HEIKIN_B[247] == (int)V_JYUSHIN[247]) &&
          ((int)HEIKIN_B[248] == (int)V_JYUSHIN[248]) &&
          ((int)HEIKIN_B[249] == (int)V_JYUSHIN[249]) &&
          ((int)HEIKIN_B[250] == (int)V_JYUSHIN[250]) &&
          ((int)HEIKIN_B[251] == (int)V_JYUSHIN[251]) &&
          ((int)HEIKIN_B[252] == (int)V_JYUSHIN[252]) &&
          ((int)HEIKIN_B[253] == (int)V_JYUSHIN[253]) &&
          ((int)HEIKIN_B[254] == (int)V_JYUSHIN[254]) &&
          ((int)HEIKIN_B[255] == (int)V_JYUSHIN[255])) {
        fprintf(stderr, "%d iter\n", p);
        break;
      }
      for (i = 0; i < IROSUU; i++) {
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
  for (i = 0; i < hsize * vsize; i++) {
    // if(index3[i]!=-1){
    MINZ = 60e60;
    for (j = 0; j < IROSUU; j++) {
      Z[j] =
          ((((*(YIN + i)))) - Y_JYUSHIN[j]) *
              ((((*(YIN + i)))) - Y_JYUSHIN[j]) +
          ((((*(UIN + i)))) - U_JYUSHIN[j]) *
              ((((*(UIN + i)))) - U_JYUSHIN[j]) +
          ((((*(VIN + i)))) - V_JYUSHIN[j]) * ((((*(VIN + i)))) - V_JYUSHIN[j]);
    }
    for (j = 0; j < IROSUU; j++) {
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
  for (i = 0; i < IROSUU; i++) {
    Y_JYUSHIN1[i] = Y_JYUSHIN[i];
    U_JYUSHIN1[i] = U_JYUSHIN[i];
    V_JYUSHIN1[i] = V_JYUSHIN[i];
  }

  // YUVパレットが完成
  //これをＲＧＢに変換
  int IREDUCE_R[IROSUU];
  int IREDUCE_G[IROSUU];
  int IREDUCE_B[IROSUU];
  for (i = 0; i < IROSUU; i++) {
    // U_JYUSHIN[i] = 255.0*bmult*pow(U_JYUSHIN[i]/(bmult*255.0),1.0/bgamma);
    // V_JYUSHIN[i] = 255.0*rmult*pow(V_JYUSHIN[i]/(rmult*255.0),1.0/rgamma);
  }
  if (cs == 4 || cs == 9) {
    for (i = 0; i < IROSUU; i++) {
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
    for (i = 0; i < IROSUU; i++) {
      KARIU[i] = U_JYUSHIN[i] / ((float)(bmult)) * cos(IRO_THIETA) +
                 V_JYUSHIN[i] / ((float)(rmult)) * sin(IRO_THIETA);
      KARIV[i] = -U_JYUSHIN[i] / ((float)(bmult)) * sin(IRO_THIETA) +
                 V_JYUSHIN[i] / ((float)(rmult)) * cos(IRO_THIETA);
      U_JYUSHIN[i] = KARIU[i];
      V_JYUSHIN[i] = KARIV[i];
      Y_JYUSHIN[i] /= gmult;
    }
  }
  for (i = 0; i < IROSUU; i++) {
    Y_JYUSHIN[i] = 255.0 * pow(Y_JYUSHIN[i] / 255.0, 1.0 / ygamma);
  }
  double KAARII;
  double *FREDUCE_R, *FREDUCE_G, *FREDUCE_B;
  FREDUCE_R = (double *)malloc(sizeof(double) * IROSUU);
  FREDUCE_G = (double *)malloc(sizeof(double) * IROSUU);
  FREDUCE_B = (double *)malloc(sizeof(double) * IROSUU);

  if (cs == 0) {
    for (i = 0; i < IROSUU; i++) {
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
    for (i = 0; i < IROSUU; i++) {
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
    for (i = 0; i < IROSUU; i++) {
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
    for (i = 0; i < IROSUU; i++) {
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
    ///n1*0.9998247134019967+n2*-0.000044452970991771824+n3* 1.7528659800326482;
    ////            GG= n1*1.0000893144198046+n2*-0.4320813565838843+
    ///n3*-0.89314419804726; /            BB=
    ///n1*1.0000000115762946+n2* 2.213731098385467+ n3*-0.00011576294529099052;
    //}
    //}
  } else if (cs == 4) {
    for (i = 0; i < IROSUU; i++) {
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
    for (i = 0; i < IROSUU; i++) {
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
    for (i = 0; i < IROSUU; i++) {
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
    for (i = 0; i < IROSUU; i++) {
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
    for (i = 0; i < IROSUU; i++) {
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

  for (i = 0; i < IROSUU; i++) {
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
    for (i = 0; i < hsize * vsize; i++) {
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
      for (i = 0; i < mx * 2; i++) {
        *(errorR + i) = 0.0;
        *(errorG + i) = 0.0;
        *(errorB + i) = 0.0;
      }
    } else if (dither == 2 || dither == 3) {
      for (i = 0; i < mx * 3; i++) {
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
        for (i = 0; i < IROSUU; i++) {
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
        for (j = 0; j < mx; j++) {
          *(errorR + j) = *(errorR + j + mx);
          *(errorG + j) = *(errorG + j + mx);
          *(errorB + j) = *(errorB + j + mx);
          *(errorR + j + mx) = 0.0;
          *(errorG + j + mx) = 0.0;
          *(errorB + j + mx) = 0.0;
        }
      } else if (dither == 2 || dither == 3) {
        for (j = 0; j < mx; j++) {
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

    for (i = 0; i < hsize * vsize; i++) {
      for (j = 0; j < IROSUU; j++) {
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
    for (i = 0; i < hsize * vsize; i++) {
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
      for (i = 0; i < mx * 2; i++) {
        *(errorR + i) = 0.0;
        *(errorG + i) = 0.0;
        *(errorB + i) = 0.0;
      }
    } else if (dither == 2 || dither == 3) {
      for (i = 0; i < mx * 3; i++) {
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
        for (i = 0; i < IROSUU; i++) {
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
        for (j = 0; j < mx; j++) {
          *(errorR + j) = *(errorR + j + mx);
          *(errorG + j) = *(errorG + j + mx);
          *(errorB + j) = *(errorB + j + mx);
          *(errorR + j + mx) = 0.0;
          *(errorG + j + mx) = 0.0;
          *(errorB + j + mx) = 0.0;
        }
      } else if (dither == 2 || dither == 3) {
        for (j = 0; j < mx; j++) {
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

    for (i = 0; i < hsize * vsize; i++) {
      for (j = 0; j < IROSUU; j++) {
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

int main(int argc, char *argv[]) {
  CNTL *cnt;
  char inputfile[128];
  char outputfile[128];
  FILE *fpr;
  FILE *fpw;
  int hsize;
  int vsize;
  int filesizebyte;
  unsigned char *vin;
  unsigned char *B;
  unsigned char *G;
  unsigned char *R;
  int i, j, k, l, m, n, o, p;
  unsigned char *BIN;
  unsigned char *GIN;
  unsigned char *RIN;
  // unsigned char *BOUT;
  // unsigned char *GOUT;
  // unsigned char *ROUT;
  unsigned char *PALETGAZOU;
  unsigned char REDUCE_R[256];
  unsigned char REDUCE_G[256];
  unsigned char REDUCE_B[256];

  int filesize;
  unsigned char temp;
  cnt = (CNTL *)malloc(sizeof(CNTL) * sizeof(unsigned char));
  cnt->edge = 6;
  cnt->et = 20;
  cnt->bai = 1.0;
  cnt->kmon = 0;
  cnt->edon = 2;
  cnt->r = 4.0;
  cnt->g = 4.0;
  cnt->b = 3.0;
  cnt->cs = 3;
  cnt->gamma = 1.1;
  cnt->rgamma = 1.1;
  cnt->bgamma = 1.1;
  cnt->ygamma = 1.1;
  cnt->omh = 4;
  cnt->hasyori = 0.0;
  cnt->dither = 0;
  cnt->per = 0.75;
  cnt->pw = 2.0;
  cnt->pw2 = 1.0;
  cnt->bun = 1;
  cnt->lpfon = 0;
  cnt->t = 5.0;
  cnt->a = 0.001;
  cnt->lpfk = 0.8;
  if (argc < 2) {
    ;
    print_help(cnt);
    exit(0);
  }
  /*******************************************/
  /* */
  /* MAIN KANSUU NO HIKISUU NO SYORI */
  /* */
  /*******************************************/

  option_set(argc, argv, cnt);

  strcpy(inputfile, cnt->InputFileName);
  strcpy(outputfile, cnt->OutputFileName);

  /*******************************************/
  /* */
  /* BITMAP FILE OPEN */
  /* */
  /*******************************************/

  if ((fpr = fopen(inputfile, "rb")) == NULL) {
    exit(0);
  }
  fseek(fpr, 18, SEEK_SET);
  fread(&hsize, sizeof(int), 1, fpr); // hsize syutoku
  fread(&vsize, sizeof(int), 1, fpr); // vsize syutoku
  fseek(fpr, 2, SEEK_SET);
  fread(&filesizebyte, sizeof(int), 1, fpr);

  /*******************************************/
  /* */
  /* PICTURE GET */
  /* */
  /*******************************************/

  vin = (unsigned char *)malloc(sizeof(unsigned char) * (filesizebyte - 54));
  B = (unsigned char *)malloc(sizeof(unsigned char) * hsize * vsize);
  G = (unsigned char *)malloc(sizeof(unsigned char) * hsize * vsize);
  R = (unsigned char *)malloc(sizeof(unsigned char) * hsize * vsize);
  fseek(fpr, 54, SEEK_SET);
  if (fread(vin, sizeof(unsigned char), filesizebyte - 54, fpr) !=
      filesizebyte - 54) {
    fprintf(stderr, "InputFileSizeError!!!\n");
    exit(0);
  }
  fclose(fpr);
  /*******************************************/
  /* */
  /* SET B,G,R ARRAYS WITH RGB PIXEL DATAS */
  /* */
  /*******************************************/

  for (i = 0, k = 0; i < vsize; i++) {
    for (j = 0; j < hsize; j++) {
      *(B + j * vsize + i) = *(vin + k);
      k++;
      *(G + j * vsize + i) = *(vin + k);
      k++;
      *(R + j * vsize + i) = *(vin + k);
      k++;
    }
    if (((hsize * 3) % 4) == 1) {
      k = k + 3;
    }
    if (((hsize * 3) % 4) == 2) {
      k = k + 2;
    }
    if (((hsize * 3) % 4) == 3) {
      k = k + 1;
    }
    if (((hsize * 3) % 4) == 0) { // kokokara 3gyou nakutemoyoi
      ;
    }
  }
  printf("hsize=%d,vsize=%d\n", hsize, vsize);
  BIN = (unsigned char *)malloc(sizeof(unsigned char) * hsize * vsize);
  GIN = (unsigned char *)malloc(sizeof(unsigned char) * hsize * vsize);
  RIN = (unsigned char *)malloc(sizeof(unsigned char) * hsize * vsize);

  /*******************************************/
  /* */
  /* RGB PIXEL SITAKARA HAIRETU KARA UEKARA */
  /* */
  /* HAIRETU NI NAOSU */
  /* */
  /*******************************************/

  for (i = 0; i < vsize; i++) {
    for (j = 0; j < hsize; j++) {
      *(RIN + j * vsize + i) = *(R + j * vsize + vsize - 1 - i);
      *(GIN + j * vsize + i) = *(G + j * vsize + vsize - 1 - i);
      *(BIN + j * vsize + i) = *(B + j * vsize + vsize - 1 - i);
    }
  }
  free(R);
  free(G);
  free(B);
  // BOUT = (unsigned char *)malloc(sizeof(unsigned char)*hsize*vsize);
  // GOUT = (unsigned char *)malloc(sizeof(unsigned char)*hsize*vsize);
  // ROUT = (unsigned char *)malloc(sizeof(unsigned char)*hsize*vsize);
  PALETGAZOU = (unsigned char *)malloc(sizeof(unsigned char) * hsize * vsize);

  // resize(hsize, vsize, cnt->sizex, cnt->sizey, RIN, ROUT);
  // resize(hsize, vsize, cnt->sizex, cnt->sizey, GIN, GOUT);
  // resize(hsize, vsize, cnt->sizex, cnt->sizey, BIN, BOUT);
  /*cnt->edge = 3*/; /*cnt->et = 25;cnt->r = 3.0;cnt->g = 5.0;cnt->b = 1.0;*/
      /*cnt->gamma = 1.0;*/ /*cnt->omh = 0;*/
  cnt->div = 1;             // cnt->hasyori=20;
  MedianCut(hsize, vsize, RIN, GIN, BIN, PALETGAZOU, REDUCE_R, REDUCE_G,
            REDUCE_B, cnt->edge, cnt->et, cnt->edon, cnt->kmon, cnt->r, cnt->g,
            cnt->b, cnt->gamma, cnt->omh, cnt->hasyori, cnt->dither, cnt->cs,
            cnt->per, cnt->div, cnt->pw, cnt->pw2, cnt->t, cnt->a, cnt->rgamma,
            cnt->bgamma, cnt->ygamma, cnt->bun, cnt->lpfon, cnt->bai,
            cnt->lpfk); // ROUT, GOUT, BOUT);
  // debug start
  // for(;;);
  // debug end
  sprintf(cnt->OutputFileName, "%s.bmp", outputfile);
  fpw = fopen(cnt->OutputFileName, "wb");
  // hsize = cnt->sizex;
  // vsize = cnt->sizey;

  /*******************************************/
  /* */
  /* BITMAP HEADER 54BYTES WRITE */
  /* */
  /*******************************************/

  if ((hsize % 4) == 0) {
    filesize = 0x36 + 0x400 + hsize * vsize;
  } else {
    filesize = 0x36 + 0x400 + hsize * vsize + (4 - hsize % 4) * vsize;
  }

  temp = 0x42; // B
  fwrite(&temp, sizeof(unsigned char), 1, fpw);
  temp = 0x4D; // M
  fwrite(&temp, sizeof(unsigned char), 1, fpw);
  fwrite(&filesize, sizeof(int), 1, fpw);
  filesize = 0;
  fwrite(&filesize, sizeof(int), 1, fpw);
  // fwrite(&filesize,sizeof(int),1,fpw);
  filesize = 0x00000436;
  fwrite(&filesize, sizeof(int), 1, fpw);
  filesize = 0x28;
  fwrite(&filesize, sizeof(int), 1, fpw);
  fwrite(&hsize, sizeof(int), 1, fpw);
  fwrite(&vsize, sizeof(int), 1, fpw);
  filesize = 0x00080001;
  fwrite(&filesize, sizeof(int), 1, fpw);
  filesize = 0x0;
  fwrite(&filesize, sizeof(int), 1, fpw);

  if ((hsize % 4) == 0) {
    filesize = hsize * vsize;
  } else {
    filesize = hsize * vsize + (4 - hsize % 4) * vsize;
  }

  fwrite(&filesize, sizeof(int), 1, fpw);
  filesize = 0x0;
  fwrite(&filesize, sizeof(int), 1, fpw);
  fwrite(&filesize, sizeof(int), 1, fpw);
  filesize = 0x100;
  fwrite(&filesize, sizeof(int), 1, fpw);
  fwrite(&filesize, sizeof(int), 1, fpw);
  temp = 0;
  for (i = 0; i < 256; i++) {
    fwrite(REDUCE_B + i, sizeof(unsigned char), 1, fpw);
    fwrite(REDUCE_G + i, sizeof(unsigned char), 1, fpw);
    fwrite(REDUCE_R + i, sizeof(unsigned char), 1, fpw);
    fwrite(&temp, sizeof(unsigned char), 1, fpw);
  }

  for (i = 0; i < vsize; i++) {
    for (j = 0; j < hsize; j++) {
      fwrite(PALETGAZOU + j * vsize + vsize - 1 - i, sizeof(unsigned char), 1,
             fpw);
      // fwrite(GOUT+j*vsize+vsize-1-i,sizeof(unsigned char),1,fpw);
      // fwrite(ROUT+j*vsize+vsize-1-i,sizeof(unsigned char),1,fpw);
    }
    l = 4 - hsize % 4;
    if (l != 4) {
      for (k = 0; k < l; k++) {
        fwrite(&temp, sizeof(unsigned char), 1, fpw);
      }
    }
  }
  // free(ROUT);
  // free(GOUT);
  // free(BOUT);
  free(PALETGAZOU);
  fclose(fpw);
  return 0;
} // main owari
