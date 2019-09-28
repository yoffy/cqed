// ./*** -i 入力24bitBMPファイル名 -o 出力bmpファイル名（.bmpはいらない）

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void MedianCut(int hsize, int vsize, unsigned char *RIN, unsigned char *GIN,
               unsigned char *BIN, unsigned char *PALETGAZOU,
               unsigned char REDUCE_R[256], unsigned char REDUCE_G[256],
               unsigned char REDUCE_B[256], int vvv, int et, int edon, int kmon,
               double rmult, double gmult, double bmult, double gamma, int omh,
               int hasyori, int dither, int cs, double per, int div, double pw,
               double pw2, double IRO_THIETA, double aaaa, double rgamma,
               double bgamma, double ygamma, int bun, int lpf, double nbai,
               double lpfk); // unsigned char *ROUT, unsigned char *GOUT,
                            // unsigned char *BOUT)

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
  int i, j, k, l;
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
