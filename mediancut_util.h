#pragma once

#include <stdint.h>
#include <stdlib.h>
#include <tuple>
#include <vector>

// ピクセルの属性
//
// 0-255 はパレット番号。
// 0 未満を特別な値として意味を持たせる。
enum PixelType {
  kRegularPixel = -2, // !edge (旧index3 == 1)
  kEdgePixel = -1,    //  edge (旧index3 == -1)
  // 0-255 は !edge
};

void JacobiSetupNotEquals(int SIZE, const int16_t *PIXTYPES, int16_t PIXTYPE,
                          const double *YIN, const double *UIN,
                          const double *VIN,
                          const std::tuple<double, double, double> &U_RGB,
                          double A[3][3]);
void JacobiSetupEquals(int SIZE, const int16_t *PIXTYPES, int16_t PIXTYPE,
                       const double *YIN, const double *UIN, const double *VIN,
                       const std::tuple<double, double, double> &U_RGB,
                       double A[3][3]);
int Jacobi(int ct, double eps, double A[3][3], double A1[3][3], double A2[3][3],
           double X1[3][3], double X2[3][3]);
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
int DetectEdge0(int SIZE, int et, const int *HHEDGER, const int *VVEDGER,
                const int *HHEDGEG, const int *VVEDGEG, const int *HHEDGEB,
                const int *VVEDGEB, int16_t *OUT);
int DetectEdge1(int SIZE, int et, const int *HHEDGER, const int *VVEDGER,
                const int *HHEDGEG, const int *VVEDGEG, const int *HHEDGEB,
                const int *VVEDGEB, int16_t *OUT);
int DetectEdge2(int SIZE, const int *EDGER, const int *VEDGER, const int *EDGEG,
                const int *VEDGEG, const int *EDGEB, const int *VEDGEB,
                int16_t *OUT);
int DetectEdge3(int SIZE, int et, const int *HHEDGER, const int *VVEDGER,
                const int *HHEDGEG, const int *VVEDGEG, const int *HHEDGEB,
                const int *VVEDGEB, int16_t *OUT);
int DetectEdge4(int SIZE, const int *HHEDGER, const int *VVEDGER,
                const int *HHEDGEG, const int *VVEDGEG, const int *HHEDGEB,
                const int *VVEDGEB, int16_t *OUT, double *EDGERASISAY);
int DetectEdge5(int SIZE, int et, const int *HHEDGER, const int *VVEDGER,
                const int *HHEDGEG, const int *VVEDGEG, const int *HHEDGEB,
                const int *VVEDGEB, int16_t *OUT);
int DetectEdge6(int SIZE, int et, const int *HHEDGER, const int *VVEDGER,
                const int *HHEDGEG, const int *VVEDGEG, const int *HHEDGEB,
                const int *VVEDGEB, int16_t *OUT);
std::tuple<double, double, double> SumEquals(int SIZE, const int16_t *PIXTYPES,
                                             int16_t PIXTYPE, const double *v1,
                                             const double *v2,
                                             const double *v3);
std::tuple<double, double, double>
SumNotEquals(int SIZE, const int16_t *PIXTYPES, int16_t PIXTYPE,
             const double *v1, const double *v2, const double *v3);
int FindMinimumDistanceIndex(int SIZE, double v1, double v2, double v3,
                             const uint8_t *a1, const uint8_t *a2,
                             const uint8_t *a3);
int FindMinimumDistanceIndex(int SIZE, double v1, double v2, double v3,
                             const double *a1, const double *a2,
                             const double *a3);
void Dithering(int IMAGE_SIZE, const double *YIN, const double *UIN,
               const double *VIN, int PALET_SIZE, const double *Y_JYUSHIN,
               const double *U_JYUSHIN, const double *V_JYUSHIN, uint8_t *OUT);
void Dithering1(int hsize, int vsize, int dither, double per,
                const uint8_t *RIN, const uint8_t *GIN, const uint8_t *BIN,
                int PALET_SIZE, const uint8_t *REDUCE_R,
                const uint8_t *REDUCE_G, const uint8_t *REDUCE_B, uint8_t *OUT);
void Dithering2(int hsize, int vsize, int dither, double per, const double *VIN,
                const double *YIN, const double *UIN, int PALET_SIZE,
                const double *V_JYUSHIN1, const double *Y_JYUSHIN1,
                const double *U_JYUSHIN1, uint8_t *OUT);