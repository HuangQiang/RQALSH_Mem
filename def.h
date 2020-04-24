#ifndef __DEF_H
#define __DEF_H

// -----------------------------------------------------------------------------
//  Macros
// -----------------------------------------------------------------------------
#define MIN(a, b)	(((a) < (b)) ? (a) : (b))
#define MAX(a, b)	(((a) > (b)) ? (a) : (b))
#define SQR(x)		((x) * (x))
#define SUM(x, y)	((x) + (y))
#define DIFF(x, y)	((y) - (x))
#define SWAP(x, y)	{int tmp=x; x=y; y=tmp;}

// -----------------------------------------------------------------------------
//  Constants
// -----------------------------------------------------------------------------
const int   TOPK[]        = { 1, 2, 5, 10 };
const int   MAX_ROUND     = 4;
const int   MAXK          = TOPK[MAX_ROUND - 1];

const int   CANDIDATES    = 100;
const int   N_THRESHOLD   = (CANDIDATES + MAXK) * 2;
const int   SCAN_SIZE     = 64;
const int   MAX_BLOCK_NUM = 10000;
const int   MAGIC         = 36553368;
const float LAMBDA        = 0.9f;

const int   SIZEBOOL      = (int) sizeof(bool);
const int   SIZECHAR      = (int) sizeof(char);
const int   SIZEINT       = (int) sizeof(int);
const int   SIZEFLOAT     = (int) sizeof(float);
const int   SIZEDOUBLE    = (int) sizeof(double);

const int   MAXINT        = 2147483647;
const int   MININT        = -MAXINT;
const float MAXREAL       = 3.402823466e+38F;
const float MINREAL       = -MAXREAL;

const float E             = 2.7182818F;
const float PI            = 3.141592654F;
const float FLOATZERO     = 1e-6F;
const float ANGLE         = PI / 8.0f;

#endif // __DEF_H
