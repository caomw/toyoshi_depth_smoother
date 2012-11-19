#include <cstdio>
#include <cstdarg>

//#include <direct.h>

#include <string>

//12/11/16 Masahiko MKLのINTは使わずに、clapackで定義されてるやつを使います。
#define MKL_INT __CLPK_integer

//MKLのLAPACKの関数をclapackのLAPACK関数に置き換えます
#define SGESVD sgesvd_
#define DGETRI dgetri_
#define DGETRF dgetrf_

#include "MathBaseLapack.h"
#include "MiscUtil.h"
#include "BmpUtil.h"

//12/11/13　かわりにmacのAccelerateを使いまする
#include <Accelerate/Accelerate.h>
#include <vecLib/vecLib.h>
#include <assert.h>
