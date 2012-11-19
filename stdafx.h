#include <cstdio>
#include <cstdarg>

//#include <direct.h>

#include <string>

//12/11/16 Masahiko MKL��INT�͎g�킸�ɁAclapack�Œ�`����Ă����g���܂��B
#define MKL_INT __CLPK_integer

//MKL��LAPACK�̊֐���clapack��LAPACK�֐��ɒu�������܂�
#define SGESVD sgesvd_
#define DGETRI dgetri_
#define DGETRF dgetrf_

#include "MathBaseLapack.h"
#include "MiscUtil.h"
#include "BmpUtil.h"

//12/11/13�@������mac��Accelerate���g���܂���
#include <Accelerate/Accelerate.h>
#include <vecLib/vecLib.h>
#include <assert.h>
