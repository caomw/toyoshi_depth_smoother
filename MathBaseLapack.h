//
// Copyright (c) 2009-2011  Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
//  and the National Institute of Advanced Industrial Science and Technology
//
// $Id: MathBaseLapack.h 5744 2012-03-15 16:53:25Z shun $
//

#pragma once

// http://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/
#ifdef _DLL // /MD
#if defined(_M_IX86)
#pragma comment(lib, "mkl_intel_c_dll.lib")
#elif defined(_M_X64)
#pragma comment(lib, "mkl_intel_lp64_dll.lib")
#endif
#pragma comment(lib, "mkl_intel_thread_dll.lib")
#pragma comment(lib, "mkl_core_dll.lib")
#pragma comment(lib, "libiomp5md.lib")
#else // /MT
#if defined(_M_IX86)
#pragma comment(lib, "mkl_intel_c.lib")
#elif defined(_M_X64)
#pragma comment(lib, "mkl_intel_lp64.lib")
#endif
#pragma comment(lib, "mkl_intel_thread.lib")
#pragma comment(lib, "mkl_core.lib")
#pragma comment(lib, "libiomp5mt.lib")
#endif
// http://software.intel.com/sites/products/documentation/hpc/composerxe/en-us/cpp/mac/optaps/common/optaps_par_compat_libs_using.htm
// /nodefaultlib:vcomp

#undef small // conflict with "rpcndr.h" used in MFC
#include <mkl_lapack.h>
#include <mkl_blas.h>

#include "MathBaseUtil.h"
#include "MiscUtil.h"

#undef min
#undef max

namespace slib
{
namespace
{

//------------------------------------------------------------
// lower, non-transposed matrix in Rectangular Full Packed format
//------------------------------------------------------------
template <typename CMatrixType>
inline
void ConvertToRFP(const CMatrixType& lower, CDynamicMatrix<double>& rfp)
{
	int dim = lower.GetNumRows();
	if (dim != lower.GetNumCols())
	{
		throw std::logic_error("invalid matrix size in "  __FUNCTION__);
	}
	if (dim % 2)
	{
		// odd
		int col = lower.GetNumCols() / 2 + 1;
		rfp.Resize(dim, col);
		for (int c = 0; c < dim; c++)
		{
			for (int r = c; r < dim; r++)
			{
				if (c < col)
				{
					rfp(r, c) = lower(r, c);
				}
				else
				{
					rfp(c - col, r - col + 1) = lower(r, c);
				}
			}
		}
	}
	else
	{
		// even
		int col = lower.GetNumCols() / 2;
		rfp.Resize(dim + 1, col);
		for (int c = 0; c < dim; c++)
		{
			for (int r = c; r < dim; r++)
			{
				if (c < col)
				{
					rfp(r + 1, c) = lower(r, c);
				}
				else
				{
					rfp(c - col, r - col) = lower(r, c);
				}
			}
		}
	}
}

template <typename CMatrixType>
inline
void ConvertFromRFP(const CDynamicMatrix<double>& rfp, CMatrixType& lower)
{
	int col = rfp.GetNumCols();
	if (col % 2)
	{
		//even
		int dim = rfp.GetNumRows() - 1;
		if (dim != lower.GetNumRows() || col != lower.GetNumCols() / 2)
		{
			throw std::logic_error("invalid matrix size in " __FUNCTION__);
		}
		lower.Fill(0);


		for (int c = 0; c < dim; c++)
		{
			for (int r = c; r < dim; r++)
			{
				if (c < col)
				{
					lower(r, c) = rfp(r + 1, c);
				}
				else
				{
					lower(r, c) = rfp(c - col, r - col);
				}
			}
		}
	}
	else
	{
		// odd
		int dim = rfp.GetNumRows();
		if (dim != lower.GetNumRows() || 2 * col - 1 != lower.GetNumCols())
		{
			throw std::logic_error("invalid matrix size in" __FUNCTION__);
		}
		for (int c = 0; c < dim; c++)
		{
			for (int r = c; r < dim; r++)
			{
				if (c < col)
				{
					lower(r, c) = rfp(r, c);
				}
				else
				{
					lower(r, c) = rfp(c - col, r - col + 1);
				}
			}
		}
	}
}

} // unnamed namespace

//------------------------------------------------------------
// Cholesky factorization of a symmetric (Hermitian) positive-definite matrix
//------------------------------------------------------------

template <typename CMatrixType>
inline
void DPOTRF_driver(
    const CMatrixType& mat, // mxm
    CMatrixType& lower // lower triangular matrix
)
{
	const char uplo = 'L';
	const MKL_INT n = mat.GetNumRows();
	if (mat.GetNumCols() != n || lower.GetNumRows() != n || lower.GetNumCols() != n)
	{
		throw std::logic_error("invalid matrix size in " __FUNCTION__);
	}

	CDynamicMatrix<double> a(n, n, mat.ptr());
	MKL_INT lda = n;
	MKL_INT info;
	DPOTRF(&uplo, &n, a.ptr(), &lda, &info);
	if (info < 0)
	{
		throw std::runtime_error(string_format("%d-th parameter of DPFTRF had an illegal value in " __FUNCTION__, -info));
	}
	else if (info > 0)
	{
		throw std::runtime_error(string_format("the leading minor of order %d is not positive-definite in " __FUNCTION__, info));
	}
	for (int c = 0; c < n; c++)
	{
		for (int r = 0; r < n; r++)
		{
			if (r < c)
			{
				lower(r, c) = 0;
			}
			else
			{
				lower(r, c) = a(r, c);
			}
		}
	}
}

//------------------------------------------------------------
// QR Decomposition
//------------------------------------------------------------

template <typename CMatrixType>
inline
void DGELS_driver(
    const CMatrixType& mat, // mxn
    const CMatrixType& vec, // mxk
    CMatrixType& sol // nxk
)
{
	const char trans = 'N';
	const MKL_INT m = mat.GetNumRows();
	const MKL_INT n = mat.GetNumCols();
	const MKL_INT nrhs = vec.GetNumCols();
	if (vec.GetNumRows() != m || sol.GetNumRows() != n || sol.GetNumCols() != nrhs)
	{
		throw std::logic_error("invalid matrix size in " __FUNCTION__);
	}

	CDynamicMatrix<double> a(m, n, mat.ptr());
	const MKL_INT lda = m;
	CDynamicMatrix<double> b(m, nrhs, vec.ptr());
	const MKL_INT ldb = m;
	double chkwork;
	MKL_INT lwork = -1;
	MKL_INT info;
	DGELS(&trans, &m, &n, &nrhs, a.ptr(), &lda, b.ptr(), &ldb, &chkwork, &lwork, &info);
	lwork = chkwork;
	CDynamicVector<double> work(lwork);
	DGELS(&trans, &m, &n, &nrhs, a.ptr(), &lda, b.ptr(), &ldb, work.ptr(), &lwork, &info);
	if (info < 0)
	{
		throw std::runtime_error(string_format("%d-th parameter of DGELS had an illegal value in " __FUNCTION__, -info));
	}
	else if (info > 0)
	{
		throw std::runtime_error(string_format("%d-th diagonal element of the triangular factor is zero in " __FUNCTION__, info));
	}
	sol.Initialize(b.ptr());
}

//------------------------------------------------------------
// Singular Value Decomposition
//------------------------------------------------------------

template <typename CMatrixType, typename CVectorType>
inline
void DGESVD_driver(
    const CMatrixType& mat, // mxn
    CMatrixType& matU, // mxm
    CVectorType& vecW, // minmn
    CMatrixType& matVt // nxn
)
{
	const char jobu = 'A';
	const char jobvt = 'A';
	const MKL_INT m = mat.GetNumRows();
	const MKL_INT n = mat.GetNumCols();
	const MKL_INT minmn = std::min(m, n);
	if (matU.GetNumRows() != m || matU.GetNumCols() != m || vecW.GetNumRows() != minmn || matVt.GetNumRows() != n || matVt.GetNumCols() != n)
	{
		throw std::logic_error("invalid matrix size in " __FUNCTION__);
	}
	CDynamicMatrix<double> a(m, n, mat.ptr());
	const MKL_INT lda = m;
	CDynamicVector<double> s(minmn);
	CDynamicMatrix<double> u(m, m);
	const MKL_INT ldu = m;
	CDynamicMatrix<double> vt(n, n);
	const MKL_INT ldvt = n;
	double chkwork;
	MKL_INT lwork = -1;
	MKL_INT info;
	DGESVD(&jobu, &jobvt, &m, &n, a.ptr(), &lda, s.ptr(), u.ptr(), &ldu, vt.ptr(), &ldvt, &chkwork, &lwork, &info);
	lwork = chkwork;
	CDynamicVector<double> work(lwork);
	DGESVD(&jobu, &jobvt, &m, &n, a.ptr(), &lda, s.ptr(), u.ptr(), &ldu, vt.ptr(), &ldvt, work.ptr(), &lwork, &info);
	if (info < 0)
	{
		throw std::runtime_error(string_format("%d-th parameter of DGESVD had an illegal value in " __FUNCTION__, -info));
	}
	else if (info > 0)
	{
		throw std::runtime_error(string_format("%d superdiagonals of the intermediate bidiagonal form did not converge to zero in " __FUNCTION__, info));
	}
	vecW.Initialize(s.ptr());
	matU.Initialize(u.ptr());
	matVt.Initialize(vt.ptr());
}

template <typename CMatrixType, typename CVectorType>
inline
void DSYEV_driver(const CMatrixType& mat, CMatrixType& vec, CVectorType& value)
{
	char jobz = 'V';
	char uplo = 'U';
	MKL_INT n = mat.GetNumRows();
	if (mat.GetNumCols() != n || vec.GetNumRows() != n || vec.GetNumCols() != n || value.GetNumRows() != n)
	{
		throw std::logic_error("invalid matrix size in " __FUNCTION__);
	}
	CDynamicMatrix<double> a(n, n, mat.ptr());
	const MKL_INT lda = n;
	CDynamicVector<double> w(n);
	double chkwork;
	MKL_INT lwork = -1;
	MKL_INT info;
	DSYEV(&jobz, &uplo, &n, a.ptr(), &lda, w.ptr(), &chkwork, &lwork, &info);
	lwork = chkwork;
	CDynamicVector<double> work(lwork);
	DSYEV(&jobz, &uplo, &n, a.ptr(), &lda, w.ptr(), work.ptr(), &lwork, &info);
	if (info < 0)
	{
		throw std::runtime_error(string_format("%d-th parameter of DSYEV had an illegal value in " __FUNCTION__, -info));
	}
	else if (info > 0)
	{
		throw std::runtime_error(string_format("DSYEV failed to converge in " __FUNCTION__, info));
	}

	vec.Initialize(a.ptr());
	value.Initialize(w.ptr());
}

//------------------------------------------------------------
// Inverse
//------------------------------------------------------------

template <typename CMatrixType>
inline
void get_inverse(const CMatrixType& mat, CMatrixType& ret)
{
	const MKL_INT m = mat.GetNumRows();
	if (mat.GetNumCols() != m || ret.GetNumRows() != m || ret.GetNumCols() != m)
	{
		throw std::logic_error("invalid matrix size in " __FUNCTION__);
	}

	const MKL_INT n = m;
	CDynamicMatrix<double> a(m, n, mat.ptr());

	const MKL_INT lda = m;
	CDynamicVector<MKL_INT> ipiv(std::min(m, n));
	MKL_INT info;
	DGETRF(&m, &n, a.ptr(), &lda, ipiv.ptr(), &info);
	if (info < 0)
	{
		throw std::runtime_error(string_format("%d-th parameter of DGETRF had an illegal value in " __FUNCTION__, -info));
	}
	else if (info > 0)
	{
		throw std::runtime_error(string_format("%d-th diagonal element of the factor U is zero in " __FUNCTION__, info));
	}

	double chkwork;
	MKL_INT lwork = -1;
	DGETRI(&n, a.ptr(), &lda, ipiv.ptr(), &chkwork, &lwork, &info);
	lwork = chkwork;
	CDynamicVector<double> work(lwork);
	DGETRI(&n, a.ptr(), &lda, ipiv.ptr(), work.ptr(), &lwork, &info);
	if (info < 0)
	{
		throw std::runtime_error(string_format("%d-th parameter of DGETRI had an illegal value in " __FUNCTION__, -info));
	}
	else if (info > 0)
	{
		throw std::runtime_error(string_format("%d-th diagonal element of the factor U is zero in " __FUNCTION__, info));
	}

	ret.Initialize(a.ptr());
}

template <typename CMatrixType>
inline
const CMatrixType inverse_of(const CMatrixType& mat)
{
	CMatrixType ret;
	// set dimension (only for CDynamicMatrix)
	if (ret.GetNumRows() != mat.GetNumRows() || ret.GetNumCols() != mat.GetNumCols())
	{
		ret = mat;
	}
	get_inverse(mat, ret);
	return ret;
}

template <typename T>
inline
T determinant_of(const CDynamicMatrix<T>& mat)
{
	const MKL_INT m = mat.GetNumRows();
	if (m != mat.GetNumCols())
	{
		throw std::logic_error("invalid matrix size in " __FUNCTION__);
	}

	const MKL_INT n = m;
	CDynamicMatrix<double> a(m, n, mat.ptr());
	const MKL_INT lda = m;
	CDynamicVector<MKL_INT> ipiv(std::min(m, n));
	const MKL_INT info;
	DGETRF(&m, &n, a.ptr(), &lda, ipiv.ptr(), &info);
	if (info < 0)
	{
		throw std::runtime_error(string_format("%d-th parameter of DGETRF had an illegal value in " __FUNCTION__, -info));
	}
	else if (info > 0)
	{
		return 0;
	}

	T det = 1;
	for (int i = 0; i < minmn; i++)
	{
		if (ipiv[i] != i + 1)
		{
			det *= -a(i, i);
		}
		else
		{
			det *= a(i, i);
		}
	}
	return det;
}

template <int nDimension, typename T>
inline
void normalize_rotation(CMatrix<nDimension, nDimension, T>& mat)
{
	STATIC_ASSERT(nDimension == 3 || nDimension == 4);

	const char jobu = 'A';
	const char jobvt = 'A';
	const MKL_INT m = 3;
	const MKL_INT n = 3;
	const MKL_INT minmn = 3;
	float a[3 * 3];
	for (int c = 0; c < 3; c++)
	{
		for (int r = 0; r < 3; r++)
		{
			a[r + 3 * c] = mat(r, c);
		}
	}
	const MKL_INT lda = 3;
	float s[3];
	float u[3 * 3];
	const MKL_INT ldu = 3;
	float vt[3 * 3];
	const MKL_INT ldvt = 3;
	MKL_INT lwork = 15;
	float work[15];
	MKL_INT info;
	SGESVD(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);
	if (info != 0)
	{
		throw std::runtime_error("error in SGESVD");
	}
	CMatrix<3, 3, float> rot = CMatrix<3, 3, float>(u) * CMatrix<3, 3, float>(vt);
	for (int c = 0; c < 3; c++)
	{
		for (int r = 0; r < 3; r++)
		{
			mat(r, c) = rot(r, c);
		}
	}
}

}