//
// Copyright (c) 2009-2011 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology
//
// $Id: MathBase.h 5744 2012-03-15 16:53:25Z shun $
//

#pragma once

#ifndef TRACE
#define TRACE printf
#endif

#ifndef _ASSERTE
#define _ASSERTE assert
#endif

#undef min
#undef max

namespace slib
{
#define STATIC_ASSERT(x) static_assert(x, "static_assert() failed")

//--------------------------------------------------------------------
// static matrix
//--------------------------------------------------------------------

template <int nNumRows, int nNumCols, typename T>
class CMatrix
{
protected:
	T m[nNumRows *nNumCols]; // row major order

public:
	// constructors
	CMatrix(void)
	{
		for (int i = 0; i < nNumRows * nNumCols; i++)
		{
			m[i] = T();
		}
	}

	explicit CMatrix(T const *const e)
	{
		Initialize(e);
	}

	// copy
	CMatrix(const CMatrix& mat)
	{
		for (int i = 0; i < nNumRows * nNumCols; i++)
		{
			m[i] = mat.m[i];
		}
	}

	const CMatrix& operator =(const CMatrix& mat)
	{
		for (int i = 0; i < nNumRows * nNumCols; i++)
		{
			m[i] = mat.m[i];
		}
		return *this;
	}

	// coercion
	template <typename T2>
	CMatrix(const CMatrix<nNumRows, nNumCols, T2>& mat)
	{
		for (int i = 0; i < nNumRows * nNumCols; i++)
		{
			m[i] = mat.ptr()[i];
		}
	}

	template <typename T2>
	const CMatrix& operator =(const CMatrix<nNumRows, nNumCols, T2>& mat)
	{
		for (int i = 0; i < nNumRows * nNumCols; i++)
		{
			m[i] = mat.ptr()[i];
		}
		return *this;
	}

	// helpers
	template <typename T2>
	const CMatrix& Initialize(T2 const *const element)
	{
		_ASSERTE(element);
		for (int i = 0; i < nNumRows * nNumCols; i++)
		{
			m[i] = element[i];
		}
		return *this;
	}

	template <typename T2>
	void Fill(T2 const e)
	{
		for (int i = 0; i < nNumRows * nNumCols; i++)
		{
			m[i] = e;
		}
	}

	//---------- accessors
	T& operator()(const int r, const int c)
	{
		_ASSERTE(0 <= r && r < nNumRows && 0 <= c && c < nNumCols);
		return m[r + nNumRows * c];
	}

	T const& operator()(const int r, const int c) const
	{
		_ASSERTE(0 <= r && r < nNumRows && 0 <= c && c < nNumCols);
		return m[r + nNumRows * c];
	}

	T const *ptr(void) const
	{
		return m;
	}

	T *ptr(void)
	{
		return m;
	}

	//---------- arithmetic operators
	template <typename T2>
	const CMatrix& operator +=(const CMatrix<nNumRows, nNumCols, T2>& mat)
	{
		for (int i = 0; i < nNumRows * nNumCols; i++)
		{
			m[i] += mat.ptr()[i];
		}
		return *this;
	}

	template <typename T2>
	const CMatrix& operator -=(const CMatrix<nNumRows, nNumCols, T2>& mat)
	{
		for (int i = 0; i < nNumRows * nNumCols; i++)
		{
			m[i] -= mat.ptr()[i];
		}
		return *this;
	}

	template <typename T2>
	const CMatrix& operator *=(const CMatrix<nNumCols, nNumCols, T2>& mat)
	{
		return *this = *this * mat;
	}

	// "s" may be self-referencing
	const CMatrix& operator *=(double s)
	{
		for (int i = 0; i < nNumRows * nNumCols; i++)
		{
			m[i] *= s;
		}
		return *this;
	}

	// "s" may be self-referencing
	const CMatrix& operator /=(double s)
	{
		_ASSERTE(s != 0);
		for (int i = 0; i < nNumRows * nNumCols; i++)
		{
			m[i] /= s;
		}
		return *this;
	}

	// unary operator
	CMatrix operator -(void) const
	{
		return *this * (-1);
	}

	const CMatrix& operator +(void) const
	{
		return *this;
	}

	// binary operator
	CMatrix operator +(const CMatrix& v2) const
	{
		return CMatrix(*this) += v2;
	}

	CMatrix operator -(const CMatrix& v2) const
	{
		return CMatrix(*this) -= v2;
	}

	template <int nNumCols2, typename T2>
	CMatrix<nNumRows, nNumCols2, T> operator *(const CMatrix<nNumCols, nNumCols2, T2>& mat2) const
	{
		CMatrix<nNumRows, nNumCols2, T> iret;
		for (int c = 0; c < nNumCols2; c++)
		{
			for (int r = 0; r < nNumRows; r++)
			{
				iret(r, c) = T(0);
				for (int i = 0; i < nNumCols; i++)
				{
					iret(r, c) += (*this)(r, i) * mat2(i, c);
				}
			}
		}
		return iret;
	}

	CMatrix operator *(double s) const
	{
		return CMatrix(*this) *= s;
	}

	CMatrix operator /(double s) const
	{
		return CMatrix(*this) /= s;
	}

	friend const CMatrix operator *(double s, const CMatrix mat)
	{
		return mat * s;
	}

	// relationship operators
	bool operator ==(const CMatrix& mat2) const
	{
		return !memcmp(m, mat2.m, sizeof(T) * nNumRows * nNumCols);
	}

	bool operator !=(const CMatrix& mat2) const
	{
		return !(*this == mat2);
	}

	// helper methods
	CMatrix<nNumRows, 1, T> GetColumn(int c) const
	{
		CMatrix<nNumRows, 1, T> iret;
		for (int r = 0; r < nNumRows; r++)
		{
			iret(r, 0) = (*this)(r, c);
		}
		return iret;
	}

	CMatrix<1, nNumCols, T> GetRow(int r) const
	{
		CMatrix<1, nNumCols, T> iret;
		for (int c = 0; c < nNumCols; c++)
		{
			iret(0, c) = (*this)(r, c);
		}
		return iret;
	}

	template <int nSubRows, int nSubCols>
	CMatrix<nSubRows, nSubCols, T> GetSubMatrix(int sr, int sc) const
	{
		STATIC_ASSERT(nSubRows <= nNumRows);
		STATIC_ASSERT(nSubCols <= nNumCols);
		CMatrix<nSubRows, nSubCols, T> iret;
		for (int c = 0; c < nSubCols; c++)
		{
			for (int r = 0; r < nSubRows; r++)
			{
				iret(r, c) = (*this)(sr + r, sc + c);
			}
		}
		return iret;
	}

	template <int nNewCols, typename T2>
	CMatrix < nNumRows, nNumCols + nNewCols, T > AppendCols(const CMatrix<nNumRows, nNewCols, T2>& mat) const
	{
		CMatrix < nNumRows, nNumCols + nNewCols, T > iret;
		T *p = iret.ptr();

		// left half
		int nelem1 = nNumRows * nNumCols;
		for (int i = 0; i < nelem1; i++)
		{
			p[i] = m[i];
		}

		// right half
		int nelem2 = nNumRows * nNewCols;
		for (int i = 0; i < nelem2; i++)
		{
			p[i + nelem1] = mat.ptr()[i];
		}
		return iret;
	}

	template <int nNewRows>
	CMatrix < nNumRows + nNewRows, nNumCols, T > AppendRows(const CMatrix<nNewRows, nNumCols, T>& mat) const
	{
		CMatrix < nNumRows + nNewRows, nNumCols, T > iret;
		for (int c = 0; c < nNumCols; c++)
		{
			int r = 0;
			for (; r < nNumRows; r++)
			{
				iret(r, c) = (*this)(r, c);
			}
			for (int r2 = 0; r < nNumRows + nNewRows; r++, r2++)
			{
				iret(r, c) = mat(r2, c);
			}
		}
		return iret;
	}

	T min(void) const
	{
		T iret = m[0];
		for (int i = 1; i < nNumRows * nNumCols; i++)
		{
			iret = std::min(iret, m[i]);
		}
		return iret;
	}

	T max(void) const
	{
		T iret = m[0];
		for (int i = 1; i < nNumRows * nNumCols; i++)
		{
			iret = std::max(iret, m[i]);
		}
		return iret;
	}

	CMatrix abs(void) const
	{
		CMatrix result;
		for (int i = 0; i < nNumRows * nNumCols; i++)
		{
			result.m[i] = std::abs(m[i]);
		}
		return result;
	}

	// Note:
	// - static const variable cannot be used in the dynamic version.
	// - non-static method cannot be used in template variables.
	static int GetNumRows(void)
	{
		return nNumRows;
	}

	static int GetNumCols(void)
	{
		return nNumCols;
	}

	void Read(const std::string& filename)
	{
		TRACE("matrix <= %s\n", filename.c_str());
		FILE *fr = fopen(filename.c_str(), "rb");
		if (!fr)
		{
			throw std::runtime_error("failed to open file in " __FUNCTION__);
		}

		int nrows, ncols;
		fscanf(fr, "%d%d", &nrows, &ncols);
		if (nrows != nNumRows || ncols != nNumCols)
		{
			throw std::runtime_error("invalid matrix size in " __FUNCTION__);
		}

		for (int r = 0; r < nNumRows; r++)
		{
			for (int c = 0; c < nNumCols; c++)
			{
				double el;
				if (fscanf(fr, "%lf", &el) != 1)
				{
					fclose(fr);
					throw std::runtime_error("failed to parse file in " __FUNCTION__);
				}
				(*this)(r, c) = el;
			}
		}

		fclose(fr);
	}

	void Write(const std::string& filename) const
	{
		TRACE("matrix => %s\n", filename.c_str());
		FILE *fw = fopen(filename.c_str(), "wb");
		if (!fw)
		{
			throw std::runtime_error("failed to open file in " __FUNCTION__);
		}

		fprintf(fw, "%d %d\n", nNumRows, nNumCols);
		for (int r = 0; r < nNumRows; r++)
		{
			for (int c = 0; c < nNumCols; c++)
			{
				fprintf(fw, "% le ", double((*this)(r, c)));
			}
			fprintf(fw, "\n");
		}

		fclose(fw);
	}
};

// construction

template <typename T>
inline
CMatrix<2, 2, T> make_matrix(
    T const& a00, T const& a01,
    T const& a10, T const& a11)
{
	T element[] =
	{
		a00, a10,
		a01, a11,
	};
	return CMatrix<2, 2, T>(element);
}

template <typename T>
inline
CMatrix<3, 3, T> make_matrix(
    T const& a00, T const& a01, T const& a02,
    T const& a10, T const& a11, T const& a12,
    T const& a20, T const& a21, T const& a22)
{
	T element[] =
	{
		a00, a10, a20,
		a01, a11, a21,
		a02, a12, a22,
	};
	return CMatrix<3, 3, T>(element);
}

template <typename T>
inline
CMatrix<4, 4, T> make_matrix(
    T const& a00, T const& a01, T const& a02, T const& a03,
    T const& a10, T const& a11, T const& a12, T const& a13,
    T const& a20, T const& a21, T const& a22, T const& a23,
    T const& a30, T const& a31, T const& a32, T const& a33)
{
	T element[] =
	{
		a00, a10, a20, a30,
		a01, a11, a21, a31,
		a02, a12, a22, a32,
		a03, a13, a23, a33,
	};
	return CMatrix<4, 4, T>(element);
}

template <int nDim, typename T>
inline
CMatrix<nDim, nDim, T> make_diagonal_matrix(const CMatrix<nDim, 1, T>& v)
{
	CMatrix<nDim, nDim, T> iret;
	iret.Fill(0);
	for (int i = 0; i < nDim; i++)
	{
		iret(i, i) = v(i, 0);
	}
	return iret;
}

template <typename T>
inline
CMatrix<2, 2, T> make_diagonal_matrix(const T& e1, const T& e2)
{
	return make_diagonal_matrix(make_vector(e1, e2));
}

template <typename T>
inline
CMatrix<3, 3, T> make_diagonal_matrix(const T& e1, const T& e2, const T& e3)
{
	return make_diagonal_matrix(make_vector(e1, e2, e3));
}

template <typename T>
inline
CMatrix<4, 4, T> make_diagonal_matrix(const T& e1, const T& e2, const T& e3, const T& e4)
{
	return make_diagonal_matrix(make_vector(e1, e2, e3, e4));
}

// helper functions

template <int nNumRows, typename T>
inline
T trace_of(const CMatrix<nNumRows, nNumRows, T>& mat)
{
	T t = T();
	for (int i = 0; i < nNumRows; i++)
	{
		t += mat(i, i);
	}
	return t;
}

template <int nNumRows, int nNumCols, typename T>
inline
CMatrix<nNumCols, nNumRows, T> transpose_of(const CMatrix<nNumRows, nNumCols, T>& mat)
{
	CMatrix<nNumCols, nNumRows, T> iret;
	for (int c = 0; c < nNumCols; c++)
	{
		for (int r = 0; r < nNumRows; r++)
		{
			iret(c, r) = mat(r, c);
		}
	}
	return iret;
}

template <typename T>
inline
T determinant_of(const CMatrix<2, 2, T>& a)
{
	return a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0);
}

template <typename T>
inline
T determinant_of(const CMatrix<3, 3, T>& a)
{
	return a(0, 0) * a(1, 1) * a(2, 2) +
	       a(0, 1) * a(1, 2) * a(2, 0) +
	       a(0, 2) * a(1, 0) * a(2, 1) -
	       a(0, 0) * a(1, 2) * a(2, 1) -
	       a(0, 1) * a(1, 0) * a(2, 2) -
	       a(0, 2) * a(1, 1) * a(2, 0);
}

template <typename T>
inline
T determinant_of(const CMatrix<4, 4, T>& a)
{
	return ((a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0)) * a(2, 2) + (a(0, 2) * a(1, 0) - a(0, 0) * a(1, 2)) * a(2, 1) + (a(0, 1) * a(1, 2) - a(0, 2) * a(1, 1)) * a(2, 0)) * a(3, 3) +
	       ((a(0, 1) * a(1, 0) - a(0, 0) * a(1, 1)) * a(2, 3) + (a(0, 0) * a(1, 3) - a(0, 3) * a(1, 0)) * a(2, 1) + (a(0, 3) * a(1, 1) - a(0, 1) * a(1, 3)) * a(2, 0)) * a(3, 2) +
	       ((a(0, 0) * a(1, 2) - a(0, 2) * a(1, 0)) * a(2, 3) + (a(0, 3) * a(1, 0) - a(0, 0) * a(1, 3)) * a(2, 2) + (a(0, 2) * a(1, 3) - a(0, 3) * a(1, 2)) * a(2, 0)) * a(3, 1) +
	       ((a(0, 2) * a(1, 1) - a(0, 1) * a(1, 2)) * a(2, 3) + (a(0, 1) * a(1, 3) - a(0, 3) * a(1, 1)) * a(2, 2) + (a(0, 3) * a(1, 2) - a(0, 2) * a(1, 3)) * a(2, 1)) * a(3, 0);
}

template <typename T>
inline
CMatrix<2, 2, T> inverse_of(const CMatrix<2, 2, T>& a)
{
	return make_matrix(a(1, 1), -a(0, 1), -a(1, 0), a(0, 0)) /
	       determinant_of(a);
}

template <typename T>
inline
CMatrix<3, 3, T> inverse_of(const CMatrix<3, 3, T>& a)
{
	T det = determinant_of(a);
	return make_matrix(
	           (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)), (a(0, 2) * a(2, 1) - a(0, 1) * a(2, 2)), (a(0, 1) * a(1, 2) - a(0, 2) * a(1, 1)),
	           (a(1, 2) * a(2, 0) - a(1, 0) * a(2, 2)), (a(0, 0) * a(2, 2) - a(0, 2) * a(2, 0)), (a(0, 2) * a(1, 0) - a(0, 0) * a(1, 2)),
	           (a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0)), (a(0, 1) * a(2, 0) - a(0, 0) * a(2, 1)), (a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0))
	       ) / det;
}

template <typename T>
inline
CMatrix<4, 4, T> inverse_of(const CMatrix<4, 4, T>& a)
{
	T det = determinant_of(a);
	return make_matrix(
	           ((a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)) * a(3, 3) + (a(1, 3) * a(2, 1) - a(1, 1) * a(2, 3)) * a(3, 2) + (a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2)) * a(3, 1)),
	           -((a(0, 1) * a(2, 2) - a(0, 2) * a(2, 1)) * a(3, 3) + (a(0, 3) * a(2, 1) - a(0, 1) * a(2, 3)) * a(3, 2) + (a(0, 2) * a(2, 3) - a(0, 3) * a(2, 2)) * a(3, 1)),
	           ((a(0, 1) * a(1, 2) - a(0, 2) * a(1, 1)) * a(3, 3) + (a(0, 3) * a(1, 1) - a(0, 1) * a(1, 3)) * a(3, 2) + (a(0, 2) * a(1, 3) - a(0, 3) * a(1, 2)) * a(3, 1)),
	           -((a(0, 1) * a(1, 2) - a(0, 2) * a(1, 1)) * a(2, 3) + (a(0, 3) * a(1, 1) - a(0, 1) * a(1, 3)) * a(2, 2) + (a(0, 2) * a(1, 3) - a(0, 3) * a(1, 2)) * a(2, 1)),
	           -((a(1, 0) * a(2, 2) - a(1, 2) * a(2, 0)) * a(3, 3) + (a(1, 3) * a(2, 0) - a(1, 0) * a(2, 3)) * a(3, 2) + (a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2)) * a(3, 0)),
	           ((a(0, 0) * a(2, 2) - a(0, 2) * a(2, 0)) * a(3, 3) + (a(0, 3) * a(2, 0) - a(0, 0) * a(2, 3)) * a(3, 2) + (a(0, 2) * a(2, 3) - a(0, 3) * a(2, 2)) * a(3, 0)),
	           -((a(0, 0) * a(1, 2) - a(0, 2) * a(1, 0)) * a(3, 3) + (a(0, 3) * a(1, 0) - a(0, 0) * a(1, 3)) * a(3, 2) + (a(0, 2) * a(1, 3) - a(0, 3) * a(1, 2)) * a(3, 0)),
	           ((a(0, 0) * a(1, 2) - a(0, 2) * a(1, 0)) * a(2, 3) + (a(0, 3) * a(1, 0) - a(0, 0) * a(1, 3)) * a(2, 2) + (a(0, 2) * a(1, 3) - a(0, 3) * a(1, 2)) * a(2, 0)),
	           ((a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0)) * a(3, 3) + (a(1, 3) * a(2, 0) - a(1, 0) * a(2, 3)) * a(3, 1) + (a(1, 1) * a(2, 3) - a(1, 3) * a(2, 1)) * a(3, 0)),
	           -((a(0, 0) * a(2, 1) - a(0, 1) * a(2, 0)) * a(3, 3) + (a(0, 3) * a(2, 0) - a(0, 0) * a(2, 3)) * a(3, 1) + (a(0, 1) * a(2, 3) - a(0, 3) * a(2, 1)) * a(3, 0)),
	           ((a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0)) * a(3, 3) + (a(0, 3) * a(1, 0) - a(0, 0) * a(1, 3)) * a(3, 1) + (a(0, 1) * a(1, 3) - a(0, 3) * a(1, 1)) * a(3, 0)),
	           -((a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0)) * a(2, 3) + (a(0, 3) * a(1, 0) - a(0, 0) * a(1, 3)) * a(2, 1) + (a(0, 1) * a(1, 3) - a(0, 3) * a(1, 1)) * a(2, 0)),
	           -((a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0)) * a(3, 2) + (a(1, 2) * a(2, 0) - a(1, 0) * a(2, 2)) * a(3, 1) + (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)) * a(3, 0)),
	           ((a(0, 0) * a(2, 1) - a(0, 1) * a(2, 0)) * a(3, 2) + (a(0, 2) * a(2, 0) - a(0, 0) * a(2, 2)) * a(3, 1) + (a(0, 1) * a(2, 2) - a(0, 2) * a(2, 1)) * a(3, 0)),
	           -((a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0)) * a(3, 2) + (a(0, 2) * a(1, 0) - a(0, 0) * a(1, 2)) * a(3, 1) + (a(0, 1) * a(1, 2) - a(0, 2) * a(1, 1)) * a(3, 0)),
	           ((a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0)) * a(2, 2) + (a(0, 2) * a(1, 0) - a(0, 0) * a(1, 2)) * a(2, 1) + (a(0, 1) * a(1, 2) - a(0, 2) * a(1, 1)) * a(2, 0))
	       ) / det;
}

#define DUMP_VECMAT(mat) \
	{\
		TRACE("["__FUNCTION__"] " #mat " =:\n"); \
		for (int i__=0; i__<(mat).GetNumRows(); i__++) \
		{ \
			TRACE("\t["); \
			for (int j__=0; j__<(mat).GetNumCols(); j__++) \
				TRACE(" % g", double((mat)(i__,j__))); \
			TRACE(" ]\n"); \
		} \
	}

//----------------------------------------------------------------------
// static vector
//----------------------------------------------------------------------

template <int nDimension, typename T>
class CVector : public CMatrix<nDimension, 1, T>
{
public:
	// initializers
	CVector(void) : CMatrix()
	{
	}

	explicit CVector(T const *const element)
		: CMatrix(element)
	{
	}

	// copy
	CVector(const CVector& vec)
		: CMatrix(vec)
	{
	}

	const CVector& operator =(const CVector& vec)
	{
		CMatrix::operator=(vec);
		return *this;
	}

	// coersion
	template <typename T2>
	CVector(const CMatrix<nDimension, 1, T2>& vec)
		: CMatrix(vec)
	{
	}

	template <typename T2>
	const CVector& operator =(const CMatrix<nDimension, 1, T2>& vec)
	{
		CMatrix::operator=(vec);
		return *this;
	}

	// operators
	T& operator [](const int n)
	{
		_ASSERTE(n < nDimension);
		return m[n];
	}

	T const& operator [](const int n) const
	{
		_ASSERTE(n < nDimension);
		return m[n];
	}
};

// construction

template <typename T>
inline
CVector<2, T> make_vector(T const& e1, T const& e2)
{
	T element[] = { e1, e2 };
	return CVector<2, T>(element);
}

template <typename T>
inline
CVector<3, T> make_vector(T const& e1, T const& e2, T const& e3)
{
	T element[] = { e1, e2, e3 };
	return CVector<3, T>(element);
}

template <typename T>
inline
CVector<4, T> make_vector(T const& e1, T const& e2, T const& e3, T const& e4)
{
	T element[] = { e1, e2, e3, e4 };
	return CVector<4, T>(element);
}

// vector utils

template <int nDimension, typename T>
inline
T norm_infty_of(const CMatrix<nDimension, 1, T>& vec)
{
	T len = abs(vec[0]);
	for (int i = 1; i < nDimension; i++)
	{
		len = std::max(len, abs(vec[i]));
	}
	return len;
}

template <int nDimension, typename T>
inline
T norm1_of(const CMatrix<nDimension, 1, T>& vec)
{
	T len = abs(vec[0]);
	for (int i = 1; i < nDimension; i++)
	{
		len += abs(vec[i]);
	}
	return len;
}

template <int nDimension, typename T>
inline
double norm2_of(const CMatrix<nDimension, 1, T>& vec)
{
	return sqrt(static_cast<double>(dot(vec, vec)));
}

template <int nDimension, typename T>
inline
T norm2_squared_of(const CMatrix<nDimension, 1, T>& vec)
{
	return dot(vec, vec);
}

template <int nDimension, typename T>
inline
CVector<nDimension, T> normalized_of(const CMatrix<nDimension, 1, T>& vec)
{
	T len = norm2_of(vec);
	_ASSERTE(len != 0);
	return vec / len;
}

template <int nDimension, typename T>
inline
T dot(const CMatrix<nDimension, 1, T>& vec0, const CMatrix<nDimension, 1, T>& vec1)
{
	T iret = T(0);
	for (int i = 0; i < nDimension; i++)
	{
		iret += vec0(i, 0) * vec1(i, 0);
	}
	return iret;
}

template <typename T>
inline
CVector<3, T> cross(const CMatrix<3, 1, T>& vec0, const CMatrix<3, 1, T>& vec1)
{
	return make_vector(
	           vec0(1, 0) * vec1(2, 0) - vec0(2, 0) * vec1(1, 0),
	           vec0(2, 0) * vec1(0, 0) - vec0(0, 0) * vec1(2, 0),
	           vec0(0, 0) * vec1(1, 0) - vec0(1, 0) * vec1(0, 0));
}

template <typename T>
inline
T crossZ(const CMatrix<2, 1, T>& vec0, const CMatrix<2, 1, T>& vec1)
{
	return vec0(0, 0) * vec1(1, 0) - vec0(1, 0) * vec1(0, 0);
}

template <typename T>
inline
CMatrix<3, 3, T> skew_symmetric_of(CMatrix<3, 1, T>& vec)
{
	return make_matrix<T>(
	           0, -vec(2, 0), vec(1, 0),
	           vec(2, 0), 0, -vec(0, 0),
	           -vec(1, 0), vec(0, 0), 0
	       );
}

template <int nDimension, typename T>
inline
CVector < nDimension + 1, T > homogeneous_of(const CMatrix<nDimension, 1, T>& vec)
{
	CVector < nDimension + 1, T > iret;
	for (int i = 0; i < nDimension; i++)
	{
		iret[i] = vec(i, 0);
	}
	iret[nDimension] = 1;
	return iret;
}

template <int nDimension, typename T>
inline
CVector < nDimension - 1, T > euclidean_of(const CMatrix<nDimension, 1, T>& vec)
{
	CVector < nDimension - 1, T > iret;
	for (int i = 0; i < nDimension - 1; i++)
	{
		iret[i] = vec(i, 0) / vec(nDimension - 1, 0);
	}
	return iret;
}

////////////////////////////////////////////////////////////////////////
// linear transformation
////////////////////////////////////////////////////////////////////////

// rigid transformation

template <typename T>
inline
CMatrix<4, 4, T> rigid_inverse_of(const CMatrix<4, 4, T>& a)
{
	CMatrix<4, 4, T> inv_rot =
	    make_matrix<T>(
	        a(0, 0), a(1, 0), a(2, 0), 0,
	        a(0, 1), a(1, 1), a(2, 1), 0,
	        a(0, 2), a(1, 2), a(2, 2), 0,
	        0, 0, 0, 1);
	CVector<4, T> inv_trans = make_vector<T>(-a(0, 3), -a(1, 3), -a(2, 3), 1);
	inv_trans = inv_rot * inv_trans;
	inv_rot(0, 3) = inv_trans [0];
	inv_rot(1, 3) = inv_trans [1];
	inv_rot(2, 3) = inv_trans [2];
	return inv_rot;
}

// translation

template <typename T>
inline
CVector<3, T> make_translation_vector(const CMatrix<4, 4, T>& mat)
//CVector<3, T> GetTranslationVector(const CMatrix<4, 4, T>& mat)
{
	return make_vector(mat(0, 3), mat(1, 3), mat(2, 3));
}

template <typename T>
inline
CMatrix<4, 4, T> make_translation_matrix(const CMatrix<4, 4, T>& mat)
//CMatrix<4, 4, T> GetTranslationMatrix(const CMatrix<4, 4, T>& mat)
{
	return make_translation_matrix(mat(0, 3), mat(1, 3), mat(2, 3));
}

template <typename T>
inline
CMatrix<4, 4, T> make_translation_matrix(const CMatrix<3, 1, T>& vec)
//CMatrix<4, 4, T> GetTranslationMatrix(const CMatrix<3, 1, T>& vec)
{
	return make_translation_matrix(vec(0, 0), vec(1, 0), vec(2, 0));
}

template <typename T>
inline
CMatrix<4, 4, T> make_translation_matrix(T const& tx, T const& ty, T const& tz)
//CMatrix<4, 4, T> GetTranslationMatrix(T const& tx, T const& ty, T const& tz)
{
	CMatrix<4, 4, T> mat = make_diagonal_matrix<T>(1, 1, 1, 1);
	mat(0, 3) = tx;
	mat(1, 3) = ty;
	mat(2, 3) = tz;
	return mat;
}

// rotation

template <int nDimension, typename T>
inline
void GetEulerAngles(const CMatrix<nDimension, nDimension, T>& mat, T& a, T& b, T& c)
{
	STATIC_ASSERT(nDimension == 3 || nDimension == 4);

	// assume b >= 0
	if (abs(mat(2, 2)) != 1)
	{
		a = atan2(mat(1, 2), mat(0, 2));
		b = atan2((T)sqrt(mat(0, 2) * mat(0, 2) + mat(1, 2) * mat(1, 2)), mat(2, 2));
		c = atan2(mat(2, 1), -mat(2, 0));
	}
	else
	{
		// assume c = 0
		a = atan2(mat(1, 0), mat(1, 1));
		b = 0;
		c = 0;
	}
}

template <int nDimension, typename T>
inline
void GetRollPitchYawAngles(const CMatrix<nDimension, nDimension, T>& mat, T& rollRadian, T& pitchRadian, T& yawRadian)
{
	STATIC_ASSERT(nDimension == 3 || nDimension == 4);

	// assume pitchRadian >= 0
	T cosPitch = sqrt(mat(0, 0) * mat(0, 0) + mat(1, 0) * mat(1, 0));
	if (cosPitch)
	{
		rollRadian = atan2(mat(1, 0), mat(0, 0));
		pitchRadian = atan2(-mat(2, 0), cosPitch);
		yawRadian = atan2(mat(2, 1), mat(2, 2));
	}
	else
	{
		rollRadian = -atan2(mat(0, 1), mat(1, 1));
		pitchRadian = -mat(2, 0) * M_PI / 2;
		yawRadian = 0;
	}
}

template <typename T>
inline
CMatrix<4, 4, T> make_rotation_matrix(const CMatrix<4, 4, T>& mat)
//CMatrix<4, 4, T> GetRotationMatrix(const CMatrix<4, 4, T>& mat)
{
	CMatrix<4, 4, T> mat2 = mat;
	mat2(0, 3) = 0;
	mat2(1, 3) = 0;
	mat2(2, 3) = 0;
	return mat2;
}

template <typename T>
inline
CMatrix<4, 4, T> make_rotation_matrix_from_roll_pitch_yaw(T rollRadian, T pitchRadian, T yawRadian)
//CMatrix<4, 4, T> GetRotationMatrixByRollPitchYaw(T rollRadian, T pitchRadian, T yawRadian)
{
	double cr = cos(rollRadian);
	double sr = sin(rollRadian);
	double cp = cos(pitchRadian);
	double sp = sin(pitchRadian);
	double cy = cos(yawRadian);
	double sy = sin(yawRadian);
	T element[] =
	{
		cr * cp, sr * cp, -sp, 0,
		cr *sp *sy - sr * cy, sr *sp *sy + cr * cy, cp * sy, 0,
		cr *sp *cy + sr * sy, sr *sp *cy - cr * sy, cp * cy, 0,
		0, 0, 0, 1,
	};
	CMatrix<4, 4, T> mat;
	return mat.Initialize(element);
}

template <typename T>
inline
//CMatrix<4, 4, T> GetRotationMatrixByEuler(T a, T b, T c)
CMatrix<4, 4, T> make_rotation_matrix_from_euler_angle(T a, T b, T c)
{
	double ca = cos(a);
	double sa = sin(a);
	double cb = cos(b);
	double sb = sin(b);
	double cc = cos(c);
	double sc = sin(c);
	T element[] =
	{
		ca *cb *cc - sa * sc, sa *cb *cc + ca * sc, -sb * cc, 0,
		-ca *cb *sc - sa * cc, -sa *cb *sc + ca * cc, sb * sc, 0,
		ca * sb, sa * sb, cb, 0,
		0, 0, 0, 1,
	};
	return CMatrix<4, 4, T>(element);
}

template <typename T>
inline
CMatrix<4, 4, T> make_rotation_matrix_from_axis(const CMatrix<3, 1, T>& axis, float cosangle, float sinangle)
//CMatrix<4, 4, T> GetRotationMatrixByAxis(const CMatrix<3, 1, T>& axis, float cosangle, float sinangle)
{
	T len = norm2_of(axis);
	if (!len)
	{
		return make_diagonal_matrix<T>(1, 1, 1, 1);
	}

	CVector<3, T> unitv = normalized_of(axis);
	T element[] =
	{
		unitv[0] *unitv[0] + (1 - unitv[0] * unitv[0]) *cosangle,
		unitv[0] *unitv[1] *(1 - cosangle) + unitv[2] *sinangle,
		unitv[0] *unitv[2] *(1 - cosangle) - unitv[1] *sinangle,
		0,
		unitv[0] *unitv[1] *(1 - cosangle) - unitv[2] *sinangle,
		unitv[1] *unitv[1] + (1 - unitv[1] * unitv[1]) *cosangle,
		unitv[1] *unitv[2] *(1 - cosangle) + unitv[0] *sinangle,
		0,
		unitv[0] *unitv[2] *(1 - cosangle) + unitv[1] *sinangle,
		unitv[1] *unitv[2] *(1 - cosangle) - unitv[0] *sinangle,
		unitv[2] *unitv[2] + (1 - unitv[2] * unitv[2]) *cosangle,
		0,
		0, 0, 0, 1,
	};
	return CMatrix<4, 4, T>(element);
}

template <typename T>
inline
CMatrix<4, 4, T> make_rotation_matrix_from_axis(const CMatrix<3, 1, T>& axis, float angle)
//CMatrix<4, 4, T> GetRotationMatrixByAxis(const CMatrix<3, 1, T>& axis, float angle)
{
	return make_rotation_matrix_from_axis(axis, cos(angle), sin(angle));
}

template <typename T>
inline
CMatrix<4, 4, T> make_rotation_matrix_from_mouse_drag(T rx, T ry, T dx, T dy)
//CMatrix<4, 4, T> GetRotationMatrixByMouseDrag(T rx, T ry, T dx, T dy)
{
	CVector<3, T> v1 = normalized_of(make_vector<T>(rx, ry, 1));
	CVector<3, T> v2 = normalized_of(make_vector<T>(rx + dx, ry + dy, 1));
	CVector<3, T> cr = cross(v1, v2);
	return make_rotation_matrix_from_axis(cr, dot(v1, v2), norm2_of(cr));
}

#if defined(MK_LBUTTON)
template <typename T>
inline
CMatrix<4, 4, T> make_rigid_matrix_from_mouse_drag(T sx, T sy, T dx, T dy, int nFlags)
//CMatrix<4, 4, T> GetTransformationMatrixByMouseDrag(T sx, T sy, T dx, T dy, int nFlags)
{
	// left button
	if ((nFlags & MK_LBUTTON) &&
	    !(nFlags & MK_MBUTTON) &&
	    !(nFlags & MK_RBUTTON))
	{
		return make_translation_matrix<T>(dx, -dy, 0);
	}
	// middle button
	else if ((!(nFlags & MK_LBUTTON) &&
	          (nFlags & MK_MBUTTON) &&
	          !(nFlags & MK_RBUTTON)) ||
	         ((nFlags & MK_LBUTTON) &&
	          !(nFlags & MK_MBUTTON) &&
	          (nFlags & MK_RBUTTON)))
	{
		return make_translation_matrix<T>(0, 0, dy);
	}
	// right button
	else if (!(nFlags & MK_LBUTTON) &&
	         !(nFlags & MK_MBUTTON) &&
	         (nFlags & MK_RBUTTON))
	{
		return make_rotation_matrix_from_mouse_drag(sx, -sy, dx, -dy);
	}
	throw std::logic_error("unknown mouse flag");
}
#endif

template <int nDimension, typename T>
inline
CMatrix<4, 4, T> make_rigid_matrix(const CMatrix<nDimension, nDimension, T>& rot, const CMatrix<nDimension, 1, T>& trans)
{
	STATIC_ASSERT(nDimension == 3 || nDimension == 4);
	T e[] =
	{
		rot(0, 0), rot(1, 0), rot(2, 0), 0,
		rot(0, 1), rot(1, 1), rot(2, 1), 0,
		rot(0, 2), rot(1, 2), rot(2, 2), 0,
		trans(0, 0), trans(1, 0), trans(2, 0), 1,
	};
	return CMatrix<4, 4, T>(e);
}

// similarity transformation

template <typename T>
inline
CMatrix<4, 4, T> make_scaling_matrix(T const& s)
//CMatrix<4, 4, T> GetScalingMatrix(T const& s)
{
	return make_scaling_matrix<T>(s, s, s);
}

template <typename T>
inline
CMatrix<4, 4, T> make_scaling_matrix(T const& sx, T const& sy, T const& sz)
//CMatrix<4, 4, T> GetScalingMatrix(T const& sx, T const& sy, T const& sz)
{
	CMatrix<4, 4, T> mat = make_diagonal_matrix<T>(1, 1, 1, 1);
	mat(0, 0) = sx;
	mat(1, 1) = sy;
	mat(2, 2) = sz;
	return mat;
}

// mat4x4 x vec3 --or-- mat3x4 x vec3
template <int nRows, typename T>
inline
const CVector<3, T> AffineTransform(const CMatrix<nRows, 4, T>& mat, const CVector<3, T>& vec)
{
	STATIC_ASSERT(nRows == 3 || nRows == 4);
	return make_vector<T>(
	           mat(0, 0) * vec[0] + mat(0, 1) * vec[1] + mat(0, 2) * vec[2] + mat(0, 3),
	           mat(1, 0) * vec[0] + mat(1, 1) * vec[1] + mat(1, 2) * vec[2] + mat(1, 3),
	           mat(2, 0) * vec[0] + mat(2, 1) * vec[1] + mat(2, 2) * vec[2] + mat(2, 3));
}

//----------------------------------------------------------------------
// dynamic matrix
//----------------------------------------------------------------------

template <typename T>
class CDynamicMatrix
{
public:
	// constructors
	CDynamicMatrix(void) :
		m_nNumRows(0),
		m_nNumCols(0),
		m(0)
	{
	}

	CDynamicMatrix(int nRows, int nCols) :
		m_nNumRows(0),
		m_nNumCols(0),
		m(0)
	{
		Resize(nRows, nCols);
	}

	CDynamicMatrix(int nRows, int nCols, T const *const element) :
		m_nNumRows(0),
		m_nNumCols(0),
		m(0)
	{
		Resize(nRows, nCols);
		Initialize(element);
	}

	// copy
	CDynamicMatrix(const CDynamicMatrix& mat) :
		m_nNumRows(0),
		m_nNumCols(0),
		m(0)
	{
		*this = mat;
	}

	const CDynamicMatrix& operator =(const CDynamicMatrix& mat)
	{
		Resize(mat.GetNumRows(), mat.GetNumCols());
		return Initialize(mat.m);
	}

	// coercing
	template <typename T2>
	CDynamicMatrix(const CDynamicMatrix<T2>& mat) :
		m_nNumRows(0),
		m_nNumCols(0),
		m(0)
	{
		*this = mat;
	}

	template <typename T2>
	CDynamicMatrix(int nRows, int nCols, T2 const *const element) :
		m_nNumRows(0),
		m_nNumCols(0),
		m(0)
	{
		Resize(nRows, nCols);
		Initialize(element);
	}

	template <typename T2>
	const CDynamicMatrix& operator =(const CDynamicMatrix<T2>& mat)
	{
		Resize(mat.GetNumRows(), mat.GetNumCols());
		for (int i = 0; i < mat.GetNumRows()*mat.GetNumCols(); i++)
		{
			m[i] = mat.ptr()[i];
		}
		return *this;
	}

	// conversion from static version
	template <int nRows, int nCols, typename T2>
	CDynamicMatrix(const CMatrix<nRows, nCols, T2>& mat) :
		m_nNumRows(0),
		m_nNumCols(0),
		m(0)
	{
		*this = mat;
	}

	template <int nRows, int nCols, typename T2>
	const CDynamicMatrix& operator =(const CMatrix<nRows, nCols, T2>& mat)
	{
		Resize(nRows, nCols);
		for (int i = 0; i < nRows * nCols; i++)
		{
			m[i] = mat.ptr()[i];
		}
		return *this;
	}

	// destructor
	~CDynamicMatrix(void)
	{
		if (m)
		{
			delete [] m;
		}
	}

	// helpers
	void Resize(int nRows, int nCols)
	{
		if (m_nNumRows == nRows && m_nNumCols == nCols)
		{
			return;
		}

		if (m)
		{
			delete [] m;
		}
		m = new T [nRows * nCols];
		if (!m_nNumRows && !m_nNumCols)
		{
			for (int i = 0; i < nRows * nCols; i++)
			{
				m[i] = T();
			}
		}
		m_nNumRows = nRows;
		m_nNumCols = nCols;
	}

	template <typename T2>
	const CDynamicMatrix& Initialize(T2 const *const e)
	{
		for (int i = 0; i < m_nNumRows * m_nNumCols; i++)
		{
			m[i] = e[i];
		}
		return *this;
	}

	template <typename T2>
	void Fill(T2 const e)
	{
		for (int i = 0; i < m_nNumRows * m_nNumCols; i++)
		{
			m[i] = e;
		}
	}

	// accessors
	T& operator()(const int nRow, const int nCol)
	{
		_ASSERTE(nRow < m_nNumRows && nCol < m_nNumCols);
		return m[nRow + m_nNumRows * nCol];
	}

	T const& operator()(const int nRow, const int nCol) const
	{
		_ASSERTE(nRow < m_nNumRows && nCol < m_nNumCols);
		return m[nRow + m_nNumRows * nCol];
	}

	// arithmetic operators
	const CDynamicMatrix& operator +=(const CDynamicMatrix& mat)
	{
		_ASSERTE(m_nNumRows == mat.GetNumRows() && m_nNumCols == mat.GetNumCols());
		for (int i = 0; i < m_nNumRows * m_nNumCols; i++)
		{
			m[i] += mat.m[i];
		}
		return *this;
	}

	const CDynamicMatrix& operator -=(const CDynamicMatrix& mat)
	{
		_ASSERTE(m_nNumRows == mat.GetNumRows() && m_nNumCols == mat.GetNumCols());
		for (int i = 0; i < m_nNumRows * m_nNumCols; i++)
		{
			m[i] -= mat.m[i];
		}
		return *this;
	}

	const CDynamicMatrix& operator *=(const CDynamicMatrix& mat)
	{
		return *this = *this * mat;
	}

	const CDynamicMatrix& operator *=(double s)
	{
		for (int i = 0; i < m_nNumRows * m_nNumCols; i++)
		{
			m[i] *= s;
		}
		return *this;
	}

	const CDynamicMatrix& operator /=(double s)
	{
		_ASSERTE(s != 0);
		for (int i = 0; i < m_nNumRows * m_nNumCols; i++)
		{
			m[i] /= s;
		}
		return *this;
	}

	// unary operator
	CDynamicMatrix operator -(void) const
	{
		return *this * (-1);
	}
	/*
		const CDynamicMatrix& operator +(void) const
		{
			return *this;
		}
		*/
	// binary operator
	CDynamicMatrix operator +(const CDynamicMatrix& v2) const
	{
		return CDynamicMatrix(*this) += v2;
	}

	CDynamicMatrix operator -(const CDynamicMatrix& v2) const
	{
		return CDynamicMatrix(*this) -= v2;
	}

	CDynamicMatrix operator *(const CDynamicMatrix& mat) const
	{
		_ASSERTE(m_nNumCols == mat.GetNumRows());
		CDynamicMatrix iret(m_nNumRows, mat.GetNumCols());
		for (int nc = 0; nc < mat.GetNumCols(); nc++)
		{
			for (int nr = 0; nr < m_nNumRows; nr++)
			{
				iret(nr, nc) = 0;
				for (int ncc = 0; ncc < mat.GetNumRows(); ncc++)
				{
					iret(nr, nc) += (*this)(nr, ncc) * mat(ncc, nc);
				}
			}
		}
		return iret;
	}

	CDynamicMatrix operator *(double s) const
	{
		return CDynamicMatrix(*this) *= s;
	}

	CDynamicMatrix operator /(double s) const
	{
		return CDynamicMatrix(*this) /= s;
	}

	friend CDynamicMatrix operator *(double s, const CDynamicMatrix& mat)
	{
		return mat * s;
	}

	// relationship operators
	bool operator ==(const CDynamicMatrix& mat2) const
	{
		return
		    m_nNumRows == mat2.GetNumRows() &&
		    m_nNumCols == mat2.GetNumCols() &&
		    !memcmp(m, mat2.m, sizeof(T) * m_nNumRows * m_nNumCols);
	}

	bool operator !=(const CDynamicMatrix& mat2) const
	{
		return !(*this == mat2);
	}

	// utilities

	T min(void) const
	{
		T iret = m[0];
		for (int i = 1; i < m_nNumRows * m_nNumCols; i++)
		{
			iret = std::min(iret, m[i]);
		}
		return iret;
	}

	T max(void) const
	{
		T iret = m[0];
		for (int i = 1; i < m_nNumRows * m_nNumCols; i++)
		{
			iret = std::max(iret, m[i]);
		}
		return iret;
	}

	CDynamicMatrix abs(void) const
	{
		CDynamicMatrix result(m_nNumRows, m_nNumCols);
		for (int i = 0; i < m_nNumRows * m_nNumCols; i++)
		{
			result.m[i] = std::abs(m[i]);
		}
		return result;
	}

	T const *ptr(void) const
	{
		return m;
	}

	T *ptr(void)
	{
		return m;
	}

	int GetNumRows(void) const
	{
		return m_nNumRows;
	}

	int GetNumCols(void) const
	{
		return m_nNumCols;
	}

	CDynamicMatrix<T> GetColumn(int c) const
	{
		throw std::logic_error("not implemented");
	}

	CDynamicMatrix<T> GetRow(int r) const
	{
		throw std::logic_error("not implemented");
	}

	void GetSubMatrix(int nr, int nc, int sr, int sc, CDynamicMatrix& mat) const
	{
		_ASSERTE(nr <= GetNumRows());
		_ASSERTE(nc <= GetNumCols());
		mat.Resize(nr, nc);
		for (int c = 0; c < nc; c++)
		{
			for (int r = 0; r < nr; r++)
			{
				mat(r, c) = (*this)(sr + r, sc + c);
			}
		}
	}

	const CDynamicMatrix& ExpandRows(int nrows)
	{
		int nNumNewRows = m_nNumRows + nrows;
		T *e = new T [nNumNewRows * m_nNumCols];
		for (int c = 0; c < m_nNumCols; c++)
		{
			for (int r = 0; r < m_nNumRows; r++)
			{
				e[r + c * nNumNewRows] = m[r + c * m_nNumRows];
			}
			for (int r = m_nNumRows; r < m_nNumRows + nrows; r++)
			{
				e[r + c * nNumNewRows] = 0;
			}
		}
		delete [] m;
		m = e;
		m_nNumRows = nNumNewRows;
		return *this;
	}

	const CDynamicMatrix& ExpandCols(int ncols)
	{
		int nNumNewCols = m_nNumCols + ncols;
		T *e = new T [m_nNumRows * nNumNewCols];
		for (int c = 0; c < m_nNumCols; c++)
		{
			for (int r = 0; r < m_nNumRows; r++)
			{
				e[r + c * m_nNumRows] = m[r + c * m_nNumRows];
			}
		}
		for (int c = m_nNumCols; c < nNumNewCols; c++)
		{
			for (int r = 0; r < m_nNumRows; r++)
			{
				e[r + c * m_nNumRows] = 0;
			}
		}
		delete [] m;
		m = e;
		m_nNumCols = nNumNewCols;
		return *this;
	}
	/*
	CDynamicMatrix AppendCols(const CDynamicMatrix& mat) const
	{
		_ASSERTE(m_nNumRows == mat.GetNumRows());
		CDynamicMatrix result(m_nNumRows, m_nNumCols + mat.GetNumCols());
		// left half
		int nelem1 = m_nNumRows * m_nNumCols;
		memcpy(result.m, m, sizeof(T) * nelem1);
		// right half
		int nelem2 = mat.GetNumRows() * mat.GetNumCols();
		memcpy(result.m + nelem1, mat.m, sizeof(T) * nelem2);
		return result;
	}

	CDynamicMatrix AppendRows(const CDynamicMatrix& mat) const
	{
		_ASSERTE(m_nNumCols == mat.GetNumCols());
		CDynamicMatrix result(m_nNumRows + mat.GetNumRows(), m_nNumCols);
		for (int c = 0; c < m_nNumCols; c++)
		{
			for (int r = 0; r < m_nNumRows; r++)
			{
				result(r, c) = (*this)(r, c);
			}
			for (int r = 0; r < mat.GetNumRows(); r++)
			{
				result(r + m_nNumRows, c) = mat(r, c);
			}
		}
		return result;
	}
	*/
	void Read(const std::string& filename)
	{
		TRACE("matrix <= %s\n", filename.c_str());
		FILE *fr = fopen(filename.c_str(), "rb");
		if (!fr)
		{
			throw std::runtime_error("failed to open file in " __FUNCTION__);
		}
		int nrows, ncols;
		fscanf(fr, "%d%d", &nrows, &ncols);
		Resize(nrows, ncols);
		for (int r = 0; r < m_nNumRows; r++)
		{
			for (int c = 0; c < m_nNumCols; c++)
			{
				double el;
				if (fscanf(fr, "%lf", &el) != 1)
				{
					fclose(fr);
					throw std::runtime_error("failed to parse file in " __FUNCTION__);
				}
				(*this)(r, c) = el;
			}
		}

		fclose(fr);
	}

	void Write(const std::string& filename) const
	{
		TRACE("matrix => %s\n", filename.c_str());
		FILE *fw = fopen(filename.c_str(), "wb");
		if (!fw)
		{
			throw std::runtime_error("failed to open file in " __FUNCTION__);
		}
		fprintf(fw, "%d %d\n", m_nNumRows, m_nNumCols);
		for (int r = 0; r < m_nNumRows; r++)
		{
			for (int c = 0; c < m_nNumCols; c++)
			{
				fprintf(fw, "% le ", double((*this)(r, c)));
			}
			fprintf(fw, "\n");
		}

		fclose(fw);
	}

protected:
	int m_nNumRows;
	int m_nNumCols;
	T *m; // row major order
};

//----------------------------------------------------------------------
// dynamic vector
//----------------------------------------------------------------------

template <typename T>
class CDynamicVector : public CDynamicMatrix<T>
{
public:
	// constructors
	CDynamicVector(void)
	{
	}

	explicit CDynamicVector(int nDim)
		: CDynamicMatrix<T>(nDim, 1)
	{
	}

	CDynamicVector(int nDim, T const *const element)
		: CDynamicMatrix<T>(nDim, 1, element)
	{
	}

	// copy
	CDynamicVector(const CDynamicVector& vec)
		: CDynamicMatrix(vec)
	{
	}

	const CDynamicVector& operator =(const CDynamicVector& vec)
	{
		return CDynamicMatrix::operator =(vec);
	}

	// coercion
	template <typename T2>
	CDynamicVector(const CDynamicMatrix<T2>& vec)
		: CDynamicMatrix(vec)
	{
	}

	template <typename T2>
	const CDynamicVector& operator =(const CDynamicMatrix<T2>& vec)
	{
		return CDynamicMatrix::operator =(vec);
	}

	// conversion from CMatrix
	template <int nRows, typename T2>
	CDynamicVector(const CMatrix<nRows, 1, T2>& vec)
		: CDynamicMatrix(vec)
	{
	}

	template <int nRows, typename T2>
	const CDynamicVector& operator =(const CMatrix<nRows, 1, T2>& vec)
	{
		return CDynamicMatrix::operator =(vec);
	}

	// helpers
	void Resize(int nRows)
	{
		CDynamicMatrix::Resize(nRows, 1);
	}

	// operators
	T& operator [](const int n)
	{
		_ASSERTE(0 <= n && n < m_nNumRows);
		return m[n];
	}

	T const& operator [](const int n) const
	{
		_ASSERTE(0 <= n && n < m_nNumRows);
		return m[n];
	}
};

// helpers

template <typename T>
inline
CDynamicMatrix<T> normalized_of(const CDynamicMatrix<T>& vec)
{
	T len = norm2_of(vec);
	if (len == 0)
	{
		TRACE("warning: attempted to normalize a zero vector.\n");
		return vec;
	}
	else
	{
		return vec / len;
	}
}

template <typename T>
inline
CDynamicMatrix<T> transpose_of(const CDynamicMatrix<T>& mat)
{
	CDynamicMatrix<T> iret(mat.GetNumCols(), mat.GetNumRows());
	for (int c = 0; c < mat.GetNumCols(); c++)
	{
		for (int r = 0; r < mat.GetNumRows(); r++)
		{
			iret(c, r) = mat(r, c);
		}
	}
	return iret;
}

// helpers for CDynamicVector

template <typename T>
inline
double norm2_of(const CDynamicMatrix<T>& vec)
{
	return sqrt(static_cast<double>(dot(vec, vec)));
}

template <typename T>
inline
T norm2_squared_of(const CDynamicMatrix<T>& vec)
{
	return dot(vec, vec);
}

template <typename T>
inline
T dot(const CDynamicMatrix<T>& v1, const CDynamicMatrix<T>& v2)
{
	_ASSERTE(v1.GetNumCols() == 1 && v2.GetNumCols() == 1);
	_ASSERTE(v1.GetNumRows() == v2.GetNumRows());
	T iret = 0;
	for (int d = 0; d < v1.GetNumRows(); d++)
	{
		iret += v1(d, 0) * v2(d, 0);
	}
	return iret;
}

template <typename T>
inline
CDynamicMatrix<T> skew_symmetric_of(CDynamicMatrix<T>& vec)
{
	if (vec.GetNumRows() != 3)
	{
		throw std::runtime_error("not a 3-vector in " __FUNCTION__);
	}

	CDynamicMatrix<T> iret(3, 3);
	T element[] =
	{
		0, vec(2, 0), -vec(1, 0),
		-vec(2, 0), 0, vec(0, 0),
		vec(1, 0), -vec(0, 0), 0
	};
	return iret.Initialize(element);
}

} // namespace slib