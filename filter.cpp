#include "stdafx.h"
#include "depth2ply.h"

//12/11/13Å@Ç©ÇÌÇËÇ…macÇÃAccelerateÇégÇ¢Ç‹Ç∑ÇÈ
#include <Accelerate/Accelerate.h>

using namespace slib;

void bilateral_filter(unsigned short *depth, int window, float depth_ratio)
{
	static unsigned short work[640 * 480 * 2];
	memcpy(work, depth, 640 * 480 * 2);
	float var1 = depth_ratio * depth_ratio;
	for (int y = 0; y < 480; y++)
	{
		printf("\rbilateral filtering %d%%", 100 * y / 479);
		//#pragma omp parallel for
		for (int x = 0; x < 640; x++)
		{
			int index = (work[x + 640 * y] & 7);
			short ref = (work[x + 640 * y] >> 3);
			if (!ref)
			{
				continue;
			}
			float nvar1 = var1 * ref * ref;
			float weight = 0;
			float sum = 0;
			for (int iy = -window; iy <= window; iy++)
			{
				for (int ix = -window; ix <= window; ix++)
				{
					int px = x + ix;
					int py = y + iy;
					if (px >= 0 && px < 640 && py >= 0 && py < 480)
					{
						short d = (work[px + 640 * py] >> 3);
						if (!d)
						{
							continue;
						}
						float w1 = exp(-(ref - d) * (ref - d) / nvar1);
						float w2 = 1; // exp(-(ix*ix+iy*iy)/var2);
						sum += w1 * w2 * d;
						weight += w1 * w2;
					}
				}
			}
			if (weight != 0)
			{
				sum /= weight;
			}
			depth[x+640*y]=((short)sum << 3 | index);
		}
	}
	printf("\n");
}

void fast_least_squares(float *a, int ndata)
{
	const char jobu = 'A';
	const char jobvt = 'N';
	const MKL_INT m = 4;
	const MKL_INT n = ndata;
	const MKL_INT minmn = std::min(m, n);
	const MKL_INT lda = m;
	float s[4];
	float u[4 * 4];
	const MKL_INT ldu = m;
	const MKL_INT ldvt = n;
	float work[1024];
	MKL_INT lwork = 1024;
	MKL_INT info;
	sgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, 0, &ldvt, work, &lwork, &info);
	//dgesv_(&m, &n, (double)a, &lda, s, u, &ldu, 0, &ldvt, work, &lwork, &info);
	for (int i = 0; i < 4; i++)
	{
		a[i] = u[12 + i];
	}
}

void least_squares_fitting(unsigned short *depth, int window, float depth_ratio)
{
	static vertex vtx[640 * 480];
	#pragma omp parallel for
	for (int y = 0; y < 480; y++)
	{
		for (int x = 0; x < 640; x++)
		{
			unsigned short w = depth[x + 640 * y];
			if (w)
			{
				depth2pos(x, y, w, vtx[x + 640 * y]);
			}
		}
	}

	static unsigned short work[640 * 480 * 2];
	memcpy(work, depth, 640 * 480 * 2);
	//int window = m_option.window;
	float var1 = depth_ratio * depth_ratio;
	for (int y = 0; y < 480; y++)
	{
		printf("\rleast squares filtering %d%%", 100 * y / 479);
		#pragma omp parallel for
		for (int x = 0; x < 640; x++)
		{
			int index = (work[x + 640 * y] & 7);
			short ref = (work[x + 640 * y] >> 3);
			if (!ref)
			{
				continue;
			}
			float nvar1 = var1 * ref * ref;
			float mat[1024];
			int ndata = 0;
			for (int iy = -window; iy <= window; iy++)
			{
				for (int ix = -window; ix <= window; ix++)
				{
					int px = x + ix;
					int py = y + iy;
					if (px >= 0 && px < 640 && py >= 0 && py < 480)
					{
						short d = (work[px + 640 * py] >> 3);
						if (!d)
						{
							continue;
						}
						float w1 = exp(-(ref - d) * (ref - d) / nvar1);
						float w2 = 1; // exp(-(ix*ix+iy*iy)/var2);
						mat[ndata * 4 + 0] = w1 * w2 * vtx[px + 640 * py].x;
						mat[ndata * 4 + 1] = w1 * w2 * vtx[px + 640 * py].y;
						mat[ndata * 4 + 2] = w1 * w2 * vtx[px + 640 * py].z;
						mat[ndata * 4 + 3] = w1 * w2;
						ndata++;
					}
				}
			}
			if (ndata > 2)
			{
				fast_least_squares(mat, ndata);
				vertex& v = vtx[x + 640 * y];
				float t = -mat[3] / dot(make_vector(v.x, v.y, v.z), make_vector(mat[0], mat[1], mat[2]));
				ref *= t;
			}
			depth[x + 640 * y] = ((ref << 3) | index);
		}
	}
	printf("\n");
}

