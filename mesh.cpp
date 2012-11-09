#include "stdafx.h"
#include "depth2ply.h"

using namespace slib;

static const float finv = 3.501e-3f / 2;
static const float cx = 320;
static const float cy = 240;

void depth2pos(int x, int y, unsigned short w, vertex& pos)
{
	float Z = (w >> 3);
	float X = (x - cx) * finv * Z;
	float Y = -(y - cy) * finv * Z;
	pos.x = X ;
	pos.y = Y ;
	pos.z = Z ;
}

bool isconnected(int i1, int i2, int i3, const int *vid, const vertex *v, float threshold)
{
	return
	    vid[i1] != -1 && vid[i2] != -1 && vid[i3] != -1 &&
	    abs(v[i1].z - v[i2].z) < threshold &&
	    abs(v[i2].z - v[i3].z) < threshold &&
	    abs(v[i3].z - v[i1].z) < threshold ;
}

void export_mesh(unsigned short *depth, int index, const std::string& meshname, float threshold)
{
	int nverts = 0;
	static vertex v [640 * 480];
	static int vid [640 * 480];
	for (int y = 0; y < 480; y++)
	{
		for (int x = 0; x < 640; x++)
		{
			int idx = x + y * 640;
			unsigned short w = depth[idx];
			int id = (w & 0x007);
			if (!w || id != index)
			{
				vid[idx] = -1;
			}
			else
			{
				depth2pos(x, y, w, v[idx]);
				vid[idx] = nverts++;
			}
		}
	}
	if (!nverts)
	{
		return;
	}

	int ntriangles = 0;
	static triangle t [640 * 480 * 2];
	for (int y = 0; y < 480 - 1; y++)
	{
		for (int x = 0; x < 640 - 1; x++)
		{
			int idx = x + y * 640;
			if (isconnected(idx, idx + 1, idx + 640, vid, v, threshold))
			{
				t[ntriangles].a = vid[idx ];
				t[ntriangles].b = vid[idx + 1];
				t[ntriangles].c = vid[idx + 640];
				t[ntriangles].d = 3;
				ntriangles++;
			}
			if (isconnected(idx + 1, idx + 640, idx + 640 + 1, vid, v, threshold))
			{
				t[ntriangles].a = vid[idx + 640];
				t[ntriangles].b = vid[idx + 1];
				t[ntriangles].c = vid[idx + 640 + 1];
				t[ntriangles].d = 3;
				ntriangles++;
			}
		}
	}
	/*
	std::string f = string_replace(m_option.filename, ".depth", string_format("-%03d.ply", frame));
	if (m_option.bilateral)
	{
		f = string_replace(m_option.filename, ".depth", string_format("-%03d-bl-%.2f-%d.ply", frame, m_option.depth_ratio, m_option.window));
	}
	if (m_option.least_squares)
	{
		f = string_replace(m_option.filename, ".depth", string_format("-%03d-ls-%.2f-%d.ply", frame, m_option.depth_ratio, m_option.window));
	}
	*/
	printf("ply => %s\n", meshname.c_str());
	FILE *fw = fopen(meshname.c_str(), "wb");
	if (!fw)
	{
		throw std::runtime_error("failed to open");
	}
	fprintf(fw,
	        "ply\n"
	        "format binary_little_endian 1.0\n"
	        "element vertex %d\n"
	        "property float x\n"
	        "property float y\n"
	        "property float z\n"
	        "element face %d\n"
	        "property list uchar int vertex_indices\n"
	        "end_header\n", nverts, ntriangles);
	for (int idx = 0; idx < 640 * 480; idx++)
	{
		if (vid[idx] != -1)
		{
			fwrite(&v[idx], sizeof(vertex), 1, fw);
		}
	}
	fwrite(t, sizeof(triangle), ntriangles, fw);
	fclose(fw);
}