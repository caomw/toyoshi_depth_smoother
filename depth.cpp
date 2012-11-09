#include "stdafx.h"
#include "depth2ply.h"

using namespace slib;

void export_mask(unsigned short *depth, int index, const std::string& f)
{
	static unsigned char mask[640 * 480];
	for (int i = 0; i < 640 * 480; i++)
	{
		if ((depth[i] & 0x007) == index)
		{
			mask[i] = 255;
		}
		else
		{
			mask[i] = 0;
		}
	}
	slib::image::FlipVertically(640, 480, mask, mask);
	slib::image::WriteBmp(f.c_str(), 640, 480, 1, mask);
	printf("bmp => %s\n", f.c_str());
}

void export_depth(unsigned short *depth, int index, const std::string& f)
{
	static unsigned char depth_byte[640 * 480];
	for (int i = 0; i < 640 * 480; i++)
	{
		if ((depth[i] & 7) == index)
		{
			depth_byte[i] = (depth[i] >> 8);
		}
		else
		{
			depth_byte[i] = 0;
		}
	}
	slib::image::FlipVertically(640, 480, depth_byte, depth_byte);
	slib::image::WriteBmp(f.c_str(), 640, 480, 1, depth_byte);
	printf("bmp => %s\n", f.c_str());
}

void export_timestamp(unsigned long *tm, int nframes, const std::string& f)
{
	FILE *fw = fopen(f.c_str(), "wb");
	if (!fw)
	{
		throw std::runtime_error("failed to open");
	}
	for (int i = 0; i < nframes; i++)
	{
		fprintf(fw, "%d\n", tm[i]);
	}
	fclose(fw);
	printf("txt => %s\n", f.c_str());
}