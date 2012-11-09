#pragma once

struct vertex
{
	float x, y, z;
};

#pragma pack(push, 1)
struct triangle
{
	unsigned char d;
	int a, b, c;
};
#pragma pack(pop)

void bilateral_filter(unsigned short *depth,  int window, float depth_ratio);
void least_squares_fitting(unsigned short *depth, int window, float depth_ratio);
void depth2pos(int x, int y, unsigned short w, struct vertex& pos);
bool isconnected(int i1, int i2, int i3, const int *vid, const vertex *v, float threshold);
void largest_component(const unsigned short *depth, slib::image::CBmpImage& bmp, float depth_ratio);

void export_mask(unsigned short *depth, int index, const std::string& f);
void export_depth(unsigned short *depth, int index, const std::string& f);
void export_mesh(unsigned short *depth, int index, const std::string& meshname, float threshold);
void export_timestamp(unsigned long *tm, int nframes, const std::string& f);
