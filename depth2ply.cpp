// depth2ply.cpp : Defines the entry point for the console application.
//


#include "stdafx.h"
#include "depth2ply.h"

using namespace slib;

//オプション構造体
struct
{
	std::string filename;
	bool depth ;
	bool mesh ;
	bool timestamp ;
	bool bilateral;
	bool least_squares;
	bool mask;
	bool pipe;
	int frame_begin;
	int frame_end; // -2: until the last frame
	float depth_ratio;
	int window;
	float threshold;
	int index;
}
m_option =
{
	"",
	false, false, false, false, false, false, false,
	0, -1,
	0.03, 5,
	50,
	0,
};
//使い方描画関数
void print_usage(char *argv1)
{
	printf("Usage: %s [-b [depth_ratio window]] [-d] [-f frames] [-i index] [-l [depth_ratio window]] [-m] [-p [threshold]] [-r] [-t] <depth>\n", argv1);
	printf("frames: n, n-m, n-, -m\n");
	exit(-1);
}

void parse_frames(char *arg)
{
	if (strchr(arg, '-'))
	{
		char *p = strchr(arg, '-');
		*p = 0;
		m_option.frame_begin = atoi(arg);
		m_option.frame_end = atoi(p + 1);
	}
	else
	{
		m_option.frame_begin = atoi(arg);
		m_option.frame_end = m_option.frame_begin + 1;
	}
	if (!m_option.frame_end)
	{
		m_option.frame_end = -1;
	}
}

void parse_option(int argc, char *argv[])
{
	for (int i = 1; i < argc; i++)
	{
		if (argv[i][0] != '-')
		{
			m_option.filename = argv[i];
			return;
		}
		else
		{
			switch (argv[i][1])
			{
			case 'b':
				m_option.bilateral = true;
				m_option.least_squares = false;
				if (i + 3 < argc && argv[i + 1][0] != '-' && argv[i + 2][0] != '-')
				{
					m_option.depth_ratio = atof(argv[i + 1]);
					m_option.window = atoi(argv[i + 2]);
					i += 2;
				}
				break;
			case 'd':
				m_option.depth = true;
				break;
			case 'f':
				if (i + 2 < argc)
				{
					parse_frames(argv[i + 1]);
					i++;
				}
				else
				{
					print_usage(argv[0]);
				}
				break;
			case 'i':
				if (i + 2 < argc)
				{
					m_option.index = atoi(argv[i + 1]);
					i++;
				}
				else
				{
					print_usage(argv[0]);
				}
				break;
			case 'l':
				m_option.bilateral = false;
				m_option.least_squares = true;
				if (i + 3 < argc && argv[i + 1][0] != '-' && argv[i + 2][0] != '-')
				{
					m_option.depth_ratio = atof(argv[i + 1]);
					m_option.window = atoi(argv[i + 2]);
					i += 2;
				}
				break;
			case 'm':
				m_option.mask = true;
				break;
			case 'p':
				m_option.mesh = true;
				if (i + 2 < argc && argv[i + 1][0] != '-')
				{
					m_option.threshold = atof(argv[i + 1]);
					i++;
				}
				break;
			case 'r':
				m_option.pipe = true;
				break;
			case 't':
				m_option.timestamp = true;
				break;
			default:
				print_usage(argv[0]);
				break;
			}
		}
	}
	if (m_option.filename.empty())
	{
		print_usage(argv[0]);
	}
//	if (!m_option.bilateral && !m_option.least_squares)
//	{
//		m_option.pipe=false;
//	}
}

int main(int argc, char *argv[])
{
	//例外処理を読み込み
	try
	{
		//オプションの読み込み
		parse_option(argc, argv);

		//オプションからファイルを開く
		FILE *fr = fopen(m_option.filename.c_str(), "rb");
		if (!fr)
		{
			print_usage(argv[1]);
		}

		//フレームを読み込み
		int nframes, d;
		fread(&nframes, sizeof(int), 1, fr);
		printf("#frames = %d\n", nframes);
		if (m_option.frame_end < 0)
		{
			m_option.frame_end = nframes;
		}

		// check file format
		fseek(fr, 0, SEEK_END);
		long filesize = ftell(fr);
		if (filesize != (640 * 480 * 2 + 4) * nframes + 4)
		{
			fclose(fr);
			printf("error: wrong file size.");
			if (filesize == (640 * 480 * 2 + 4) * nframes + 12)
			{
				printf("perhaps in old format.");
			}
			printf("\n");
			print_usage(argv[1]);
		}
		fseek(fr, 4, SEEK_SET);

		FILE *fw = 0;
		if (m_option.pipe && m_option.frame_begin == 0 && m_option.frame_end == nframes)
		{
			if (m_option.bilateral)
			{
				std::string file = string_replace(m_option.filename, ".depth", "-bilateral.depth");
				fw = fopen(file.c_str(), "wb");
			}
			else if (m_option.least_squares)
			{
				std::string file = string_replace(m_option.filename, ".depth", "-leastsquares.depth");
				fw = fopen(file.c_str(), "wb");
			}
		}
		if (fw)
		{
			fwrite(&nframes, 4, 1, fw);
		}

		unsigned long *tm = new unsigned long [nframes];
		for (int frame = 0; frame < nframes; frame++)
		{
			fread(&tm[frame], sizeof(unsigned long), 1, fr);
			static unsigned short depth[640 * 480];
			fread(depth, sizeof(unsigned short), 640 * 480, fr);

			if (frame < m_option.frame_begin || frame >= m_option.frame_end)
			{
				continue;
			}

			// filter
			if (m_option.bilateral)
			{
				bilateral_filter(depth, m_option.window, m_option.depth_ratio);
			}
			else if (m_option.least_squares)
			{
				least_squares_fitting(depth, m_option.window, m_option.depth_ratio);
			}

			if (fw)
			{
				printf("depth <= frame #%d\n", frame);
				fwrite(&tm[frame], 4, 1, fw);
				fwrite(depth, 2, 640 * 480, fw);
			}
			// depth
			if (m_option.depth)
			{
				export_depth(depth, m_option.index, slib::string_replace(m_option.filename, ".depth", slib::string_format("-depth-%03d.bmp", frame)));
			}

			if (m_option.mask)
			{
				export_mask(depth, m_option.index, slib::string_replace(m_option.filename, ".depth", slib::string_format("-mask-%03d.bmp", frame)));
			}

			// mesh
			if (m_option.mesh)
			{
				std::string meshname = string_replace(m_option.filename, ".depth", string_format("-%03d.ply", frame));;
				if (m_option.bilateral)
				{
					meshname = string_replace(meshname, ".ply", string_format("-bl-%.2f-%d.ply", m_option.depth_ratio, m_option.window));
				}
				if (m_option.least_squares)
				{
					meshname = string_replace(meshname, ".ply", string_format("-ls-%.2f-%d.ply", m_option.depth_ratio, m_option.window));
				}
				export_mesh(depth, m_option.index,  meshname, m_option.threshold);
			}
		}

		// timestamp
		if (m_option.timestamp)
		{
			export_timestamp(tm, nframes, string_replace(m_option.filename, ".depth", "-time.txt"));
		}

		if (fw)
		{
			fclose(fw);
		}
		fclose(fr);
	}
	catch (const std::exception& e)
	{
		printf("exception: %s\n", e.what());
	}
	return 0;
}

