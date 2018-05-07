#ifndef METADATA_H
#define METADATA_H
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
using namespace std;
/* Viewport-related params */
struct Viewport{
	double Fh;				// Horizontal FoV (rad)
	double Fv;				// Vertical FoV (rad)
	int vp_W;				// Viewport's width (Number of pixels)
	int vp_H;				// Viewport's height (Number of pixels)
	int erp_W;				// ERP's width
	int erp_H;				// ERP's height
	int tile_W;				// tile's width
	int tile_H;				// tile's height
	int No_face;			// Number of faces
	int face_W;				// Face's width
	int face_H; 			// Face's height
	int No_tile_h;			// Number of vertical tiles per faces
	int No_tile_v;			// Number of horizontal tiles per faces
	int No_tile;			// Number of tiles
	string vmask_data_loc;	// Folder containing vmask data
	int*** vmask;			// Visible mask
	int*** pixel;			// Number of visible pixels
};
	/* Video-related params*/
struct Video{
	string name;
	string data_loc;		// folder containing video data
	double SD; 				// segment duration in secs
	double SESS_DUR; 		// session duration in secs
	int FPS;				// Frame rate
	int NO_SEG_FULL;		// Number of total segments
	int NO_SEG;				// Number of origin segments
	int BUFF;				// buffer size as number of frames
	int INTERVAL;			// Number of frames per interval
	int NO_FRAME;			// Total number of frames
	int NO_FRAME_ORIGIN;	// Total number of origin frames
	int NO_VER;				// Number of quality versions
	double*** TILE_BR;		// Tiles' bitrates
	double*** TILE_SIZE;	// Tiles' sizes
	double*** TILE_MSE;		// Tiles' MSEs
	double*** TILE_PSNR;	// Tiles' PSNRs
	double** DASH_BR;		// DASH's version bitrates (ERP)
	double** DASH_MSE;		// DASH's version MSEs (ERP)
	double** DASH_PSNR;		// DASH's version PSNRs (ERP)
	double** CMP_DASH_BR;		// DASH's version bitrates (CMP)
	double** CMP_DASH_MSE;		// DASH's version MSEs (CMP)
	double** CMP_DASH_PSNR;		// DASH's version PSNRs (CMP)
	double* BR;				// Tilings' bitrates (ERP)
	double* PSNR;			// Tilings' PSNRs (ERP)
};
	/* BellLab */
struct BellLab{
	double speed; 			// head moving speed (rad/sec)
	double delay; 			// Network delay
	double*** Uti; 			// utility
	double*** Cost; 		// Cost
};
/* Throughput estimation */
struct Adapt
{
	double alpha;				// weight values
	double thrp_est_margin;		// safety margin 
	int THRP_EST_METHOD;
	int VP_EST_METHOD;
	int TILE_SELECT_METHOD;	
};

/* head trace */
struct Headtrace
{
	int OFFSET; 				// Number of frames to skip at beginning
	int NO_TRACE; 				// Number of traces
	int TRACE_ID;
	string trace_loc;			// Traces' location
	int*** frame_vp;			// Frames' viewport location
};

class Metadata
{
public:
	Viewport vp;
	Video video;
	BellLab bellLab;
	Adapt adapt;
	Headtrace trace;
	Metadata();
	Metadata(const char* log);
	// Metadata& operator=(Metadata& other);
	~Metadata();
	void load_config_info(const char* filename);
	void print(void);
	void load_video_info(string data_loc);
	void load_head_trace(string trace_loc);
	void load_visible_mask(char *filename);
	void import_matrix_from_txt_file(const char* filename_X, vector <double>& v, int& rows, int& cols);
	int ReadNumbers( const string & s, vector <double> & v );
};
#endif
