#include "Metadata.h"
#include <stdio.h>
#include <cmath>
#define _USE_MATH_DEFINES
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include "common.h"
using namespace std;

Metadata::Metadata(const char* cfg){
	char fname[200];
	load_config_info(cfg);
	/* load tile info */
	load_video_info(video.data_loc);
	// /* load head trace info */
	// loadHeadTrace(filename);
	/* load visible mask */
	// sprintf(fname, "%svmask_%d_face_%dx%d_FoV_90.txt", vp.vmask_data_loc.c_str(), vp.No_face, vp.No_tile_h, vp.No_tile_v);
	sprintf(fname, "%svmask_%d_face_%dx%d_FoV_90.txt", vp.vmask_data_loc.c_str(), vp.No_face, vp.No_tile_h, vp.No_tile_v);
	load_visible_mask(fname);
	load_head_trace(trace.trace_loc);
	// exit(1);
}
Metadata::Metadata(){

}
// Metadata& Metadata::operator=(Metadata& other){
// 	using std::swap;
// 	// swap(vp, other.vp);
// 	// swap(video, other.video);
// 	// swap(bellLab, other.bellLab);
// 	// swap(adapt, other.adapt);
// 	// swap(trace, other.trace);
// 	vp = other.vp;
// 	video = other.video;
// 	bellLab = other.bellLab;
// 	adapt = other.adapt;
// 	trace = other.trace;
// 	return *this;
// }
Metadata::~Metadata(){
}
void Metadata::import_matrix_from_txt_file(const char* filename_X, vector <double>& v, int& rows, int& cols){
	ifstream file_X;
    string line;
    // erase all current elements
    v.erase(v.begin(), v.end());
    cout << "open text file: " << filename_X << endl;
    file_X.open(filename_X);
    if (file_X.is_open())
    {
        int i=0;
        getline(file_X, line);
        cols = ReadNumbers( line, v );
        rows = 1;
        // cout << "cols:" << cols << endl;
		while(getline(file_X, line) != 0){
			ReadNumbers(line, v);
			rows ++;
		}        
        
     
        file_X.close();
    }
    else{
        cout << "#[import_matrix_from_txt_file] Cannot open file " << filename_X << endl;
        exit(-1);
    }
}
int Metadata::ReadNumbers( const string & s, vector <double> & v ){
	istringstream is( s );
    double n;
    while( is >> n ) {
        v.push_back( n );
    }
    return v.size();
}
void Metadata::load_video_info(string data_loc){
	int f,t,tid = 0;
	int i,j,k;
	int seg_id;
	printf("#[load_video_info]:\n");
	int*** tile_frame_size = init3dArrayInt(vp.No_tile, video.NO_VER, video.NO_FRAME);		// frames' sizes in bits
	double*** tile_frame_psnr = init3dArrayDouble(vp.No_tile, video.NO_VER, video.NO_FRAME);	// frames' psnr
	double** frame_size = init2dArrayDouble(video.NO_VER, video.NO_FRAME);
	double** frame_psnr = init2dArrayDouble(video.NO_VER, video.NO_FRAME);
	// CUBE
	double** DASH_BR2 = init2dArrayDouble(video.NO_VER, video.NO_SEG_FULL);
	double** DASH_PSNR2 = init2dArrayDouble(video.NO_VER, video.NO_SEG_FULL);
	double** DASH_MSE2 = init2dArrayDouble(video.NO_VER, video.NO_SEG_FULL);
	// Cube
	double** CMP_DASH_BR2 = init2dArrayDouble(video.NO_VER, video.NO_SEG_FULL);
	double** CMP_DASH_PSNR2 = init2dArrayDouble(video.NO_VER, video.NO_SEG_FULL);
	double** CMP_DASH_MSE2 = init2dArrayDouble(video.NO_VER, video.NO_SEG_FULL);
	// init
	video.TILE_BR = init3dArrayDouble(video.NO_SEG_FULL, vp.No_tile, video.NO_VER);
	video.TILE_PSNR = init3dArrayDouble(video.NO_SEG_FULL, vp.No_tile, video.NO_VER);
	video.TILE_MSE = init3dArrayDouble(video.NO_SEG_FULL, vp.No_tile, video.NO_VER);
	video.TILE_SIZE = init3dArrayDouble(video.NO_SEG_FULL, vp.No_tile, video.NO_VER);
	bellLab.Uti = init3dArrayDouble(video.NO_SEG_FULL, vp.No_tile, video.NO_VER);
	bellLab.Cost = init3dArrayDouble(video.NO_SEG_FULL, vp.No_tile, video.NO_VER);
	// DASH-ERP
	video.DASH_BR = init2dArrayDouble(video.NO_SEG_FULL, video.NO_VER);
	video.DASH_PSNR = init2dArrayDouble(video.NO_SEG_FULL, video.NO_VER);
	video.DASH_MSE = init2dArrayDouble(video.NO_SEG_FULL, video.NO_VER);
	// DASH-CUBE
	video.CMP_DASH_BR = init2dArrayDouble(video.NO_SEG_FULL, video.NO_VER);
	video.CMP_DASH_PSNR = init2dArrayDouble(video.NO_SEG_FULL, video.NO_VER);
	video.CMP_DASH_MSE = init2dArrayDouble(video.NO_SEG_FULL, video.NO_VER);
	//
	video.BR = new double[video.NO_VER];
	video.PSNR = new double[video.NO_VER];
	//
	char tile_info_path[200];
	sprintf(tile_info_path, "%s%df_%dx%d/%dframe/", video.data_loc.c_str(), vp.No_face, vp.No_tile_h, vp.No_tile_v, video.INTERVAL);
	string s;
	char buff[100];
	vector <double> v;
	int rows;
	int cols;
	int N = trace.OFFSET / video.INTERVAL;
	/* Load tiles' info */
	for(tid=0; tid < vp.No_tile; tid++){
		get_face_tid(vp.No_face, vp.No_tile_h, vp.No_tile_v, tid, &f, &t);
		sprintf(buff, "%sf%d_t%d.txt", tile_info_path, f, t);
		s = buff;
		cout << s << endl;
		import_matrix_from_txt_file(buff, v, rows, cols);
		for(i=0; i < video.NO_SEG; i++){
			for(j=0; j < video.NO_VER; j++){
				// cout << i << " " << j << " " << tid <<endl;
				video.TILE_BR[i][tid][j] = v[(N+i)* cols + 3*j];
				video.TILE_PSNR[i][tid][j] = v[(N+i)* cols + 3*j + 1];
				video.TILE_MSE[i][tid][j] = v[(N+i)* cols + 3*j + 2];
				video.TILE_SIZE[i][tid][j] = video.TILE_BR[i][tid][j] * video.SD * video.INTERVAL;
				//
				if(j==0){
					bellLab.Uti[i][tid][j] = video.TILE_MSE[i][tid][j] - 65025;
					bellLab.Cost[i][tid][j] = video.TILE_BR[i][tid][j];
				}else{
					bellLab.Uti[i][tid][j] = video.TILE_MSE[i][tid][j-1] - video.TILE_MSE[i][tid][j];
					bellLab.Cost[i][tid][j] = video.TILE_BR[i][tid][j] - video.TILE_BR[i][tid][j-1];
				}
			}
		}

		// repeat video sequence
		if(video.NO_SEG < video.NO_SEG_FULL){
			for(i=video.NO_SEG; i < video.NO_SEG_FULL; i++){
				// cout << "i:" << i << endl;
				for(k=1; k <= video.NO_VER; k++){
					//System.out.printf("k=%d\n", k);
					video.TILE_BR[i][tid][k-1] = video.TILE_BR[i%video.NO_SEG][tid][k-1];
					video.TILE_PSNR[i][tid][k-1] = video.TILE_PSNR[i%video.NO_SEG][tid][k-1];
					video.TILE_MSE[i][tid][k-1] = video.TILE_MSE[i%video.NO_SEG][tid][k-1];
					video.TILE_SIZE[i][tid][k-1] = video.TILE_SIZE[i%video.NO_SEG][tid][k-1];
					//
					bellLab.Uti[i][tid][k-1] = bellLab.Uti[i%video.NO_SEG][tid][k-1];
					bellLab.Cost[i][tid][k-1] = bellLab.Cost[i%video.NO_SEG][tid][k-1];
				}
			}
		}
	}
	/* load DASH segments info */
	// ERP
	sprintf(buff,"%snon-tile/DASH_seg_%dframe.txt", video.data_loc.c_str(), video.INTERVAL);
	import_matrix_from_txt_file(buff, v, rows, cols);
	for(i=0; i < video.NO_SEG; i++){
		for(j=0; j < video.NO_VER; j++){
			video.DASH_BR[i][j] = v[(N+i) * cols + 3*j];
			video.DASH_PSNR[i][j] = v[(N+i) * cols + 3*j + 1];
			video.DASH_MSE[i][j] = v[(N+i) * cols + 3*j + 2];
			//
			DASH_BR2[j][i] = v[(N+i) * cols + 3*j];
			DASH_PSNR2[j][i] = v[(N+i) * cols + 3*j + 1];
			DASH_MSE2[j][i] = v[(N+i) * cols + 3*j + 2];
		}
	}
	printf("DASH_BR[0][0]=%.2f\tDASH_BR[5][4]=%.2f\n", video.DASH_BR[0][0], video.DASH_BR[5][4]);
	// exit(1);
	// repeat until the end of the streaming session
	if(video.NO_SEG < video.NO_SEG_FULL){
		for(i=video.NO_SEG; i < video.NO_SEG_FULL; i++){
			for(j=0; j < video.NO_VER; j++){
				video.DASH_BR[i][j] = video.DASH_BR[i%video.NO_SEG][j];
				video.DASH_PSNR[i][j] = video.DASH_PSNR[i%video.NO_SEG][j];
				video.DASH_MSE[i][j] = video.DASH_MSE[i%video.NO_SEG][j];
			}
		}
	}
	// CMP
	sprintf(buff,"%s6f_DASH/%dframe/DASH.txt", video.data_loc.c_str(), video.INTERVAL);
	import_matrix_from_txt_file(buff, v, rows, cols);
	for(i=0; i < video.NO_SEG; i++){
		for(j=0; j < video.NO_VER; j++){
			video.CMP_DASH_BR[i][j] = v[(N+i) * cols + 3*j];
			video.CMP_DASH_PSNR[i][j] = v[(N+i) * cols + 3*j + 1];
			video.CMP_DASH_MSE[i][j] = v[(N+i) * cols + 3*j + 2];
			//
			CMP_DASH_BR2[j][i] = v[(N+i) * cols + 3*j];
			CMP_DASH_PSNR2[j][i] = v[(N+i) * cols + 3*j + 1];
			CMP_DASH_MSE2[j][i] = v[(N+i) * cols + 3*j + 2];
		}
	}
	printf("DASH_BR[0][0]=%.2f\tDASH_BR[5][4]=%.2f\n", video.CMP_DASH_BR[0][0], video.CMP_DASH_BR[5][4]);
	// exit(1);
	// repeat until the end of the streaming session
	if(video.NO_SEG < video.NO_SEG_FULL){
		for(i=video.NO_SEG; i < video.NO_SEG_FULL; i++){
			for(j=0; j < video.NO_VER; j++){
				video.CMP_DASH_BR[i][j] = video.CMP_DASH_BR[i%video.NO_SEG][j];
				video.CMP_DASH_PSNR[i][j] = video.CMP_DASH_PSNR[i%video.NO_SEG][j];
				video.CMP_DASH_MSE[i][j] = video.CMP_DASH_MSE[i%video.NO_SEG][j];
			}
		}
	}
	//
	/* load frames' info and calculate R-D values */
	for(tid = 0; tid < vp.No_tile; tid++){
		get_face_tid(vp.No_face, vp.No_tile_h, vp.No_tile_v, tid, &f, &t);
		sprintf(buff,"%s%df_%dx%d/tile_size_per_frame/f%d_t%d.txt", video.data_loc.c_str(), vp.No_face, vp.No_tile_h, vp.No_tile_v, f, t);
		import_matrix_from_txt_file(buff, v, rows, cols);
		for(i=0; i < video.NO_VER; i++){
			for(j=0; j < video.NO_FRAME_ORIGIN; j++){
				// printf("# %d %d \n", j, video.NO_FRAME);
				tile_frame_size[tid][i][j] = v[3 * video.NO_VER * (j+trace.OFFSET) + 3*i];
				tile_frame_psnr[tid][i][j] = v[3 * video.NO_VER * (j+trace.OFFSET) + 3*i + 1];
			}
			if(video.NO_FRAME_ORIGIN < video.NO_FRAME){
				for(j=video.NO_FRAME_ORIGIN; j < video.NO_FRAME; j++){
					tile_frame_size[tid][i][j] = tile_frame_size[tid][i][j%video.NO_FRAME_ORIGIN];
					tile_frame_psnr[tid][i][j] = tile_frame_psnr[tid][i][j%video.NO_FRAME_ORIGIN];
				}
			}
		}
	}
	printf("#[load_video_info]: %d %.2f\n", tile_frame_size[0][0][0], tile_frame_psnr[0][0][0]);
	/* calculate R-D values */
	sprintf(buff,"result/RD_%df%dx%d.txt", vp.No_face, vp.No_tile_h, vp.No_tile_v);
	FILE* fout = fopen(buff, "w");
	fprintf(fout, "Bitrate\tPSNR\n");
	for(i=0; i < video.NO_VER; i++){
		for(j=0; j < video.NO_FRAME; j++){
			frame_size[i][j] = 0;
			frame_psnr[i][j] = 0;
			for(tid = 0; tid < vp.No_tile; tid++){
				frame_size[i][j] += tile_frame_size[tid][i][j];
				frame_psnr[i][j] += (255*255)/pow(10, tile_frame_psnr[tid][i][j]/10.0);
			}
			frame_psnr[i][j] = 10 * log10(255*255*vp.No_tile/frame_psnr[i][j]);
		}
		video.BR[i] = sum(frame_size[i], video.NO_FRAME)/1000.0/(video.NO_FRAME*1.0/video.FPS);
		video.PSNR[i] = avg(frame_psnr[i], video.NO_FRAME);
		fprintf(fout, "%.2f\t%.2f\n", video.BR[i], video.PSNR[i]);
		printf("%.2f\t%.2f\t%.2f\t%.2f\n", avg(DASH_BR2[i], video.NO_SEG_FULL), avg(DASH_PSNR2[i], video.NO_SEG_FULL),avg(CMP_DASH_BR2[i], video.NO_SEG_FULL),avg(CMP_DASH_PSNR2[i], video.NO_SEG_FULL));
		// fprintf(fout, "%.2f\t%.2f\t%.2f\t%.2f\n", video.BR[i], video.PSNR[i], avg(DASH_BR2[i], video.NO_SEG_FULL), avg(DASH_PSNR2[i], video.NO_SEG_FULL));
	}
	fclose(fout);
	// exit(1);
}
void Metadata::load_head_trace(string trace_loc){
	printf("[load_head_trace]:\n");
	int i,j,k;
	char fname[200];
	vector <double> v;
	int rows;
	int cols;
	trace.frame_vp = init3dArrayInt(trace.NO_TRACE, video.NO_FRAME, 2);
	for(k=0; k < trace.NO_TRACE; k++){
		sprintf(fname, "%sxyz_vid_5_uid_%d.txt", trace_loc.c_str(), trace.TRACE_ID);
		import_matrix_from_txt_file(fname, v, rows, cols);
		printf("[load_head_trace]: finished reading file\n");
		for(i=0; i < video.NO_FRAME; i++){
			for(j=0; j < cols; j++){
				trace.frame_vp[k][i][j] = v[(i+trace.OFFSET) * cols + j];
			}
		}
	}
}
void Metadata::load_visible_mask(char* filename){
	int i,j,k;
	vector <double> v;
	int rows;
	int cols;
	int phi_num = 360;
	int theta_num = 181;
	printf("#[load_visible_mask]:\n");
	import_matrix_from_txt_file(filename, v, rows, cols);
	int NO_VMASK_SAMPLE = rows;
	vp.vmask = init3dArrayInt(phi_num, theta_num, vp.No_tile);
	vp.pixel = init3dArrayInt(phi_num, theta_num, vp.No_tile);
	for(i=0; i < NO_VMASK_SAMPLE; i++){
		// printf("%d %d\n", int(v[i*cols]), int(v[i*cols + 1]) + 90);
		for(j=0; j < vp.No_tile; j++){
			if(vp.No_face == 1 || vp.No_face == 2) // ERP
				vp.pixel[int(v[i*cols])][int(v[i*cols + 1]) + 90][j] = v[i * cols + 2*j + 3];
			else // CMP
				vp.pixel[int(v[i*cols])][int(v[i*cols + 1]) + 90][j] = v[i * cols + j + 2];
			if(vp.pixel[int(v[i*cols])][int(v[i*cols + 1]) + 90][j] > 0)
				vp.vmask[int(v[i*cols])][int(v[i*cols + 1]) + 90][j] = 1;
			else
				vp.vmask[int(v[i*cols])][int(v[i*cols + 1]) + 90][j] = 0;
		}
	}
}
void Metadata::load_config_info(const char* filename){
	printf("#[read_config_info]:\n");
	string s;
	string delimiter = "=";
	string comment = "#";
	string key;
	string val_str;
	double val;
	size_t pos_deli = 0;
	size_t pos_comm;
	ifstream infile(filename);
	if(infile == NULL){
		cout << "Cannot open file " << filename << endl;
	}else{
		cout << "Reading file" << filename << endl;
	}
	while(std::getline(infile, s)){
		// cout << s << endl;
		if((pos_deli=s.find(delimiter)) != std::string::npos){
			key = s.substr(0, pos_deli-1);
			pos_comm = s.find(comment);
			if(pos_comm != std::string::npos){
				val_str = s.substr(pos_deli + 2, pos_comm-pos_deli-2);
			}else{
				val_str = s.substr(pos_deli + 2, s.length());
			}
		}
		// printf("%s\n", val_str.c_str());
		// viewport
		if(key.compare("Fh")==0) vp.Fh = std::stod(val_str) * M_PI / 180;
		if(key.compare("Fv")==0) vp.Fv= std::stod(val_str) * M_PI / 180;
		if(key.compare("vp_W")==0) vp.vp_W = (int) std::stod(val_str);
		if(key.compare("vp_H")==0) vp.vp_H = (int) std::stod(val_str);
		if(key.compare("erp_W")==0) vp.erp_W = (int) std::stod(val_str);
		if(key.compare("erp_H")==0) vp.erp_H = (int) std::stod(val_str);
		if(key.compare("No_face")==0) vp.No_face = (int) std::stod(val_str);
		if(key.compare("face_W") == 0) vp.face_W = (int) std::stod(val_str);
		if(key.compare("face_H") == 0) vp.face_H = (int) std::stod(val_str);
		if(key.compare("No_tile_h")==0) vp.No_tile_h = (int) std::stod(val_str);
		if(key.compare("No_tile_v")==0) vp.No_tile_v = (int) std::stod(val_str);
		if(key.compare("vmask_data_loc") == 0) vp.vmask_data_loc = val_str;
		// Video
		if(key.compare("NO_VER")==0) video.NO_VER = (int) std::stod(val_str);
		if(key.compare("FPS")==0) video.FPS = (int) std::stod(val_str);
		if(key.compare("BUFF")==0) video.BUFF = std::stod(val_str);
		if(key.compare("INTERVAL")==0) video.INTERVAL = (int) std::stod(val_str);	
		if(key.compare("NO_FRAME") == 0) video.NO_FRAME = (int) std::stod(val_str);
		if(key.compare("NO_FRAME_ORIGIN") == 0) video.NO_FRAME_ORIGIN = (int) std::stod(val_str);
		if(key.compare("data_loc") == 0) video.data_loc = val_str;
		if(key.compare("name") == 0) video.name = val_str;
		// Adapt
		if(key.compare("thrp_est_margin")==0) adapt.thrp_est_margin = std::stod(val_str);
		if(key.compare("alpha")==0) adapt.alpha = std::stod(val_str);
		if(key.compare("THRP_EST_METHOD")==0) adapt.THRP_EST_METHOD = (int)std::stod(val_str);
		if(key.compare("VP_EST_METHOD")==0) adapt.VP_EST_METHOD = (int)std::stod(val_str);
		if(key.compare("TILE_SELECT_METHOD")==0) adapt.TILE_SELECT_METHOD = (int)std::stod(val_str);
		// BellLab
		if(key.compare("speed") == 0) bellLab.speed = std::stod(val_str);
		if(key.compare("delay") == 0) bellLab.delay = std::stod(val_str);
		// Head trace
		if(key.compare("OFFSET") == 0) trace.OFFSET = (int) std::stod(val_str);
		if(key.compare("NO_TRACE") == 0) trace.NO_TRACE = (int) std::stod(val_str);
		if(key.compare("TRACE_ID") == 0) trace.TRACE_ID = (int) std::stod(val_str);
		if(key.compare("trace_loc") == 0) trace.trace_loc = val_str;
	}
	// vp
	// if(vp.No_face == 1 || vp.No_face == 6){
		vp.tile_W = vp.face_W / vp.No_tile_h;
		vp.tile_H = vp.face_H / vp.No_tile_v;
		vp.No_tile = vp.No_face * vp.No_tile_h * vp.No_tile_v;
	// }
	if(vp.No_face == 2){
		vp.No_tile = 2 + vp.No_tile_h * vp.No_tile_v;
	}
	//
	video.NO_SEG = video.NO_FRAME_ORIGIN / video.INTERVAL;
	video.NO_SEG_FULL = video.NO_FRAME / video.INTERVAL;
	video.SD = video.INTERVAL * 1.0 / video.FPS;
	video.SESS_DUR = video.NO_SEG * video.SD;
}
void Metadata::print(void){
	printf("[Viewport]:\n");
	printf("Fh:%.2f\nFv:%.2f\nvp_W:%d\nvp_H:%d\nerp_W:%d\nerp_H:%d\n", vp.Fh, vp.Fv, vp.vp_W, vp.vp_H, vp.erp_W, vp.erp_H);
	printf("tile_W:%d\ntile_H:%d\nNo_face:%d\nface_W:%d\nface_H:%d\n", vp.tile_W, vp.tile_H, vp.No_face, vp.face_W, vp.face_H);
	printf("No_tile_h:%d\nNo_tile_v:%d\nNo_tile:%d\nvmask_loc:%s\n", vp.No_tile_h, vp.No_tile_v, vp.No_tile, vp.vmask_data_loc.c_str());
	printf("[Video]:\n");
	printf("name:%s\ndata_loc:%s\nSD:%.2f\nSESS_DUR:%.2f\nFPS:%d\n", video.name.c_str(), video.data_loc.c_str(), video.SD, video.SESS_DUR, video.FPS);
	printf("NO_SEG:%d\nNO_SEG_ORIGIN:%d\nBUFF:%d\nINTERVAL:%d\nNO_FRAME:%d\nNO_FRAME_ORIGIN:%d\n", video.NO_SEG_FULL, video.NO_SEG, video.BUFF, video.INTERVAL, video.NO_FRAME, video.NO_FRAME_ORIGIN);
	printf("NO_VER:%d\n", video.NO_VER);
	printf("[BellLab]:\n");
	printf("speed:%.2f\ndelay:%.2f\n", bellLab.speed, bellLab.delay);
	printf("[Adapt]:\n");
	printf("thrp_est_margin:%.2f\n", adapt.thrp_est_margin);
}