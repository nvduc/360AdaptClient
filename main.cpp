#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "Projection.h"
#include "Metadata.h"
#include "AdaptLogic.h"
#include "common.h"
using namespace std;
struct run_cfg
{
	std::vector<int> bw;
	std::vector<int> method;
	std::vector<int> trace;
	int no_trace;
	int no_method;
	int no_bw;
};
void write_result(AdaptLogic, int, int);
void write_result2(AdaptLogic adaptLogic, int N, int BW);
run_cfg load_run_cfg(const char* cfg);

#ifdef DEBUG1
int main(int argc, char const *argv[])
{
	if(argc != 6){
		printf("usage: ./main NO_FACE NO_TILE_H NO_TILE_W tid ext_w\n");
		exit(1);
	}
	int NO_FACE = atoi(argv[1]);
	int NO_TILE_H = atoi(argv[2]);
	int NO_TILE_W = atoi(argv[3]);
	int tid = atoi(argv[4]);
	int ext_w = atoi(argv[5]);
	int N;
	int *adjTile;
	int* vmask_ext;
	int vmask[] = { /*Face 0*/ 0, 0, 0, 0, 0, 1, 1, 0, 0,
					/*Face 1*/ 1, 1, 0, 0, 0, 0, 0, 0, 0,
					/*Face 2*/ 0, 0, 0, 0, 0, 0, 0, 0, 0,
					/*Face 3*/ 0, 0, 0, 0, 0, 0, 0, 0, 0,
					/*Face 4*/ 0, 0, 0, 0, 0, 0, 0, 0, 0,
					/*Face 5*/ 0, 0, 0, 0, 0, 0, 0, 0, 0};
	// adjTile = getAdjacentTile(NO_FACE, NO_TILE_H, NO_TILE_W, tid, &N);
	// showArray(adjTile, N);
	vmask_ext = extVmask(vmask, NO_FACE, NO_TILE_H, NO_TILE_W, ext_w);
	showVmask(vmask_ext, NO_FACE, NO_TILE_H, NO_TILE_W);
	return 0;
}
#endif
#ifdef DEBUG2
int main(int argc, char const *argv[])
{
	int index;
	int N = 3;
	int trace_id;
	int method_id;
	int bw_id;
	int BW, METHOD;
	FILE* fout;
	char fname[1024];
	// int BW_ARR[] = {2100, 4100, 7600, 13500,23200, 37900};
	run_cfg cfg;
	if(argc != 3){
		printf("usage: ./main cfg_file run.cfg\n");
		exit(1);
	}
	cfg = load_run_cfg(argv[2]);
	Metadata metadata (argv[1]);
	// metadata.print();
	AdaptLogic adaptLogic(metadata);
	// adaptLogic.metadata.print();
	N = metadata.video.NO_SEG_FULL;
	// for(trace_id = 0; trace_id < cfg.no_trace; trace_id ++)
	// 	for(method_id = 0; method_id < cfg.no_method; method_id ++)
	// 		printf("trace: %d method: %d\n", cfg.trace[trace_id], cfg.method[method_id]);
	//
	/* open file to log the results */
	sprintf(fname, "result/video_%s/%df_%dx%d_trace_0.txt", metadata.video.name.c_str(), metadata.vp.No_face, metadata.vp.No_tile_h, metadata.vp.No_tile_v);
	fout = fopen(fname, "w");
	fprintf(fout, "BW\t");
	for(method_id = 0; method_id < cfg.no_method; method_id++){
		fprintf(fout, "%s\t%s\t%s\t", adaptLogic.get_method_name(cfg.method[method_id]), adaptLogic.get_method_name(cfg.method[method_id]), adaptLogic.get_method_name(cfg.method[method_id]));
	}
	fprintf(fout, "\n");
	//
	for(bw_id = 0; bw_id < cfg.no_bw; bw_id++){
		BW = cfg.bw[bw_id];
		fprintf(fout, "%d\t", BW);
		for(method_id = 0; method_id < cfg.no_method; method_id++){
			adaptLogic.metadata.adapt.TILE_SELECT_METHOD = cfg.method[method_id];
			for(index=0; index < N; index++){
				adaptLogic.thrp.seg_thrp[index] = BW;
				adaptLogic.get_next_segment(index);
				// printf("#[main] segment #%d\n", index);
				// showVmask(adaptLogic.tile_ver[index], adaptLogic.metadata.vp.No_face, adaptLogic.metadata.vp.No_tile_h, adaptLogic.metadata.vp.No_tile_v);
			}
			/* calculate metrics */
			adaptLogic.calc_result();
			/* write results to files */
			write_result(adaptLogic, N, BW);
			//
			fprintf(fout, "%.2f\t%.2f\t%.2f\t", avg(adaptLogic.vpsnr, adaptLogic.metadata.video.NO_FRAME), avg(adaptLogic.seg_br, adaptLogic.metadata.video.NO_SEG_FULL),100*avg(adaptLogic.seg_br, adaptLogic.metadata.video.NO_SEG_FULL)/BW);
		}
		fprintf(fout, "\n");
	}
	fclose(fout);
	return 0;
}
#endif
#ifdef DEBUG3
int main(int argc, char const *argv[])
{
	int index, i;
	int N = 3;
	int trace_id;
	int method_id;
	int bw_id, phi, theta;
	int BW, METHOD;
	FILE* fout;
	char fname[1024];
	// int BW_ARR[] = {2100, 4100, 7600, 13500,23200, 37900};
	run_cfg cfg;
	if(argc != 3){
		printf("usage: ./main cfg_file run.cfg\n");
		exit(1);
	}
	cfg = load_run_cfg(argv[2]);
	Metadata metadata (argv[1]);
	// metadata.print();
	AdaptLogic adaptLogic(metadata);
	// adaptLogic.metadata.print();
	N = metadata.video.NO_SEG_FULL;
	// for(trace_id = 0; trace_id < cfg.no_trace; trace_id ++)
	// 	for(method_id = 0; method_id < cfg.no_method; method_id ++)
	// 		printf("trace: %d method: %d\n", cfg.trace[trace_id], cfg.method[method_id]);
	//
	for(bw_id = 0; bw_id < cfg.no_bw; bw_id++){
		//				
		BW = cfg.bw[bw_id];
		/* open file to log the results */
		sprintf(fname, "result/video_%s/%df_%dx%d_trace_0_fix_BW_%d", metadata.video.name.c_str(), metadata.vp.No_face, metadata.vp.No_tile_h, metadata.vp.No_tile_v, BW);
		fout = fopen(fname, "w");
		fprintf(fout, "phi\ttheta\t");
		for(method_id = 0; method_id < cfg.no_method; method_id++){
			fprintf(fout, "%s\t%s\t%s\t", adaptLogic.get_method_name(cfg.method[method_id]), adaptLogic.get_method_name(cfg.method[method_id]), adaptLogic.get_method_name(cfg.method[method_id]));
		}
		fprintf(fout, "\n");
		for(phi = 0; phi < 360; phi += 15){
		// for(phi = 0; phi < 1; phi += 15){
			for(theta = 0; theta < 1; theta ++){
					fprintf(fout, "%d\t%d\t", phi - 180, theta);
					// overwrite head trace
					for(i=0; i < adaptLogic.metadata.video.NO_FRAME; i++){
						adaptLogic.metadata.trace.frame_vp[0][i][0] = phi - 180;
						// adaptLogic.metadata.trace.frame_vp[0][i][0] = phi;
						adaptLogic.metadata.trace.frame_vp[0][i][1] = theta;
					}
					for(method_id = 0; method_id < cfg.no_method; method_id++){
						adaptLogic.metadata.adapt.TILE_SELECT_METHOD = cfg.method[method_id];
						for(index=0; index < N; index++){
							adaptLogic.thrp.seg_thrp[index] = BW;
							adaptLogic.get_next_segment(index);
							// printf("#[main] segment #%d\n", index);
							// showVmask(adaptLogic.tile_ver[index], adaptLogic.metadata.vp.No_face, adaptLogic.metadata.vp.No_tile_h, adaptLogic.metadata.vp.No_tile_v);
						}
						/* calculate metrics */
						adaptLogic.calc_result();
						/* write results to files */
						// write_result(adaptLogic, N, BW);
						//
						fprintf(fout, "%.2f\t%.2f\t%.2f\t", avg(adaptLogic.vpsnr, adaptLogic.metadata.video.NO_FRAME), avg(adaptLogic.seg_br, adaptLogic.metadata.video.NO_SEG_FULL),100*avg(adaptLogic.seg_br, adaptLogic.metadata.video.NO_SEG_FULL)/BW);
					}
					fprintf(fout, "\n");
				}
			}
		}
	// fclose(fout);
	return 0;
}
#endif
#ifdef DEBUG4
int main(int argc, char const *argv[])
{
	int index, i;
	int N = 3;
	int trace_id;
	int method_id;
	int bw_id, phi, theta;
	int BW, METHOD;
	// double veloc = 1; /* angular velocity (degree/frame) */
	double veloc = 0.5; /* angular velocity (degree/frame) */
	// double veloc = 1; /* angular velocity (degree/frame) */
	FILE* fout;
	char fname[1024];
	// int BW_ARR[] = {2100, 4100, 7600, 13500,23200, 37900};
	run_cfg cfg;
	if(argc != 4){
		printf("usage: ./main cfg_file run.cfg velocity\n");
		exit(1);
	}
	cfg = load_run_cfg(argv[2]);
	Metadata metadata (argv[1]);
	veloc = atof(argv[3]);
	// metadata.print();
	AdaptLogic adaptLogic(metadata);
	// adaptLogic.metadata.print();
	N = metadata.video.NO_SEG_FULL;
	// for(trace_id = 0; trace_id < cfg.no_trace; trace_id ++)
	// 	for(method_id = 0; method_id < cfg.no_method; method_id ++)
	// 		printf("trace: %d method: %d\n", cfg.trace[trace_id], cfg.method[method_id]);
	//
	for(bw_id = 0; bw_id < cfg.no_bw; bw_id++){
		//				
		BW = cfg.bw[bw_id];
		/* open file to log the results */
		sprintf(fname, "result/video_%s/%df_%dx%d_trace_0_simple_trace_veloc_%.2f_BW_%d", metadata.video.name.c_str(), metadata.vp.No_face, metadata.vp.No_tile_h, metadata.vp.No_tile_v, veloc, BW);
		fout = fopen(fname, "w");
		fprintf(fout, "phi\ttheta\t");
		for(method_id = 0; method_id < cfg.no_method; method_id++){
			fprintf(fout, "%s\t%s\t%s\t", adaptLogic.get_method_name(cfg.method[method_id]), adaptLogic.get_method_name(cfg.method[method_id]), adaptLogic.get_method_name(cfg.method[method_id]));
		}
		fprintf(fout, "\n");
		for(phi = 0; phi < 1; phi += 15){
			for(theta = 0; theta < 1; theta ++){
					fprintf(fout, "%d\t%d\t", phi - 180, theta);
					// overwrite head trace
					int i_0 = adaptLogic.metadata.video.NO_FRAME/2 - 1;
					for(i=0; i < adaptLogic.metadata.video.NO_FRAME; i++){
						if( i < adaptLogic.metadata.video.NO_FRAME/2)
							adaptLogic.metadata.trace.frame_vp[0][i][0] = phi - 180;
						else{
							adaptLogic.metadata.trace.frame_vp[0][i][0] = (int) (adaptLogic.metadata.trace.frame_vp[0][i_0][0] + (i - i_0) * veloc);
							if(adaptLogic.metadata.trace.frame_vp[0][i][0] > 180)
								adaptLogic.metadata.trace.frame_vp[0][i][0] -= 360;
						}
						adaptLogic.metadata.trace.frame_vp[0][i][1] = theta;
					}
					for(method_id = 0; method_id < cfg.no_method; method_id++){
						adaptLogic.metadata.adapt.TILE_SELECT_METHOD = cfg.method[method_id];
						for(index=0; index < N; index++){
							adaptLogic.thrp.seg_thrp[index] = adaptLogic.metadata.video.DASH_BR[index][3] + 500;
							adaptLogic.get_next_segment(index);
							// printf("#[main] segment #%d\n", index);
							// showVmask(adaptLogic.tile_ver[index], adaptLogic.metadata.vp.No_face, adaptLogic.metadata.vp.No_tile_h, adaptLogic.metadata.vp.No_tile_v);
						}
						/* calculate metrics */
						adaptLogic.calc_result();
						/* write results to files */
						write_result2(adaptLogic, N, BW);
						//
						fprintf(fout, "%.2f\t%.2f\t%.2f\t", avg(adaptLogic.vpsnr, adaptLogic.metadata.video.NO_FRAME), avg(adaptLogic.seg_br, adaptLogic.metadata.video.NO_SEG_FULL),100*avg(adaptLogic.seg_br, adaptLogic.metadata.video.NO_SEG_FULL)/BW);
					}
					fprintf(fout, "\n");
				}
			}
		}
	fclose(fout);
	return 0;
}
#endif
void write_result(AdaptLogic adaptLogic, int N, int BW){
	int jj=0; // frame id
	double BR = 0, avg_vp_psnr = 0;
	int tid;
	int i,j,k;
	char fname[1024];
	FILE* log_tile_ver, *log_tile_ver_linear, *log_dec, *log_frame_psnr;
	FILE* fVPSNR;
	// generate log files
	sprintf(fname, "result/video_%s/log_tile_ver_TRACE_%d_BW_%d_METHOD_%d_VP_EST_METHOD_%d_INTER_%d_BUFF_%d_%df_%dx%d.txt",adaptLogic.metadata.video.name.c_str(), 0, BW, adaptLogic.metadata.adapt.TILE_SELECT_METHOD, adaptLogic.metadata.adapt.VP_EST_METHOD,adaptLogic.metadata.video.INTERVAL, adaptLogic.metadata.video.BUFF, adaptLogic.metadata.vp.No_face, adaptLogic.metadata.vp.No_tile_h, adaptLogic.metadata.vp.No_tile_v);
	printf("#%s\n", fname);
	log_tile_ver = fopen(fname, "w");
	sprintf(fname, "result/video_%s/log_tile_ver_linear_TRACE_%d_BW_%d_METHOD_%d_VP_EST_METHOD_%d_INTER_%d_BUFF_%d_%df_%dx%d.txt",adaptLogic.metadata.video.name.c_str(), 0, BW, adaptLogic.metadata.adapt.TILE_SELECT_METHOD, adaptLogic.metadata.adapt.VP_EST_METHOD,adaptLogic.metadata.video.INTERVAL, adaptLogic.metadata.video.BUFF, adaptLogic.metadata.vp.No_face, adaptLogic.metadata.vp.No_tile_h, adaptLogic.metadata.vp.No_tile_v);
	log_tile_ver_linear = fopen(fname, "w");
	printf("#%s\n", fname);
	sprintf(fname, "result/video_%s/log_dec_TRACE_%d_BW_%d_METHOD_%d_VP_EST_METHOD_%d_INTER_%d_BUFF_%d_%df_%dx%d.txt",adaptLogic.metadata.video.name.c_str(), 0, BW, adaptLogic.metadata.adapt.TILE_SELECT_METHOD, adaptLogic.metadata.adapt.VP_EST_METHOD,adaptLogic.metadata.video.INTERVAL, adaptLogic.metadata.video.BUFF, adaptLogic.metadata.vp.No_face, adaptLogic.metadata.vp.No_tile_h, adaptLogic.metadata.vp.No_tile_v);
	log_dec = fopen(fname, "w");
	sprintf(fname, "result/video_%s/log_frame_TRACE_%d_BW_%d_METHOD_%d_VP_EST_METHOD_%d_INTER_%d_BUFF_%d_%df_%dx%d.txt",adaptLogic.metadata.video.name.c_str(), 0, BW, adaptLogic.metadata.adapt.TILE_SELECT_METHOD,adaptLogic.metadata.adapt.VP_EST_METHOD, adaptLogic.metadata.video.INTERVAL, adaptLogic.metadata.video.BUFF, adaptLogic.metadata.vp.No_face, adaptLogic.metadata.vp.No_tile_h, adaptLogic.metadata.vp.No_tile_v);
	printf("#%s\n", fname);
	log_frame_psnr = fopen(fname, "w");
	//
	if(log_tile_ver == NULL || log_tile_ver_linear == NULL || log_dec == NULL){
	printf("#[init] Cannot open log files\n");
	exit(1);
	}
	if(adaptLogic.metadata.adapt.TILE_SELECT_METHOD > 1){
			fprintf(log_dec, "id\testThrp(kbps)\ttext_width\tbitrate(kbps)\n");
			fprintf(log_frame_psnr, "fid\tdecid\test_vp_psnr\tphi\ttheta\test_phi\test_theta\terr_phi\terr_theta\text_width\n");
			for(int ii=0; ii < N; ii++){
			fprintf(log_tile_ver, "\nseg #%d calcTime: %d(ms)\n", ii, adaptLogic.proc_time[ii]);
			// MxN
			if(adaptLogic.metadata.vp.No_face == 1 || adaptLogic.metadata.vp.No_face == 6){
				for(int i=0; i < adaptLogic.metadata.vp.No_tile_v; i++){
					// version
					for(int k=0; k < adaptLogic.metadata.vp.No_face; k++){
					for(int j=0; j < adaptLogic.metadata.vp.No_tile_h; j++){
					fprintf(log_tile_ver, "%d ",adaptLogic.tile_ver[ii][k*adaptLogic.metadata.vp.No_tile_h*adaptLogic.metadata.vp.No_tile_v + i * adaptLogic.metadata.vp.No_tile_h + j]);
					fprintf(log_tile_ver_linear, "%d ",adaptLogic.tile_ver[ii][k*adaptLogic.metadata.vp.No_tile_h*adaptLogic.metadata.vp.No_tile_v + i * adaptLogic.metadata.vp.No_tile_h + j]);
					}
					fprintf(log_tile_ver, "\t");    
					}
					fprintf(log_tile_ver, "\n");
				}
				fprintf(log_tile_ver, "\n");
				// psnr
				for(int i=0; i < adaptLogic.metadata.vp.No_tile_v; i++){
					for(int k=0; k < adaptLogic.metadata.vp.No_face; k++){
					for(int j=0; j < adaptLogic.metadata.vp.No_tile_h; j++){
					fprintf(log_tile_ver, "%.2f ",adaptLogic.metadata.video.TILE_PSNR[ii][k*adaptLogic.metadata.vp.No_tile_h*adaptLogic.metadata.vp.No_tile_v + i * adaptLogic.metadata.vp.No_tile_h + j][adaptLogic.tile_ver[ii][k*adaptLogic.metadata.vp.No_tile_h*adaptLogic.metadata.vp.No_tile_v + i * adaptLogic.metadata.vp.No_tile_h + j]]);
					}
					fprintf(log_tile_ver, "\t");    
					}
					fprintf(log_tile_ver, "\n");
				}	  
				fprintf(log_tile_ver, "\n");
				// bitrate
				for(int i=0; i < adaptLogic.metadata.vp.No_tile_v; i++){
					for(int k=0; k < adaptLogic.metadata.vp.No_face; k++){
					for(int j=0; j < adaptLogic.metadata.vp.No_tile_h; j++){
					fprintf(log_tile_ver, "%.2f ",adaptLogic.metadata.video.TILE_BR[ii][k*adaptLogic.metadata.vp.No_tile_h*adaptLogic.metadata.vp.No_tile_v + i * adaptLogic.metadata.vp.No_tile_h + j][adaptLogic.tile_ver[ii][k*adaptLogic.metadata.vp.No_tile_h*adaptLogic.metadata.vp.No_tile_v + i * adaptLogic.metadata.vp.No_tile_h + j]]);
					}
					fprintf(log_tile_ver, "\t");    
					}	  
					fprintf(log_tile_ver, "\n");
				}
				fprintf(log_tile_ver, "\n");
				fprintf(log_tile_ver_linear, "\n");
			}
			// 2+MxN (ERP)
			if(adaptLogic.metadata.vp.No_face == 2){
				// version
				fprintf(log_tile_ver, "%d\n",adaptLogic.tile_ver[ii][0]);
				for(int i=0; i < adaptLogic.metadata.vp.No_tile_v; i++){
					for(int j=0; j < adaptLogic.metadata.vp.No_tile_h; j++){
						fprintf(log_tile_ver, "%d\t",adaptLogic.tile_ver[ii][2+i*adaptLogic.metadata.vp.No_tile_h + j]);
					}
					fprintf(log_tile_ver, "\n");
				}
				fprintf(log_tile_ver, "%d\n",adaptLogic.tile_ver[ii][1]);
				fprintf(log_tile_ver, "\n");
				// PSNR
				fprintf(log_tile_ver, "%.2f\n",adaptLogic.metadata.video.TILE_PSNR[ii][0][adaptLogic.tile_ver[ii][0]]);
				for(int i=0; i < adaptLogic.metadata.vp.No_tile_v; i++){
					for(int j=0; j < adaptLogic.metadata.vp.No_tile_h; j++){
						fprintf(log_tile_ver, "%.2f\t",adaptLogic.metadata.video.TILE_PSNR[ii][2+i*adaptLogic.metadata.vp.No_tile_h + j][adaptLogic.tile_ver[ii][2+i*adaptLogic.metadata.vp.No_tile_h + j]]);
					}
					fprintf(log_tile_ver, "\n");
				}
				fprintf(log_tile_ver, "%.2f\n",adaptLogic.metadata.video.TILE_PSNR[ii][1][adaptLogic.tile_ver[ii][1]]);
				fprintf(log_tile_ver, "\n");
				// Bitrate
				fprintf(log_tile_ver, "%.2f\n",adaptLogic.metadata.video.TILE_BR[ii][0][adaptLogic.tile_ver[ii][0]]);
				for(int i=0; i < adaptLogic.metadata.vp.No_tile_v; i++){
					for(int j=0; j < adaptLogic.metadata.vp.No_tile_h; j++){
						fprintf(log_tile_ver, "%.2f\t",adaptLogic.metadata.video.TILE_BR[ii][2+i*adaptLogic.metadata.vp.No_tile_h + j][adaptLogic.tile_ver[ii][2+i*adaptLogic.metadata.vp.No_tile_h + j]]);
					}
					fprintf(log_tile_ver, "\n");
				}
				fprintf(log_tile_ver, "%.2f\n",adaptLogic.metadata.video.TILE_BR[ii][1][adaptLogic.tile_ver[ii][1]]);
				fprintf(log_tile_ver, "\n");
			}
			//
			fprintf(log_dec, "%d\t%.2f\t%.2f\n", ii, adaptLogic.thrp.est_seg_thrp[ii], adaptLogic.proc_time[ii]/1000.0/1000.0);
			// calculate viewport psnr of frames in this interval
			// decide the adaptation result of this segment
			int *tmpVisi;
			// printf("#[write_result] seg: %d\n", ii);
			// printf("est_vp: (%d, %d)\n", adaptLogic->est_vp[ii][0], adaptLogic->est_vp[ii][1]);
			// tmpVisi = adaptLogic->getVisibleTile(adaptLogic->est_vp[ii]);
			// for(i=0; i < rows; i++){
			//   for(j=0; j < cols; j++){
			//   	printf("%d ", tmpVisi[i*cols + j]);
			//   }
			//  printf("\n");
			// }
			for(int k=0; k < adaptLogic.metadata.video.INTERVAL; k++){
			fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\n",ii * adaptLogic.metadata.video.INTERVAL + k, ii, adaptLogic.vpsnr[ii*adaptLogic.metadata.video.INTERVAL +k], adaptLogic.metadata.trace.frame_vp[0][ii*adaptLogic.metadata.video.INTERVAL +k][0], adaptLogic.metadata.trace.frame_vp[0][ii*adaptLogic.metadata.video.INTERVAL +k][1], adaptLogic.vp2.est_frame_vp[ii*adaptLogic.metadata.video.INTERVAL + k][0], adaptLogic.vp2.est_frame_vp[ii*adaptLogic.metadata.video.INTERVAL + k][1], adaptLogic.vp2.est_err[ii][0], adaptLogic.vp2.est_err[ii][1]);
			}

			fflush(log_tile_ver);
			fflush(log_tile_ver_linear);
			fflush(log_dec);
			fflush(log_frame_psnr);
			// if(ii==5)
			// 	exit(1);
			}

	}else{
		fprintf(log_frame_psnr, "fid\tdecid\tver\tbitrate\test_vp_psnr\n");
		for(int ii=0; ii < N; ii++){
		for(int k=0; k < adaptLogic.metadata.video.INTERVAL; k++){
		fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%.2f\n", ii * adaptLogic.metadata.video.INTERVAL + k, ii, adaptLogic.vpsnr[ii*adaptLogic.metadata.video.INTERVAL +k], adaptLogic.tile_ver[ii][0], adaptLogic.metadata.video.DASH_BR[ii][adaptLogic.tile_ver[ii][0]]);
		}
		}
		fflush(log_frame_psnr);
		fclose(log_frame_psnr);
	}
	printf("#[print_result] BW:%d METHOD:%d\t%.2f\n", BW, adaptLogic.metadata.adapt.TILE_SELECT_METHOD, avg_vp_psnr);
}
void write_result2(AdaptLogic adaptLogic, int N, int BW){
	int jj=0; // frame id
	double BR = 0, avg_vp_psnr = 0;
	int tid;
	int i,j,k;
	char fname[1024];
	FILE* log_tile_ver, *log_tile_ver_linear, *log_dec, *log_frame_psnr;
	FILE* fVPSNR;
	// generate log files
	sprintf(fname, "result/video_%s/log_tile_ver_simple_trace_BW_%d_METHOD_%d_VP_EST_METHOD_%d_INTER_%d_BUFF_%d_%df_%dx%d.txt",adaptLogic.metadata.video.name.c_str(), BW, adaptLogic.metadata.adapt.TILE_SELECT_METHOD, adaptLogic.metadata.adapt.VP_EST_METHOD,adaptLogic.metadata.video.INTERVAL, adaptLogic.metadata.video.BUFF, adaptLogic.metadata.vp.No_face, adaptLogic.metadata.vp.No_tile_h, adaptLogic.metadata.vp.No_tile_v);
	printf("#%s\n", fname);
	log_tile_ver = fopen(fname, "w");
	sprintf(fname, "result/video_%s/log_tile_ver_linear_simple_trace_BW_%d_METHOD_%d_VP_EST_METHOD_%d_INTER_%d_BUFF_%d_%df_%dx%d.txt",adaptLogic.metadata.video.name.c_str(), BW, adaptLogic.metadata.adapt.TILE_SELECT_METHOD, adaptLogic.metadata.adapt.VP_EST_METHOD,adaptLogic.metadata.video.INTERVAL, adaptLogic.metadata.video.BUFF, adaptLogic.metadata.vp.No_face, adaptLogic.metadata.vp.No_tile_h, adaptLogic.metadata.vp.No_tile_v);
	log_tile_ver_linear = fopen(fname, "w");
	printf("#%s\n", fname);
	sprintf(fname, "result/video_%s/log_dec_simple_trace_BW_%d_METHOD_%d_VP_EST_METHOD_%d_INTER_%d_BUFF_%d_%df_%dx%d.txt",adaptLogic.metadata.video.name.c_str(), BW, adaptLogic.metadata.adapt.TILE_SELECT_METHOD, adaptLogic.metadata.adapt.VP_EST_METHOD,adaptLogic.metadata.video.INTERVAL, adaptLogic.metadata.video.BUFF, adaptLogic.metadata.vp.No_face, adaptLogic.metadata.vp.No_tile_h, adaptLogic.metadata.vp.No_tile_v);
	log_dec = fopen(fname, "w");
	sprintf(fname, "result/video_%s/log_frame_simple_trace_BW_%d_METHOD_%d_VP_EST_METHOD_%d_INTER_%d_BUFF_%d_%df_%dx%d.txt",adaptLogic.metadata.video.name.c_str(), BW, adaptLogic.metadata.adapt.TILE_SELECT_METHOD,adaptLogic.metadata.adapt.VP_EST_METHOD, adaptLogic.metadata.video.INTERVAL, adaptLogic.metadata.video.BUFF, adaptLogic.metadata.vp.No_face, adaptLogic.metadata.vp.No_tile_h, adaptLogic.metadata.vp.No_tile_v);
	printf("#%s\n", fname);
	log_frame_psnr = fopen(fname, "w");
	//
	if(log_tile_ver == NULL || log_tile_ver_linear == NULL || log_dec == NULL){
	printf("#[init] Cannot open log files\n");
	exit(1);
	}
	if(adaptLogic.metadata.adapt.TILE_SELECT_METHOD > 1){
	fprintf(log_dec, "id\testThrp(kbps)\ttext_width\tbitrate(kbps)\n");
	fprintf(log_frame_psnr, "fid\tdecid\test_vp_psnr\tphi\ttheta\test_phi\test_theta\terr_phi\terr_theta\text_width\n");
	for(int ii=0; ii < N; ii++){
	fprintf(log_tile_ver, "\nseg #%d calcTime: %d(ms)\n", ii, adaptLogic.proc_time[ii]);
	for(int i=0; i < adaptLogic.metadata.vp.No_tile_v; i++){
	// versiona
	for(int k=0; k < adaptLogic.metadata.vp.No_face; k++){
	for(int j=0; j < adaptLogic.metadata.vp.No_tile_h; j++){
	fprintf(log_tile_ver, "%d ",adaptLogic.tile_ver[ii][k*adaptLogic.metadata.vp.No_tile_h*adaptLogic.metadata.vp.No_tile_v + i * adaptLogic.metadata.vp.No_tile_h + j]);
	fprintf(log_tile_ver_linear, "%d ",adaptLogic.tile_ver[ii][k*adaptLogic.metadata.vp.No_tile_h*adaptLogic.metadata.vp.No_tile_v + i * adaptLogic.metadata.vp.No_tile_h + j]);
	}
	fprintf(log_tile_ver, "\t");    
	}
	fprintf(log_tile_ver, "\n");
	}
	fprintf(log_tile_ver, "\n");
	// psnr
	for(int i=0; i < adaptLogic.metadata.vp.No_tile_v; i++){
	for(int k=0; k < adaptLogic.metadata.vp.No_face; k++){
	for(int j=0; j < adaptLogic.metadata.vp.No_tile_h; j++){
	fprintf(log_tile_ver, "%.2f ",adaptLogic.metadata.video.TILE_PSNR[ii][k*adaptLogic.metadata.vp.No_tile_h*adaptLogic.metadata.vp.No_tile_v + i * adaptLogic.metadata.vp.No_tile_h + j][adaptLogic.tile_ver[ii][k*adaptLogic.metadata.vp.No_tile_h*adaptLogic.metadata.vp.No_tile_v + i * adaptLogic.metadata.vp.No_tile_h + j]]);
	}
	fprintf(log_tile_ver, "\t");    
	}
	fprintf(log_tile_ver, "\n");
	}	  
	fprintf(log_tile_ver, "\n");
	// bitrate
	for(int i=0; i < adaptLogic.metadata.vp.No_tile_v; i++){
	for(int k=0; k < adaptLogic.metadata.vp.No_face; k++){
	for(int j=0; j < adaptLogic.metadata.vp.No_tile_h; j++){
	fprintf(log_tile_ver, "%.2f ",adaptLogic.metadata.video.TILE_BR[ii][k*adaptLogic.metadata.vp.No_tile_h*adaptLogic.metadata.vp.No_tile_v + i * adaptLogic.metadata.vp.No_tile_h + j][adaptLogic.tile_ver[ii][k*adaptLogic.metadata.vp.No_tile_h*adaptLogic.metadata.vp.No_tile_v + i * adaptLogic.metadata.vp.No_tile_h + j]]);
	}
	fprintf(log_tile_ver, "\t");    
	}	  
	fprintf(log_tile_ver, "\n");
	}
	fprintf(log_tile_ver, "\n");
	fprintf(log_tile_ver_linear, "\n");
	fprintf(log_dec, "%d\t%.2f\t%.2f\n", ii, adaptLogic.thrp.est_seg_thrp[ii], adaptLogic.proc_time[ii]/1000.0/1000.0);
	// calculate viewport psnr of frames in this interval
	// decide the adaptation result of this segment
	int *tmpVisi;
	// printf("#[write_result] seg: %d\n", ii);
	// printf("est_vp: (%d, %d)\n", adaptLogic->est_vp[ii][0], adaptLogic->est_vp[ii][1]);
	// tmpVisi = adaptLogic->getVisibleTile(adaptLogic->est_vp[ii]);
	// for(i=0; i < rows; i++){
	//   for(j=0; j < cols; j++){
	//   	printf("%d ", tmpVisi[i*cols + j]);
	//   }
	//  printf("\n");
	// }
	for(int k=0; k < adaptLogic.metadata.video.INTERVAL; k++){
	fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\n",ii * adaptLogic.metadata.video.INTERVAL + k, ii, adaptLogic.vpsnr[ii*adaptLogic.metadata.video.INTERVAL +k], adaptLogic.metadata.trace.frame_vp[0][ii*adaptLogic.metadata.video.INTERVAL +k][0], adaptLogic.metadata.trace.frame_vp[0][ii*adaptLogic.metadata.video.INTERVAL +k][1], adaptLogic.vp2.est_frame_vp[ii*adaptLogic.metadata.video.INTERVAL + k][0], adaptLogic.vp2.est_frame_vp[ii*adaptLogic.metadata.video.INTERVAL + k][1], adaptLogic.vp2.est_err[ii][0], adaptLogic.vp2.est_err[ii][1]);
	}

	fflush(log_tile_ver);
	fflush(log_tile_ver_linear);
	fflush(log_dec);
	fflush(log_frame_psnr);
	// if(ii==5)
	// 	exit(1);
	}

	}else{
	fprintf(log_frame_psnr, "fid\tdecid\tver\tbitrate\test_vp_psnr\n");
	for(int ii=0; ii < N; ii++){
	for(int k=0; k < adaptLogic.metadata.video.INTERVAL; k++){
	fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%.2f\n", ii * adaptLogic.metadata.video.INTERVAL + k, ii, adaptLogic.vpsnr[ii*adaptLogic.metadata.video.INTERVAL +k], adaptLogic.tile_ver[ii][0], adaptLogic.metadata.video.DASH_BR[ii][adaptLogic.tile_ver[ii][0]]);
	}
	}
	fflush(log_frame_psnr);
	fclose(log_frame_psnr);
	}
	printf("#[print_result] BW:%d METHOD:%d\t%.2f\n", BW, adaptLogic.metadata.adapt.TILE_SELECT_METHOD, avg_vp_psnr);
}
run_cfg load_run_cfg(const char* cfg){
	printf("#[load_run_cfg]:\n");
	string s;
	string delimiter = "=";
	string comment = "#";
	string key;
	string val_str;
	//
	string deli = ",";
	string token;
	size_t pos = 0;
	std::vector<int> method;
	std::vector<int> head_trace;
	std::vector<int> bw;
	//
	double val;
	size_t pos_deli = 0;
	size_t pos_comm;
	ifstream infile(cfg);
	run_cfg config;
	if(infile == NULL){
		cout << "Cannot open file " << cfg << endl;
	}else{
		cout << "Reading file " << cfg << endl;
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
		if(key.compare("TILE_SELECT_METHOD")==0){
			// printf("%sa\n", val_str.c_str());
			config.no_method = 0;
			while((pos=val_str.find(deli)) != std::string::npos){
				token = val_str.substr(0, pos);
				method.push_back((int) std::stod(token));
				// cout << std::stod(token) << endl;
				config.no_method ++;
				val_str.erase(0, pos + deli.length());
			}
		}
		if(key.compare("HEAD_TRACE")==0){
			// printf("%sa\n", val_str.c_str());
			config.no_trace = 0;
			while((pos=val_str.find(deli)) != std::string::npos){
				token = val_str.substr(0, pos);
				head_trace.push_back((int) std::stod(token));
				// cout << std::stod(token) << endl;
				config.no_trace++;
				val_str.erase(0, pos + deli.length());
			}
		}
		if(key.compare("BANDWIDTH")==0){
			// printf("%sa\n", val_str.c_str());
			config.no_bw = 0;
			while((pos=val_str.find(deli)) != std::string::npos){
				token = val_str.substr(0, pos);
				bw.push_back((int) std::stod(token));
				// cout << std::stod(token) << endl;
				config.no_bw++;
				val_str.erase(0, pos + deli.length());
			}
		}
	}
	config.method = method;
	config.trace = head_trace;
	config.bw = bw;
	return config;
}