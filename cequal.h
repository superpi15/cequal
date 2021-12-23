/*
 * Cequal: C-based equalizer for educational purpose 
 */

#ifndef CEQUAL_H
#define CEQUAL_H

#include <string>
#include <iostream>
#include <vector>
#include <assert.h>
#include <algorithm>
#include "math.h"

#include <thread>
#include <mutex>

#include <ncurses.h>

#include "cequal_al.h"

namespace ceq {
#include "wave.h"
#define CEQUAL_PI 3.1415926 

class Cpx_t {
public:
	typedef float data_t;
	data_t real, img;
	Cpx_t(data_t _real, data_t _img): real(_real), img(_img) {}
	Cpx_t(): real(0), img(0){}
	void print(){printf("(%f,%f)",real,img);}

};
inline Cpx_t operator*(const Cpx_t& n1, const Cpx_t& n2){
	return Cpx_t( n1.real * n2.real - n1.img * n2.img, n1.real * n2.img + n1.img * n2.real );
}
inline Cpx_t operator+(const Cpx_t& n1, const Cpx_t& n2){
	return Cpx_t( n1.real + n2.real, n1.img + n2.img );
}
inline Cpx_t operator-(const Cpx_t& n1, const Cpx_t& n2){
	return Cpx_t( n1.real - n2.real, n1.img - n2.img );
}

class Ceq_Man_t {
	bool al_play;
	bool nl_play;
	std::mutex vol_mutex;
	std::mutex vol_panel_mutex;
	int cpx_mult;
public:
	std::string filename;
	WavHeader_t header;
	int nBand;
	Ceq_Man_t(char * nfilename){
		cpx_mult = 0;
		setDefault();
		filename = nfilename;
		read_wav_header(nfilename,header,0);
	}
	~Ceq_Man_t(){
	}
	Oal_Man_t oal;
	void setDefault();
	void dumpSample();
	void scanSample();
	void writeHeader( std::ostream& ostr );

	double time_stamp;                // keep progress for dynamic adjustment 
	std::vector<Cpx_t> vDFTW;         // DFT weight
	std::vector<int> vBandFreq, vBandIndex;
	std::vector<double> vBandVol, vBandGain, vBandGainMul;
	std::vector<bool> vBandMute;
	std::vector<int> vFreq2Band;
	double totalVolumn, totalVolumnMul;
	void config_al_play(bool flag){ al_play = flag;}
	void config_nl_play(bool flag){ nl_play = flag;}

	void runFilter(std::ostream * postr=NULL, int fPlot=0);
	void run_prefetch_filter(std::ostream * postr=NULL, int fPlot=0, const bool fft=true);

private:
	void run_fft(std::vector<Cpx_t>& vCpxIn, std::vector<Cpx_t>& vResult, const std::vector<Cpx_t>& vDataIn, std::vector<Cpx_t> * pvCpx=NULL, bool fifft=false);
	int band_vol_buf_top;
	std::vector<std::vector<double> > vBandVolBuf;
	void compute_band_vol( std::vector<std::vector<Cpx_t> > * pvFreqCh, int vol_buf_index, int band_step, int min_freq, int size_of_window );

	void vol_update_panel( int vol_buf_index, int size_of_window );

	int base_band_power;
	int volumn_scale;
	int nl_band_focus;
};


// fllr: https://stackoverflow.com/questions/63614171/avaudiorecorder-generates-strange-wav-filewrong-header
inline void Ceq_Man_t::scanSample(){
	FILE * fptr = fopen(filename.c_str(),"r");
	if( !fptr ){
		printf("Cannot open \'%s\'\n", filename.c_str());
		return;
	}
	if( fseek(fptr, 44, SEEK_SET) ){
		printf("Invalid header of \'%s\'\n", filename.c_str());
		return;
	}
	long num_samples = (8 * header.data_size) / (header.channels * header.bits_per_channel);
	long bytes_per_sample = (header.channels * header.bits_per_channel) / 8;
	long bytes_per_channel = (bytes_per_sample / header.channels);
	char * data_buffer = (char*) malloc(bytes_per_sample);

	//if(postr && !fPlot) writeHeader(*postr);

	int nRead = 1;
	while(fread(data_buffer,4,1,fptr)){
		printf("%3d %c%c%c%c \n", nRead, data_buffer[0], data_buffer[1], data_buffer[2], data_buffer[3]);
		if( data_buffer[0]=='d' && data_buffer[1]=='a' && data_buffer[2]=='t' && data_buffer[3]=='a' ){
			fread(data_buffer,4,1,fptr);
			int data_size = *((int*)data_buffer);
			printf("size = %d bytes\n", data_size);
			fseek(fptr, data_size, SEEK_CUR);
		}
		nRead += 4;
		if(nRead >= header.data_size+50) break;
	}
//	for(int i = 0; i < num_samples; i ++){
//		if(!fread(data_buffer,bytes_per_sample,1,fptr)){
//			printf("Data loading mismatch\n");
//			goto FINALIZE;
//		}
//		printf("%c%c%c%c\n", data_buffer[0], data_buffer[1], data_buffer[2], data_buffer[3]);
//		for(int j = 0; j < header.channels; j ++){
//			int channel_data_in = 0;
//			for(int k = 0; k < bytes_per_channel; k ++ )
//				channel_data_in |= (data_buffer[ j*bytes_per_channel + k ]&(k==bytes_per_channel-1? ~0: 0x00ff)) << (8*k);
//			if( 1 == bytes_per_channel )
//				channel_data_in -= 128;
//		}
//	}
FINALIZE:
	free(data_buffer);
	fclose(fptr);
}

inline void Ceq_Man_t::compute_band_vol( std::vector<std::vector<Cpx_t> > * pvFreqCh, int vol_buf_index, int band_step, int min_freq, int size_of_window ){
	while(!vol_mutex.try_lock())
		;
	double decay = 0.3;
	std::vector<std::vector<Cpx_t> >& vFreqCh = *pvFreqCh;
	int num_ff_freq = vFreqCh[0].size();
	std::vector<double> vBandVolTmp( vBandVol.size(), 0 );

	for(int ch = 0; ch < vFreqCh.size(); ch ++){
		int bid = 0;
		double accu_vol = 0;
		int ffband_num = 0;
		std::vector<Cpx_t>& vFreq = vFreqCh[ch];
		for(int i = 0; i < vFreq.size(); i ++, ffband_num ++ ){
			if( bid != vFreq2Band[i] ){
				if( ffband_num )
					vBandVolTmp[bid] += accu_vol / ffband_num;
				bid = vFreq2Band[i];
				accu_vol = 0;
				ffband_num = 0;
			}
			//assert(vFreq[i].real<(1<<16));
			accu_vol += (double) 
				//vFreq[i].real / size_of_window
				sqrt((vFreq[i] * Cpx_t( vFreq[i].real, - vFreq[i].img )).real)/size_of_window
				;
		}
	}
//
	for(int i = 0; i < vBandVol.size(); i ++)
		vBandVol[i] = (double) decay * vBandVol[i] + (1.0f-decay) * vBandVolTmp[i];

	vBandVolBuf[vol_buf_index] = vBandVol;
//
//
//	{
//		//int max_vol = 1<<(header.bits_per_channel-1);
//
//		for(int i = 0; i < vBandVol.size(); i ++){
			//
//			int j;
//			double vol_ratio = vBandVol[i] / size_of_window;
//			int vol_portion = volumn_scale * vol_ratio;
//			mvprintw(6+i,1,"%6.2f ", vol_ratio);
//			for(j = 0; j < volumn_scale; j ++)
//				mvprintw(6+i,8+j,"%c", j<vol_portion? '|': ' ');
//		}
//	}

	vol_mutex.unlock();
}

inline void Ceq_Man_t::vol_update_panel( int vol_buf_index, int size_of_window ){
	int dy = 7, dx = 6;
	while(!vol_panel_mutex.try_lock())
		;

	for(int i = 0; i < nBand; i ++){
		int j;
		double vol_ratio = vBandVolBuf[vol_buf_index][i] / size_of_window;
		int vol_portion = volumn_scale * vol_ratio;

		mvprintw(dy+nBand-i-1,dx,"%c", i == nl_band_focus? '*': ' ');
		mvprintw(dy+nBand-i-1,dx+1," %c %6d %6.2f [", vBandMute[i]? 'm': ' ', vBandFreq[i], vol_ratio);
		for(j = 0; j < volumn_scale; j ++)
			mvprintw(dy+nBand-i-1,dx+20+j,"%c", j<vol_portion? '|': ' ');

		mvprintw(dy+nBand-i-1,dx+30+volumn_scale,"] x %7.4f x %7.4f ", vBandGain[i], totalVolumn);
	}

	vol_panel_mutex.unlock();
}

inline void Ceq_Man_t::setDefault(){
	al_play = true;
	nl_play = true;
	cpx_mult = 1;
	vDFTW.resize(20);
	for(int i = 0; i < vDFTW.size(); i ++){
		double angular_unit = 2*CEQUAL_PI / (1<<i);
		vDFTW[i] = Cpx_t(cpx_mult * cos(angular_unit),cpx_mult * sin(angular_unit));
		//printf(" %d %d\n", vDFTW[i].real, vDFTW[i].img);
	}
	nBand = 10;
	volumn_scale = 20;
	base_band_power = 6;
	vBandFreq.resize(nBand);
	vBandGain.resize(nBand,0);
	vBandGainMul.resize(nBand,1.0f);
	totalVolumn = 0;
	totalVolumnMul = 1.0f;
	vBandVol.resize(nBand,0);
	vBandMute.resize(nBand,0);
	for(int i = 0; i < nBand; i ++){
		vBandFreq[i] = 1<<(i+base_band_power);
	}
	vBandFreq.back() = 20000;
}

inline void Ceq_Man_t::run_fft(std::vector<Cpx_t>& vCpxIn, std::vector<Cpx_t>& vResult, const std::vector<Cpx_t>& vDataIn, std::vector<Cpx_t> * pvCpx, bool fifft){
	//std::vector<Cpx_t> vCpxIn, vResult;
	int size_of_window = vDataIn.size();
	
//	for(int i = size_of_window; i > 0 ; i >>=1 )
//		assert( !(i&1) || 1==i );

	int half_window_size = size_of_window/2;
	if( vCpxIn.size() < size_of_window )
		vCpxIn .resize(size_of_window);
	if( vResult.size() < size_of_window )
		vResult.resize(size_of_window);
	for(int i = 0; i < half_window_size; i +=2){
		vCpxIn[ i ]    = vDataIn[ i ];
		vCpxIn[ i + 1] = vDataIn[ i + half_window_size ];
		vCpxIn[ i + half_window_size ]     = vDataIn[i + 1 ];
		vCpxIn[ i + half_window_size + 1 ] = vDataIn[i + 1 + half_window_size ];
	}


	int ffstep, fflevel;
	for(ffstep = 2, fflevel = 1; ffstep <= size_of_window; ffstep <<= 1, fflevel ++ ){
		//printf("%4d %4d %7.4f %7.4f \n", ffstep, fflevel, vDFTW[fflevel].real, vDFTW[fflevel].img);
		//#pragma omp parallel for num_threads (2) if(256<size_of_window)
		//#pragma omp parallel for num_threads (2)
		for(int i = 0; i < size_of_window; i += ffstep ){ // iterate groups of nodes on current layer 
			Cpx_t w(1,0);
			// use symmetry property 
			for(int j = 0; j < (ffstep>>1); j ++){             // iterate nodes in a window  
				int idx0, idx1;
				if( j < (ffstep>>1) )
					idx0 = i+j, idx1 = i+j+(ffstep>>1);
				else
					idx1 = i+j, idx0 = i+j-(ffstep>>1);
				Cpx_t mult = vCpxIn[idx1] * w;
				//vResult[i+j] = vCpxIn[idx0] + vCpxIn[idx1] * w;
				vResult[i+j] = vCpxIn[idx0] + mult;
				vResult[i+j+(ffstep>>1)] = vCpxIn[idx0] - mult;
				w = w * (fifft? Cpx_t(vDFTW[fflevel].real, -vDFTW[fflevel].img) :vDFTW[fflevel]);
				//w = vDFTW[fflevel];
			}
		}
//		if(ffstep!=8) continue;
//		for(int i = 0; i < size_of_window; i += ffstep ){ // iterate groups of nodes on current layer 
//			for(int j = 0; j < (ffstep>>1); j ++){             // iterate nodes in a window  
//				if( (int)vResult[i+j].real == (int)vResult[i+j+(ffstep>>1)].real )
//					continue;
////				else
////				if( vResult[i+j].real == 0 && vResult[i+j].img == 0 ){
////					continue;
////				}
////				if( vResult[i+j].real == vResult[i+j+(ffstep>>1)-1].real && j < (ffstep>>1) ){
////					assert(false);
////					break;
////				}
//				for(int k = 0; k < ffstep; k ++)
//					printf("%8d %3d: %6.2f %6.2f\n",i , k, vResult[i+k].real, vResult[i+k].img);
				//
//				printf("%8d %3d: %6d vs. %6d\n", i, j, (int)vResult[i+j].real, (int)vResult[i+j+(ffstep>>1)].real );
//				printf("%8d %3d: %f vs. %f\n", i, j, vResult[i+j].real, vResult[i+j+(ffstep>>1)].real );
//				//std::cout << vResult[i+j].real <<" "<< vResult[i+j+(ffstep>>1)].real << std::endl;
//				printf("\n");
//				break;
//			}
//		}
		vCpxIn.swap(vResult);
	//	if(fflevel==2)\
		break;
	}
	//exit(0);

	if(pvCpx){
		if(!fifft)
			for(int i = 0; i < vCpxIn.size(); i ++){
				vCpxIn[i].real *= vBandMute[ vFreq2Band[i] ]? 0: vBandGainMul[ vFreq2Band[i] ];
				vCpxIn[i].img  *= vBandMute[ vFreq2Band[i] ]? 0: vBandGainMul[ vFreq2Band[i] ];
			}
		pvCpx->swap(vCpxIn);
	}
}


inline void Ceq_Man_t::run_prefetch_filter(std::ostream * postr, int fPlot, const bool fft){
	std::string base_fname = basename(filename.c_str());
	FILE * fptr = fopen(filename.c_str(),"r");
	if( !fptr ){
		printf("Cannot open \'%s\'\n", filename.c_str());
		return;
	}
	if( fseek(fptr, 44, SEEK_SET) ){
		printf("Invalid header of \'%s\'\n", filename.c_str());
		return;
	}
	if(postr && !fPlot) writeHeader(*postr);

	long bytes_per_sample = (header.channels * header.bits_per_channel) / 8;
	char * data_buffer = (char*) malloc(bytes_per_sample);
	long bytes_per_channel = (bytes_per_sample / header.channels);
	int data_size = header.data_size;
	if( header.data_chunk_header[0]=='L' ){
		int nRead = 1;
		while(fread(data_buffer,4,1,fptr)){
			//printf("%3d %c%c%c%c \n", nRead, data_buffer[0], data_buffer[1], data_buffer[2], data_buffer[3]);
			if( data_buffer[0]=='d' && data_buffer[1]=='a' && data_buffer[2]=='t' && data_buffer[3]=='a' ){
				if(!fread(data_buffer,4,1,fptr)){
					printf("Invalid header of \'%s\'\n", filename.c_str());
					free(data_buffer);
					return;
				}
				data_size = *((int*)data_buffer);
				//printf("size = %d bytes\n", data_size);
				//fseek(fptr, data_size, SEEK_CUR);
				break;
			}
			//nRead += 4;
			//if(nRead >= header.data_size+50) break;
		}
		//printf("%c%c%c%c \n", data_buffer[0], data_buffer[1], data_buffer[2], data_buffer[3]);
		//data_size = *((int*) data_buffer);
		//printf("get data size = %ld\n", data_size);
	}
	long num_samples = (8 * data_size) / (header.channels * header.bits_per_channel);
	
	std::vector< std::vector<Cpx_t::data_t> > vData, vDataOut;


	int al_top = 0, al_buf_size = 10;
	void * al_buf = NULL;
	if( this->al_play ){
		al_buf = malloc( al_buf_size * bytes_per_sample );
	}

	long ubound, lbound;
	if( 8 == header.bits_per_channel ){
		ubound = 128, lbound = -127;
	} else
	if( 16 == header.bits_per_channel ){
		ubound = 32768, lbound = -32767;
	} else
	if( 32 == header.bits_per_channel ){
		ubound = 2147483648;
		lbound = -2147483647;
	} else
		assert(false);

	int size_of_window = 512;
	int nl_update_per_sample = 2048;
	int al_sample_per_batch  = header.sample_rate/24; // al ext buffer size
	int max_freq = header.sample_rate/2; // The smallest wave crest and through are matched a input sample
	int half_window_size = size_of_window/2;
	int min_freq = max_freq / half_window_size;
	double angular_unit = 2 * CEQUAL_PI / size_of_window;
	int band_step = (max_freq-min_freq) / half_window_size;

	std::vector< std::vector<Cpx_t> > vCpxInCh(header.channels), vResultCh(header.channels);
	std::vector< std::vector<Cpx_t> > vDataTmpCh(header.channels);
	std::vector< std::vector<Cpx_t> > vFreqTmpCh(header.channels), vFreqTmpChTh;
	vBandVolBuf.clear();
	vBandVolBuf.resize( 1 + (num_samples+size_of_window) / nl_update_per_sample );
	for(int j = 0; j < header.channels; j ++)
		vDataTmpCh[j].resize(size_of_window);

	printf("num samples %6ld\n", num_samples);
	printf("max freq    %6d\n", max_freq);
	printf("window size %6d\n", size_of_window);
	printf("min freq    %6d\n", min_freq);
	printf("step        %6d\n", band_step);

	vFreq2Band.resize(size_of_window,0);
	vBandIndex.resize(nBand,0);
	{
		int bid = 0;
		for(int i = 0; i < half_window_size; i ++){
			if( bid < vBandFreq.size()-1 && vBandFreq[bid] < (int)min_freq + band_step * i ){
				vBandIndex[bid] = i;
				bid ++ ;
			}
			vFreq2Band[i] = bid;
			//vFreq2Band[i+half_window_size] = bid;
			vFreq2Band[size_of_window-i-1] = vBandFreq.size()-1;
			//printf("%d %d %d %d\n",bid, nBand, (int)min_freq + band_step * i, vBandFreq[bid]);\
			assert(bid<nBand);
		}
//		for(int i = 0; i < half_window_size; i ++){
//			printf("%7d %7d %7d :  %7d %7d\n",bid, nBand, (int)min_freq + band_step * i, vBandFreq[vFreq2Band[i]]
//				//, vBandFreq[vFreq2Band[i+half_window_size]]
//				, vBandFreq[vFreq2Band[size_of_window-i-1]]
//				);
//		}
	}
	int nl_char = -1;
	if( this->nl_play ){
		// http://www.cplusplus.com/forum/beginner/248599/
		initscr();
		start_color();
		clear();
		curs_set(0);
  		noecho();
  		cbreak();
		keypad(stdscr, true);
		refresh();
		nl_band_focus = -1;
		int dy = 20;
		init_pair(3, COLOR_BLACK, COLOR_WHITE);
		attron(COLOR_PAIR(3));
		mvprintw(dy+0,41,"clip name   %35s", base_fname.c_str());
		mvprintw(dy+1,41,"num samples %35ld", num_samples);
		mvprintw(dy+2,41,"max freq    %35d", max_freq);
		mvprintw(dy+3,41,"window size %35d", size_of_window);
		mvprintw(dy+4,41,"min freq    %35d", min_freq);
		mvprintw(dy+5,41,"step        %35d", band_step);
		attroff(COLOR_PAIR(3));
		mvprintw(29+1,2,"Cequal - Terminal DJ");
		mvprintw(29+2,2,"Author: Nero Horus");
		mvprintw(29+3,2,"Source:     %35s", "https://github.com/superpi15/cequal");


//		mvprintw(0,22,"  _____________________  _______  __________  _____   ___     ______      ___");
//		mvprintw(1,22," /       /      /  __  \/       \/__/  ___  \/  _  \ /  /    /  _   \    /__/");
//		mvprintw(2,22," ===  ==/  ====/  /_/__/  /  /  /  /  /  /  /  /_\  /  /    /  / /  /   /  / ");
//		mvprintw(3,22,"  /  / /  ____/  __  \/  /  /  /  /  /  /  /  __   /  /__  /  /_/  /___/  /  ");
//		mvprintw(4,22," /__/ /______/__/ /__/__/__/__/__/__/  /__/__/ /__/_____/ /_______/______/   ");
//		mvprintw(5,22,"-----------------------------------------------------------------------------");


		init_pair(2, COLOR_CYAN, COLOR_BLACK);
		attron(COLOR_PAIR(2));
		mvprintw(0,4,"  _____________________  _______  __________  _____   ___     ______      ___");
		mvprintw(1,4," /       /      /  __  \\/       \\/__/  ___  \\/  _  \\ /  /    /  _   \\    /__/");
		mvprintw(2,4," ===  ==/  ====/  /_/__/  /  /  /  /  /  /  /  /_\\  /  /    /  / /  /   /  / ");
		mvprintw(3,4,"  /  / /  ____/  __  \\/  /  /  /  /  /  /  /  __   /  /__  /  /_/  /___/  /  ");
		mvprintw(4,4," /__/ /______/__/ /__/__/__/__/__/__/  /__/__/ /__/_____/ /_______/______/   ");
		mvprintw(5,4,"-----------------------------------------------------------------------------");
		attroff(COLOR_PAIR(2));

		init_pair(1, COLOR_BLACK, COLOR_YELLOW);
		attron(COLOR_PAIR(1));
		mvprintw(20,1,"%6s %-30s","right","switch to overall (default)");
		mvprintw(21,1,"%6s %-30s","left","switch to bands");
		mvprintw(22,1,"%6s %-30s","up","previous band");
		mvprintw(23,1,"%6s %-30s","down","next band");
		mvprintw(24,1,"%6s %-30s","-","reduce band voloumn");
		mvprintw(25,1,"%6s %-30s","+","increase band voloumn");
		mvprintw(26,1,"%6s %-30s","2","double band voloumn");
		mvprintw(27,1,"%6s %-30s","/","half band voloumn");
		mvprintw(27,1,"%6s %-30s","m","mute band");
		mvprintw(28,1,"%6s %-30s","q","exit");
		attroff(COLOR_PAIR(1));
		refresh();
	}



	vData.resize( header.channels );
	vDataOut.resize( header.channels );
	for(int i = 0; i < header.channels; i ++){
		vData   [i].resize( num_samples + size_of_window, 0.0f );
		vDataOut[i].resize( num_samples, 0.0f );
	}

	std::vector<std::thread> vTh, vVolPanelTh;
	std::vector<Cpx_t::data_t> vInvSin, vInvCos;
	vInvSin.resize(size_of_window);
	vInvCos.resize(size_of_window);
	for(int i = 0; i < size_of_window; i ++){
		vInvSin[i] = cpx_mult * sin( - angular_unit * i * half_window_size );
		vInvCos[i] = cpx_mult * cos( - angular_unit * i * half_window_size );
		//std::cout << vInvSin[i] <<" "<< vInvCos[i] << "\n";
	}


	for(int i = 0; i < num_samples; i ++){
		if(!fread(data_buffer,bytes_per_sample,1,fptr)){
			printf("Data loading mismatch\n");
			goto FINALIZE;
		}
		for(int j = 0; j < header.channels; j ++){
			int channel_data_in = 0;
			for(int k = 0; k < bytes_per_channel; k ++ )
				channel_data_in |= (data_buffer[ j*bytes_per_channel + k ]&(k==bytes_per_channel-1? ~0: 0x00ff)) << (8*k);
			if( 1 == bytes_per_channel )
				channel_data_in -= 128;

			vData[j][i+half_window_size] = channel_data_in;
		}
	}




	if( this->al_play ){
		ALenum al_format;
		if( header.channels == 1 ){
			if( 8 == header.bits_per_channel )
				al_format = AL_FORMAT_MONO8;
			else
			if( 16 == header.bits_per_channel )
				al_format = AL_FORMAT_MONO16;
			else {
				assert(false);
			}
		} else
		if( header.channels == 2 ){
			if( 8 == header.bits_per_channel )
				al_format = AL_FORMAT_STEREO8;
			else
			if( 16 == header.bits_per_channel )
				al_format = AL_FORMAT_STEREO16;
			else {
				assert(false);
			}
		} else
			assert(false);

		oal.init();
		oal.init_extern_queue( al_sample_per_batch, al_format, header.sample_rate );
		oal.wait();
	}

	for(int i = 0; i < num_samples; i ++){
		if( this->nl_play ){
			timeout(0);
			nl_char = getch();
		}
		#pragma omp parallel for num_threads (2) if(1<header.channels)
		for(int j = 0; j < header.channels; j ++){
			std::vector<Cpx_t>& vDataTmp = vDataTmpCh[j];
			for(int l = 0; l < size_of_window; l ++){
				vDataTmp[l] = Cpx_t(vData[j][l+i], 0);
			}
			std::vector<Cpx_t>& vCpx = vFreqTmpCh[j];
			run_fft(vCpxInCh[j], vResultCh[j], vDataTmp, &vCpx);
			//run_fft(vCpx, &vDataTmp, 1);

			Cpx_t::data_t real = 0;
			for(int k = 0; k < size_of_window; k ++) // k-th freq 
				real += vCpx[k].real * vInvCos[k] - vCpx[k].img * vInvSin[k];

			real *= totalVolumnMul;
			int channel_data_out = real/size_of_window/cpx_mult;
			if( ubound < channel_data_out )
				channel_data_out = ubound;
			else
			if( lbound > channel_data_out )
				channel_data_out = lbound;
			for(int k = 0; k < bytes_per_channel; k ++ )
				data_buffer[ j*bytes_per_channel + k ] = (channel_data_out >> (8*k)) & 0x00ff;
			vDataOut[j][i] = channel_data_out;
			//if(postr && fPlot) (*postr) << channel_data_out << " ";

		}
		for(int j = 0; j < header.channels; j ++){
			if( this->al_play ){
				for(int k = 0; k < bytes_per_channel; k ++ )
					oal.push_byte(data_buffer[ j*bytes_per_channel + k ]);
			}
		}
		if(postr && fPlot) (*postr) << "\n";
		//if(postr && !fPlot) postr->write( data_buffer, bytes_per_sample );

		if( this->nl_play && 0==(i%nl_update_per_sample) ){
			int time_scale = 30;
			int progress = time_scale * i/num_samples;
			int second = i / header.sample_rate;
			int minute = second/60;
			int dx = 16;
			int dy = 18;
			mvprintw(dy,0,"        mm:ss");
			mvprintw(dy,dx+2,"%2d:%2d", second/60, second-minute*60);
			mvprintw(dy,dx+8,"[ ");
			int j;
			for(j=0; j<progress; j++)
				mvprintw(dy,dx+10+j,"=");
			for(; j<time_scale; j++)
				mvprintw(dy,dx+10+j," ");
			mvprintw(dy,dx+10+time_scale,"]");
			if( vol_mutex.try_lock() )
			{
				vol_mutex.unlock();
				band_vol_buf_top = i/nl_update_per_sample;
				vBandVolBuf[band_vol_buf_top].clear();
				vFreqTmpChTh = vFreqTmpCh;
				vTh.push_back(std::thread( &Ceq_Man_t::compute_band_vol, this, &vFreqTmpChTh, band_vol_buf_top, band_step, min_freq, size_of_window));
			}
			int vol_panel_match_index = (oal.num_al_push()) * al_sample_per_batch / nl_update_per_sample;
			if( 0 <= vol_panel_match_index && vBandVolBuf[vol_panel_match_index].size() && vol_panel_mutex.try_lock() ){
				vol_panel_mutex.unlock();
				vVolPanelTh.push_back(std::thread( &Ceq_Man_t::vol_update_panel, this, vol_panel_match_index, size_of_window));
			}
			
			refresh();
		}


		if( -1 != nl_char ){
			if( 'q' == nl_char ){
				mvprintw(40,1,"received terminating sequence (flushing buffer) ");
				refresh();
				break;
			} else 
				;//mvprintw(5,1,"last input: %2c ", nl_char);
			switch(nl_char){
				case '-':
					if( -1 == nl_band_focus ){
						totalVolumn -= -10.0f<totalVolumn? 0.0125: 0;
						totalVolumnMul = pow(2,totalVolumn);
					}
					else {
						vBandGain[nl_band_focus] -= -10.0f<vBandGain[nl_band_focus]? 0.0125: 0;
						vBandGainMul[nl_band_focus] = pow(2,vBandGain[nl_band_focus]);
					}
					break;
				case '=':
					if( -1 == nl_band_focus ){
						totalVolumn += 5.0f>totalVolumn? 0.0125: 0;
						totalVolumnMul = pow(2,totalVolumn);
					}
					else {
						vBandGain[nl_band_focus] += 10.0f>vBandGain[nl_band_focus]? 0.0125: 0;
						vBandGainMul[nl_band_focus] = pow(2,vBandGain[nl_band_focus]);
					}
					
					break;
				case 'm':
					if( -1 == nl_band_focus )
						;//totalVolumnMul = 0;
					else 
						vBandMute[nl_band_focus] = !vBandMute[nl_band_focus];
					break;
				case '2':
					if( -1 == nl_band_focus ){
						totalVolumn += 5.0f>totalVolumn? 1: 0;
						totalVolumnMul = pow(2,totalVolumn);
					} else {
						vBandGain[nl_band_focus] += 10.0f>vBandGain[nl_band_focus]? 1: 0;
						vBandGainMul[nl_band_focus] = pow(2,vBandGain[nl_band_focus]);
					}
					break;
				case '/':
					if( -1 == nl_band_focus ){
						totalVolumn -= -10.0f<totalVolumn? 1: 0;
						totalVolumnMul = pow(2,totalVolumn);
					} else {
						vBandGain[nl_band_focus] -= -10.0f<vBandGain[nl_band_focus]? 1: 0;
						vBandGainMul[nl_band_focus] = pow(2,vBandGain[nl_band_focus]);
					}
					break;
				case KEY_UP:
					nl_band_focus += nBand-1 > nl_band_focus? 1: 0;
					break;
				case KEY_DOWN:
					nl_band_focus -= 0 < nl_band_focus? 1: 0;
					break;
				case KEY_LEFT:
					if(-1 == nl_band_focus)
						nl_band_focus = 0;
					break;
				case KEY_RIGHT:
					nl_band_focus = -1;
					break;
			}
			//clrtoeol();
			nl_char = -1;
		}

	}
	for(int i = 0; i < vTh.size(); i ++)
		if(vTh[i].joinable()) vTh[i].join();
	for(int i = 0; i < vVolPanelTh.size(); i ++)
		if(vVolPanelTh[i].joinable()) vVolPanelTh[i].join();
	oal.flush();
	oal.wait();

	//vDataOut = vData;
/**/
	if( postr ){
		if( !fPlot && header.data_chunk_header[0]=='L' ){
			postr->write( header.list_data, header.data_size );
			(*postr)<<"data"<<data_size;
		}
		for(int i = 0; i < num_samples; i ++){
			for(int j = 0; j < header.channels; j ++){
				int channel_data_out = vDataOut[j][i];
				for(int k = 0; k < bytes_per_channel; k ++ )
					data_buffer[ j*bytes_per_channel + k ] = (channel_data_out >> (8*k)) & 0x00ff;
				if(fPlot) (*postr) << channel_data_out << " ";
			}
			if( fPlot) (*postr) << "\n";
			if(!fPlot) postr->write( data_buffer, bytes_per_sample );
		}
	}
	/**/

FINALIZE:
	if( al_play ){
		oal.finalize();
		free(al_buf);
	}

	if( this->nl_play ){
		endwin();
	}
	free(data_buffer);
	fclose(fptr);
}

inline void Ceq_Man_t::writeHeader( std::ostream& ostr ){
	unsigned int i;
	int j = 0;
	for(i = 0; i < 4; i ++, j++) ostr << header.riff[i];
	for(i = 0; i < 4; i ++, j++) ostr << (unsigned char)( (header.overall_size >> i*8) & 0xff );
	//WAV
	for(i = 0; i < 4; i ++, j++) ostr << header.wave[i];
	//fmt
	for(i = 0; i < 4; i ++, j++) ostr << header.fmt_chunk_marker[i];
	//17-20
	for(i = 0; i < 4; i ++, j++) ostr << (unsigned char)( (header.length_of_fmt >> i*8) & 0xff );
	//21-22
	for(i = 0; i < 2; i ++, j++) ostr << (unsigned char)( (header.format_type >> i*8) & 0xff );
	//23-24
	for(i = 0; i < 2; i ++, j++) ostr << (unsigned char)( (header.channels >> i*8) & 0xff );
	//25-28
	for(i = 0; i < 4; i ++, j++) ostr << (unsigned char)( (header.sample_rate >> i*8) & 0xff );
	//29-32
	for(i = 0; i < 4; i ++, j++) ostr << (unsigned char)( (header.byterate >> i*8) & 0xff );
	//33-34
	for(i = 0; i < 2; i ++, j++) ostr << (unsigned char)( (header.block_align >> i*8) & 0xff );
	//35-36
	for(i = 0; i < 2; i ++, j++) ostr << (unsigned char)( (header.bits_per_channel >> i*8) & 0xff );
	//37-40
	for(i = 0; i < 4; i ++, j++) ostr << header.data_chunk_header[i];
	//41-44
	for(i = 0; i < 4; i ++, j++) ostr << (unsigned char)( (header.data_size >> i*8) & 0xff );
}
}
#endif