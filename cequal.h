/*
 * Cequal: C-based equalizer for educational purpose 
 */

#ifndef CEQUAL_H
#define CEQUAL_H

#include <string>
#include <iostream>
#include <vector>
#include <assert.h>
#include "math.h"


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
	void writeHeader( std::ostream& ostr );

	double time_stamp;                // keep progress for dynamic adjustment 
	std::vector<Cpx_t> vDFTW;         // DFT weight
	void config_play(bool flag){ al_play = flag;}

	void runFilter(std::ostream * postr=NULL, int fPlot=0);
	void run_prefetch_filter(std::ostream * postr=NULL, int fPlot=0, const bool fft=true);
	void run_fft(const std::vector<Cpx_t::data_t>& vDataIn, std::vector<Cpx_t> * pvCpx=NULL);
	int window_size;
};

inline void Ceq_Man_t::setDefault(){
	al_play = true;
	cpx_mult = 1;
	vDFTW.resize(20);
	for(int i = 0; i < vDFTW.size(); i ++){
		double angular_unit = 2*CEQUAL_PI / (1<<i);
		vDFTW[i] = Cpx_t(cpx_mult * cos(angular_unit),cpx_mult * sin(angular_unit));
		//printf(" %d %d\n", vDFTW[i].real, vDFTW[i].img);
	}
	nBand = 11;
	window_size = nBand;
	int max_freq = 44100;
}

inline void Ceq_Man_t::run_fft(const std::vector<Cpx_t::data_t>& vDataIn, std::vector<Cpx_t> * pvCpx){
	std::vector<Cpx_t> vCpxIn, vResult;
	int size_of_window = vDataIn.size();
	
	for(int i = size_of_window; i > 0 ; i >>=1 )
		assert( !(i&1) || 1==i );

	int half_window_size = size_of_window/2;
	vCpxIn .resize(size_of_window);
	vResult.resize(size_of_window);
	for(int i = 0; i < half_window_size; i +=2){
		vCpxIn[ i ]    = Cpx_t( vDataIn[ i ], 0 );
		vCpxIn[ i + 1] = Cpx_t( vDataIn[ i + half_window_size ], 0 );
		vCpxIn[ i + half_window_size ]     = Cpx_t( vDataIn[i + 1 ], 0 );
		vCpxIn[ i + half_window_size + 1 ] = Cpx_t( vDataIn[i + 1 + half_window_size ], 0 );
	}


	int ffstep, fflevel;
	for(ffstep = 2, fflevel = 1; ffstep <= size_of_window; ffstep <<= 1, fflevel ++ ){
		//printf("%4d %4d %7.4f %7.4f \n", ffstep, fflevel, vDFTW[fflevel].real, vDFTW[fflevel].img);
		for(int i = 0; i < size_of_window; i += ffstep ){ // iterate groups of nodes on current layer 
			Cpx_t w(1,0);
			for(int j = 0; j < ffstep; j ++){             // iterate nodes in a window  
				int idx0, idx1;
				if( j < (ffstep>>1) )
					idx0 = i+j, idx1 = i+j+(ffstep>>1);
				else
					idx1 = i+j, idx0 = i+j-(ffstep>>1);
				vResult[i+j] = vCpxIn[idx0] + vCpxIn[idx1] * w;
				w = w * vDFTW[fflevel];
				//w = vDFTW[fflevel];
			}
		}
		vCpxIn.swap(vResult);
	//	if(fflevel==2)\
		break;
	}
	//exit(0);

	if(pvCpx)
		pvCpx->swap(vCpxIn);
}


inline void Ceq_Man_t::run_prefetch_filter(std::ostream * postr, int fPlot, const bool fft){
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
	long bytes_per_channel = (bytes_per_sample / header.channels);
	long num_samples = (8 * header.data_size) / (header.channels * header.bits_per_channel);
	char * data_buffer = (char*) malloc(bytes_per_sample);
	std::vector< std::vector<Cpx_t::data_t> > vData, vDataOut;


	int al_top = 0, al_buf_size = 22050;
	void * al_buf = NULL;
	if( this->al_play ){
		al_buf = malloc( al_buf_size * bytes_per_sample );
	}


	int half_window_size = this->window_size/2;

	int max_freq = header.sample_rate / 2; // The smallest wave crest and through are matched a input sample
	int size_of_window = 512;
	double min_freq = (double) max_freq / (size_of_window / 2);
	double angular_unit = 2 * CEQUAL_PI / size_of_window;
	
	std::vector< std::vector<Cpx_t::data_t> > vDataTmpCh(header.channels);
	for(int j = 0; j < header.channels; j ++)
		vDataTmpCh[j].resize(size_of_window);

	printf("max freq    %6d\n", max_freq);
	printf("window size %6d\n", size_of_window);
	printf("min freq    %7.3f\n", min_freq);
	printf("step        %7.3f\n", (double) (max_freq-min_freq)/size_of_window);


	vData.resize( header.channels );
	vDataOut.resize( header.channels );
	for(int i = 0; i < header.channels; i ++){
		vData   [i].resize( num_samples + size_of_window, 0.0f );
		vDataOut[i].resize( num_samples, 0.0f );
	}

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
		oal.init_extern_queue( header.sample_rate/4, al_format, header.sample_rate );
		oal.play();
	}

	for(int i = 0; i < num_samples; i ++){

		#pragma omp parallel for num_threads (2) if(1<header.channels)
		for(int j = 0; j < header.channels; j ++){
			std::vector<Cpx_t::data_t>& vDataTmp = vDataTmpCh[j];
			for(int l = 0; l < size_of_window; l ++){
				vDataTmp[l] = vData[j][l+i];
			}
			std::vector<Cpx_t> vCpx;
			run_fft(vDataTmp, &vCpx);

			Cpx_t::data_t real = 0;
			for(int k = 0; k < size_of_window; k ++) // k-th freq 
				real += vCpx[k].real * vInvCos[k] - vCpx[k].img * vInvSin[k];


			int channel_data_out = real/size_of_window/cpx_mult;
			for(int k = 0; k < bytes_per_channel; k ++ )
				data_buffer[ j*bytes_per_channel + k ] = (channel_data_out >> (8*k)) & 0x00ff;
			//if(postr && fPlot) (*postr) << channel_data_out << " ";

		}
		for(int j = 0; j < header.channels; j ++){
			if( this->al_play ){
				for(int k = 0; k < bytes_per_channel; k ++ )
					oal.push_byte(data_buffer[ j*bytes_per_channel + k ]);
			}
		}
		if(postr && fPlot) (*postr) << "\n";
		if(postr && !fPlot) postr->write( data_buffer, bytes_per_sample );
	}
	oal.wait();

	//vDataOut = vData;
/**
	for(int i = 0; i < num_samples; i ++){
		for(int j = 0; j < header.channels; j ++){
			int channel_data_out = vDataOut[j][i];
			for(int k = 0; k < bytes_per_channel; k ++ )
				data_buffer[ j*bytes_per_channel + k ] = (channel_data_out >> (8*k)) & 0x00ff;
			if(postr && fPlot) (*postr) << channel_data_out << " ";
		}
		if(postr && fPlot) (*postr) << "\n";
		if(postr && !fPlot) postr->write( data_buffer, bytes_per_sample );
	}/**/

FINALIZE:
	if( al_play ){
		oal.finalize();
		free(al_buf);
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