#ifndef CEQUAL_H
#define CEQUAL_H

#include <string>
#include <iostream>
#include <vector>
#include <assert.h>
#include "wave.h"
#include "math.h"

#define CEQUAL_PI 3.1415926 

class Ceq_Man_t {
public:
	std::string filename;
	WavHeader_t header;
	int nBand;
	Ceq_Man_t(char * nfilename){
		memset(this,0,sizeof(Ceq_Man_t));
		setDefault();
		filename = nfilename;
		read_wav_header(nfilename,header,0);
	}
	void setDefault();
	void dumpSample();
	void writeHeader( std::ostream& ostr );

	double time_stamp;                // keep progress for dynamic adjustment 
	std::vector<double> vBandGain;    // gain of each band 
	std::vector<double> vBandOmega;   // angular frequency of a band 

	void runFilter(std::ostream * postr=NULL, int fPlot=0);
	int window_size;
};

inline void Ceq_Man_t::setDefault(){
	nBand = 10;
	window_size = 10;
	int max_freq = 44100;
	vBandGain .resize(nBand,(double)1.0f/nBand);
	//vBandGain .resize(nBand,(double)1.0f/log2(nBand));
	//vBandGain .resize(nBand,1.0f);
	vBandOmega.resize(nBand);
	for(int i = 0; i < nBand; i ++)
		vBandOmega[i] = (32<<i) * CEQUAL_PI;
}

inline void Ceq_Man_t::runFilter(std::ostream * postr, int fPlot){

	printf("Band gain: ");
	for(int i = 0; i < nBand; i ++)
		printf("%9.3f ", vBandGain[i]);
	printf("\n");
	printf("Band freq: ");
	for(int i = 0; i < nBand; i ++)
		printf("%9.3f ", vBandOmega[i]/CEQUAL_PI);
	printf("\n");
	FILE * fptr = fopen(filename.c_str(),"r");
	if( !fptr ){
		printf("Cannot open \'%s\'\n", filename.c_str());
		return;
	}
	if( fseek(fptr, 44, SEEK_SET) ){
		printf("Invalid header of \'%s\'\n", filename.c_str());
		return;
	}
	long bytes_per_sample = (header.channels * header.bits_per_channel) / 8;
	long bytes_per_channel = (bytes_per_sample / header.channels);
	long num_samples = (8 * header.data_size) / (header.channels * header.bits_per_channel);
	char * data_buffer = (char*) malloc(bytes_per_sample);

	int ** data_ring = (int**) malloc(header.channels*sizeof(void*));
	for(int i = 0; i < header.channels; i++ )
		data_ring[i] = (int*) malloc(this->window_size*sizeof(int));

	if(postr && !fPlot) writeHeader(*postr);

	std::vector<double> vFreqAmpR(this->nBand,0.0f);
	std::vector<double> vFreqAmpI(this->nBand,0.0f);
	std::vector<double> vTimeAmp(this->window_size,0.0f);

	int i;
	//double eplasedTime = 0;
	for(i = 0; i < num_samples; i ++){
		//double eplasedTime = (double) ((int)i%header.sample_rate)/header.sample_rate;
		//eplasedTime = (double) i/header.sample_rate;
		//printf("t %8.6f: ", eplasedTime);

		int data_top, window_begin;
		//double time_begin;
		data_top = i % this->window_size;
		window_begin = (i+this->window_size+1) % this->window_size;
		//time_begin = (double) (i-this->window_size+1)/ header.sample_rate;
		if(!fread(data_buffer,bytes_per_sample,1,fptr)){
			printf("Data loading mismatch\n");
			goto FINALIZE;
		}
		for(int j = 0; j < header.channels; j ++){
			int channel_data_in = 0, channel_data_out = 0;
			for(int k = 0; k < bytes_per_channel; k ++ )
				channel_data_in |= (data_buffer[ j*bytes_per_channel + k ]&(k==bytes_per_channel-1? ~0: 0x00ff)) << (8*k);
			if( 1 == bytes_per_channel )
				channel_data_in -= 128;
			data_ring[j][data_top] = channel_data_in;


			for(int k = 0; k < nBand; k ++){
				vFreqAmpI[k] = 0;
				vFreqAmpR[k] = 0;
				for(int l = 0; l < this->window_size; l ++ ){
					int idx = (window_begin+l) % this->window_size;
					int cur_time = (double)(i-this->window_size+1+l) / header.sample_rate;
					vFreqAmpI[k] += (double) data_ring[j][idx] * sin(cur_time * vBandOmega[k]);// * vBandGain[k];
					vFreqAmpR[k] += (double) data_ring[j][idx] * cos(cur_time * vBandOmega[k]);// * vBandGain[k];
				}
			}
			
			for(int k = 0; k < this->window_size; k ++){
				//int idx = (window_begin+k) % this->window_size;
				int cur_time = (double)(i-this->window_size+1+k) / header.sample_rate;
				//vTimeAmp[k] = 0;
				double Real = 0, Img = 0;
				for(int l = 0; l < nBand; l ++){
					Real += (double) 
						vFreqAmpR[l] * cos( - vBandOmega[l] * cur_time ) +
						vFreqAmpI[l] * sin( - vBandOmega[l] * cur_time ) ;
//					Img += (double) 
//						vFreqAmpI[l] * cos( - vBandOmega[l] * cur_time ) +
//						vFreqAmpR[l] * sin( - vBandOmega[l] * cur_time ) ;
				}
				//vTimeAmp[k] = (double) (0<Real?1:-1) * sqrt(Real*Real + Img*Img) / nBand;
				vTimeAmp[k] = (double)Real/nBand;
				if(postr && fPlot) (*postr)<< vTimeAmp[k] << " ";
			}
			if(postr && fPlot) (*postr) << "\n";
			channel_data_out = vTimeAmp[this->window_size-1];
			//if(postr && fPlot) (*postr)<< channel_data_out << "\n";
			for(int k = 0; k < bytes_per_channel; k ++ )
				data_buffer[ j*bytes_per_channel + k ] = (channel_data_out >> (8*k)) & 0x00ff;

		}
		if(postr && !fPlot) postr->write( data_buffer, bytes_per_sample );
/**
		for(int j = 0; j < nBand; j ++){
			//vBandAmpCache[j] = (sin( eplasedTime * vBandOmega[j] )+1.0f)/2.0f;
			vBandAmpCacheF[j] = sin( eplasedTime * vBandOmega[j] );
			vBandAmpCacheB[j] = sin( -eplasedTime * vBandOmega[j] );
		}

		for(int j = 0; j < header.channels; j ++){
			int channel_data_in = 0, channel_data_out = 0;
			for(int k = 0; k < bytes_per_channel; k ++ ){
//				printf(" %4d ", data_buffer[ j*bytes_per_channel + k ] );
				channel_data_in |= (data_buffer[ j*bytes_per_channel + k ]&(k==bytes_per_channel-1? ~0: 0x00ff)) << (8*k);
			}
			if( 1 == bytes_per_channel )
				channel_data_in -= 128;
			
//			printf("(");
			double tmp_data_out = 0, tmp_data_out2 = 0;;
			for(int k = 0; k < nBand; k ++){
//				printf("  %5.2f  ", vBandAmpCache[k] * vBandGain[k] );
				vBandChannelCache[k] = (double)channel_data_in * vBandAmpCacheF[k];// * vBandGain[k];
			}
			for(int k = 0; k < nBand; k ++){
//				printf("  %5.2f  ", vBandAmpCache[k] * vBandGain[k] );
				tmp_data_out2 += (double)tmp_data_out * vBandAmpCacheB[k];// * vBandGain[k];
			}
			tmp_data_out = (double) tmp_data_out2/this->nBand;
			//printf(" %7d %10.6f \n",channel_data_in, tmp_data_out);
			channel_data_out = tmp_data_out;
			if(postr && fPlot) (*postr)<< channel_data_out << "\n";
//			(*postr) << channel_data_out << "\n";

			//printf(")");
//			printf(": in %6d out %6d :",channel_data_in, channel_data_out);

			for(int k = 0; k < bytes_per_channel; k ++ ){
				data_buffer[ j*bytes_per_channel + k ] = (channel_data_out >> (8*k)) & 0x00ff;
//				printf(" %4d ", data_buffer[ j*bytes_per_channel + k ] );
			}
//			printf("\n");
		}
		if(postr && !fPlot) postr->write( data_buffer, bytes_per_sample );
		/**/
	}

FINALIZE:
	printf("Scanned samples: %d\n", i);
	for(int i = 0; i < header.channels; i++ )
		free(data_ring[i]);
	free(data_ring);
	free(data_buffer);
	fclose(fptr);
}


inline void Ceq_Man_t::dumpSample(){
	FILE * fptr = fopen(filename.c_str(),"r");
	if( !fptr ){
		printf("Cannot open \'%s\'\n", filename.c_str());
		return;
	}
	if( fseek(fptr, 44, SEEK_SET) ){
		printf("Invalid header of \'%s\'\n", filename.c_str());
		return;
	}
	long bytes_per_sample = (header.channels * header.bits_per_channel) / 8;
	long bytes_in_each_channel = (bytes_per_sample / header.channels);
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

#endif