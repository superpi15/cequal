/*
 * Cequal: C-based equalizer for educational purpose 
 */

#ifndef CEQUAL_AL_H
#define CEQUAL_AL_H


#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

// https://www.openal.org/documentation/OpenAL_Programmers_Guide.pdf 
#ifdef __APPLE__
#include <OpenAL/al.h>
#include <OpenAL/alc.h>
#elif __linux
#include <AL/al.h>
#include <AL/alc.h>
#include <unistd.h>
#endif

namespace ceq {


class Oal_Man_t {
	int ext_buf_top, sample_per_buf, byte_per_sample;
	int sample_rate;
	int al_push_num;
	ALenum format;
public:
	Oal_Man_t(){
		output_device   = NULL;
		output_context  = NULL;
		al_push_num     = 0;
		sample_rate     = 0;
		init_buf_top    = 0;
		nbuffer = 4;
		sample_per_buf = 22050;
		ext_buf = NULL;
		buffer_list = NULL;
	}
	~Oal_Man_t(){
		if( buffer_list )
			free(buffer_list);
	}
	int nbuffer, init_buf_top;
	int num_al_push() const { return al_push_num; }
	int num_al_buffer() const { return nbuffer; }
	ALCdevice  * output_device;
	ALCcontext * output_context;
	//ALuint internal_buffer;
	ALuint * buffer_list;
	std::string device_name;
	ALuint stream_out;
	bool init();
	void flush();
	void finalize();

	void check_error(const char * str=""){
		ALenum al_error;
		if( AL_NO_ERROR != (al_error = alGetError()) ){
			printf("%s: %s\n", str, alGetString(al_error));
		}
	}

	void play(){
		alSourcePlay(stream_out);
		check_error("play");
	}
	void wait(){
		ALenum playing_state;
		do {
			alGetSourcei(stream_out, AL_SOURCE_STATE, &playing_state);
		} while(AL_PLAYING == playing_state);
	}
	void init_extern_queue( int extern_buf_size, ALenum al_format, ALsizei sample_rate );
	void push_byte( unsigned char sample );
	unsigned char * ext_buf;
private:
	void push_al_buffer();
};

inline void Oal_Man_t::flush(){
	if( 0!=ext_buf_top ){
		memset(ext_buf+ext_buf_top,0, sample_per_buf * byte_per_sample -ext_buf_top);
		push_al_buffer();
	}
}

inline void Oal_Man_t::push_al_buffer(){
	ALuint bid;
	ALenum current_playing_state;
//printf("enter\n");
//		ALenum current_playing_state;
//	alGetSourcei(stream_out, AL_SOURCE_STATE, & current_playing_state);
//	    check_error("alGetSourcei AL_SOURCE_STATE");
//	    if(AL_PLAYING == current_playing_state)
//	    	printf("playing \n");
//	    else {
//	    	printf("not playing \n");
//	    }

    ALenum al_error;
    if( init_buf_top < nbuffer ){
    	bid = buffer_list[init_buf_top++];
    } else 
	    do {
	    	alSourceUnqueueBuffers(stream_out, 1, &bid);
	    } while( AL_NO_ERROR != (al_error = alGetError()) );
	
	check_error("1 unqueue");

	alBufferData( bid, format, ext_buf, sample_per_buf * byte_per_sample, sample_rate);
	//int dummy=0;\
	alBufferData( bid, format, &dummy, 2, this->sample_rate);

	check_error("2 push");
	alSourceQueueBuffers(stream_out, 1, &bid);
	check_error("3 enqueue");
    
    alGetSourcei(stream_out, AL_SOURCE_STATE, & current_playing_state);
    check_error("alGetSourcei AL_SOURCE_STATE");
    if(AL_PLAYING == current_playing_state)
    	;//printf("playing \n");
    else {
    	play();
    }
    //printf("done\n");

	al_push_num++;
	
	check_error();
}

inline void Oal_Man_t::push_byte( unsigned char sample ){
	ext_buf[ext_buf_top++] = sample;

	if( ext_buf_top >= sample_per_buf * byte_per_sample ){
		ALuint bid;
		push_al_buffer();
		ext_buf_top = 0;
	}
}

inline void Oal_Man_t::init_extern_queue( int extern_buf_size, ALenum al_format, int nsample_rate ){
	sample_per_buf = extern_buf_size;
	ext_buf_top = 0;
	sample_rate = nsample_rate;
	format = al_format;
	if( ext_buf )
		free(ext_buf);

	if( AL_FORMAT_MONO8 == format )
		byte_per_sample = 1;
	else
	if( AL_FORMAT_STEREO8 == format || AL_FORMAT_MONO16 == format)
		byte_per_sample = 2;
	else
	if( AL_FORMAT_STEREO16 == format )
		byte_per_sample = 4;
	else {
		assert(false);
	}
	al_push_num = 0;
	ext_buf = (unsigned char*)malloc(sample_per_buf * byte_per_sample);
	memset(ext_buf,0,sizeof(char)*sample_per_buf * byte_per_sample);

//	for(int i = 0; i < nbuffer; i ++)
//		alBufferData( buffer_list[i], al_format, ext_buf, sample_per_buf * byte_per_sample, nsample_rate);
//	alSourceQueueBuffers(stream_out, nbuffer, buffer_list);
//	play();
}

inline bool Oal_Man_t::init(){
	const char * devname = alcGetString(NULL, ALC_DEFAULT_DEVICE_SPECIFIER);
	if( !devname ){
		printf("no default device found\n");
		return false;
	}

	this->device_name = devname;

	output_device  = alcOpenDevice(devname);
	output_context = alcCreateContext(output_device, NULL);
	alcMakeContextCurrent(output_context);

	buffer_list = (ALuint*) malloc(sizeof(ALuint)*nbuffer);
	alGenBuffers(nbuffer, buffer_list);
	check_error("get buffer");
	alGenSources(1, &stream_out);
	check_error("get stream_out");
	ALenum al_error;
	if( AL_NO_ERROR != (al_error = alGetError()) ){
		printf("%s\n", alGetString(al_error));
	}

	return true;
}

inline void Oal_Man_t::finalize(){
	alSourceStopv(1, &stream_out);
	alDeleteSources(1, &stream_out);
	//alDeleteBuffers(nbuffer, buffer_list);

    alcMakeContextCurrent(NULL);
    if(output_context)
    	alcDestroyContext(output_context);
    if(output_device)
    	alcCloseDevice(output_device);
	if( ext_buf )
		free(ext_buf);

}

};

#endif