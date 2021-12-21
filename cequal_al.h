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
	int ext_buf_top, ext_buf_size, ext_buf_unit_size;
	int sample_rate;
	int al_push_num;
	ALenum format;
public:
	Oal_Man_t(){
		output_device   = NULL;
		output_context  = NULL;
		al_push_num     = 0;
		sample_rate     = 0;
		nbuffer = 4;
		ext_buf_size = 22050;
		ext_buf = NULL;
		buffer_list = (ALuint*) malloc(sizeof(ALuint)&nbuffer);
	}
	~Oal_Man_t(){
		free(buffer_list);
	}
	int nbuffer;
	ALCdevice  * output_device;
	ALCcontext * output_context;
	//ALuint internal_buffer;
	ALuint * buffer_list;
	std::string device_name;
	ALuint stream_out;
	bool init();
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
};

inline void Oal_Man_t::push_byte( unsigned char sample ){
	ext_buf[ext_buf_top++] = sample;

	if( ext_buf_top >= ext_buf_size * ext_buf_unit_size ){
		ALuint bid;
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
	    do {
	    	alSourceUnqueueBuffers(stream_out, 1, &bid);
	    } while( AL_NO_ERROR != (al_error = alGetError()) );
		
		check_error("1 unqueue");

		alBufferData( bid, format, ext_buf, ext_buf_size * ext_buf_unit_size, sample_rate);
		//int dummy=0;\
		alBufferData( bid, format, &dummy, 2, this->sample_rate);

		check_error("2 push");
		alSourceQueueBuffers(stream_out, 1, &bid);
		check_error("3 enqueue");
	    
//	    alGetSourcei(stream_out, AL_SOURCE_STATE, & current_playing_state);
//	    check_error("alGetSourcei AL_SOURCE_STATE");
//	    if(AL_PLAYING == current_playing_state)
//	    	printf("playing \n");
//	    else {
//	    	play();
//	    }
//	    printf("done\n");
		ext_buf_top = 0;
		al_push_num++;
	}
	check_error();
}

inline void Oal_Man_t::init_extern_queue( int extern_buf_size, ALenum al_format, int nsample_rate ){
	ext_buf_size = extern_buf_size;
	ext_buf_top = 0;
	sample_rate = nsample_rate;
	format = al_format;
	if( ext_buf )
		free(ext_buf);

	if( AL_FORMAT_MONO8 == format )
		ext_buf_unit_size = 1;
	else
	if( AL_FORMAT_STEREO8 == format || AL_FORMAT_MONO16 == format)
		ext_buf_unit_size = 2;
	else
	if( AL_FORMAT_STEREO16 == format )
		ext_buf_unit_size = 4;
	else {
		assert(false);
	}
	al_push_num = 0;
	ext_buf = (unsigned char*)malloc(ext_buf_size * ext_buf_unit_size);
	memset(ext_buf,0,sizeof(char)*ext_buf_size * ext_buf_unit_size);


	for(int i = 0; i < nbuffer; i ++)
		alBufferData( buffer_list[i], al_format, ext_buf, ext_buf_size * ext_buf_unit_size, nsample_rate);
	alSourceQueueBuffers(stream_out, nbuffer, buffer_list);
	play();
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