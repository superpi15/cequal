/*
 * Cequal: C-based equalizer for educational purpose 
 */

#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include "cequal.h"

int main(int argc, char **argv) {
	if( argc < 2 ){
		printf("Usage: ./program <*.WAV>\n");
		return 0;
	}
//	ceq::Oal_Man_t oal;
//	oal.init();
//	oal.finalize();
//	return 0;
	ceq::Ceq_Man_t ceq(argv[1]);
	//ceq::print_wav_header(ceq.header);
	//ceq.scanSample();
	//return 0;
	//ceq.config_al_play(0);
	//ceq.config_nl_play(0);
	std::ofstream ostr1("tmp.wav");\
	ceq.run_prefetch_filter(&ostr1);\
	ostr1.close();
	//std::ofstream ostr2("tmp.plot");\
	ceq.run_prefetch_filter(&ostr2,1);\
	ostr2.close();
	
	return 0;
}
