/*
 * Cequal: C-based equalizer for educational purpose 
 */

#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "cequal.h"

#include <fstream>

int main(int argc, char **argv) {
	if( argc < 2 ){
		printf("Usage: ./program <*.WAV>\n");
		return 0;
	}
	Ceq_Man_t ceq(argv[1]);
	print_wav_header(ceq.header);

	std::ofstream ostr1("tmp.wav");\
	ceq.run_prefetch_filter(&ostr1);\
	ostr1.close();
	//std::ofstream ostr2("tmp.plot");\
	ceq.run_prefetch_filter(&ostr2,1);\
	ostr2.close();
	//std::ofstream ostr("tmp.plot");\
	ceq.runFilter(&ostr,1);
	//ceq.writeHeader(ostr);
	
	return 0;
}
