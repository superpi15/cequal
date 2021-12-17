#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "cequal.h"

#include <fstream>

int main(int argc, char **argv) {
	Ceq_Man_t ceq(argv[1]);
	print_wav_header(ceq.header);

	std::ofstream ostr1("tmp.wav");
	std::ofstream ostr2("tmp.plot");
	ceq.runFilter(&ostr1);
	ceq.runFilter(&ostr2,1);
	ostr1.close();
	ostr2.close();
	//std::ofstream ostr("tmp.plot");\
	ceq.runFilter(&ostr,1);
	//ceq.writeHeader(ostr);
	
	return 0;
}
