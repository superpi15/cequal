#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "cequal.h"

#include <fstream>

int main(int argc, char **argv) {
	Ceq_Man_t ceq(argv[1]);
	print_wav_header(ceq.header);

	//std::ofstream ostr("tmp.wav");
	//ceq.writeHeader(ostr);
	return 0;

}
