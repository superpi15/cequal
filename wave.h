/*
 * Cequal: C-based equalizer for educational purpose 
 * The header definition is modified from http://truelogic.org/wordpress/2015/09/04/parsing-a-wav-file-in-c/ 
 */
#ifndef WAVE_H
#define WAVE_H

// WAVE file header format
struct WavHeader_t {
	unsigned char riff[4];						// RIFF string
	unsigned int overall_size	;				// overall size of file in bytes
	unsigned char wave[4];						// WAVE string
	unsigned char fmt_chunk_marker[4];			// fmt string with trailing null char
	unsigned int length_of_fmt;					// length of the format data
	unsigned int format_type;					// format type. 1-PCM, 3- IEEE float, 6 - 8bit A law, 7 - 8bit mu law
	unsigned int channels;						// no.of channels
	unsigned int sample_rate;					// sampling rate (blocks per second)
	unsigned int byterate;						// SampleRate * NumChannels * BitsPerSample/8
	unsigned int block_align;					// NumChannels * BitsPerSample/8
	unsigned int bits_per_channel;				// bits per sample, 8- 8bits, 16- 16 bits etc
	unsigned char data_chunk_header [4];		// DATA string or FLLR string
	unsigned int data_size;						// NumSamples * NumChannels * BitsPerSample/8 - size of the next chunk that will be read
};

inline void print_wav_header( WavHeader_t& header ){
	const char * format_name;
	switch(header.format_type){
		case 1: format_name = "PCM"; break;
		case 6: format_name = "A-law"; break;
		case 7: format_name = "Mu-law"; break;
		default: format_name = "unknown";
	}
	#define PRINT_CHAR_DATA(B1,B2,NAME,CHAR0) printf("(%d-%d) %s: ",B1,B2,NAME);for(int i = 0; i < B2-B1+1; i ++) putchar(CHAR0[i]);printf("\n");
	PRINT_CHAR_DATA(1,4,"",header.riff);
	printf("(5-8) Overall size: bytes:%u, Kb:%u \n", header.overall_size, header.overall_size/1024);
	PRINT_CHAR_DATA( 9,12,"Wave marker",header.wave);
	PRINT_CHAR_DATA(13,16,"Fmt marker",header.fmt_chunk_marker);
	printf("(17-20) Length of Fmt header: %u \n", header.length_of_fmt);
	PRINT_CHAR_DATA(21,22,"Format type", format_name);
	printf("(23-24) Channels: %u \n", header.channels);
	printf("(25-28) Sample rate: %u\n", header.sample_rate);
	printf("(29-32) Byte Rate: %u , Bit Rate:%u\n", header.byterate, header.byterate*8);
	printf("(33-34) Block Alignment: %u \n", header.block_align);
	printf("(35-36) Bits per channel: %u \n", header.bits_per_channel);
	PRINT_CHAR_DATA(37,40, "Data Marker", header.data_chunk_header);
	printf("(41-44) Size of data chunk: %u \n", header.data_size);
	long num_samples = (8 * header.data_size) / (header.channels * header.bits_per_channel);
	printf("Number of samples:%lu \n", num_samples);

	long size_of_each_sample = (header.channels * header.bits_per_channel) / 8;
	printf("Size of each sample:%ld bytes\n", size_of_each_sample);

	// calculate duration of file
	float duration_in_seconds = (float) header.overall_size / header.byterate;
	printf("Approx.Duration in seconds=%f\n", duration_in_seconds);

}




inline int read_wav_header(char * filename, WavHeader_t& header, int fVerbose ) {

	unsigned char buffer4[4];
	unsigned char buffer2[2];
//	// get file path
//	char cwd[1024];
//	if (getcwd(cwd, sizeof(cwd)) != NULL) {
//
//		strcpy(filename, cwd);
//
//		// get filename from command line
//		if (argc < 2) {
//		  printf("No wave file specified\n");
//		  return 0;
//		}
//
//		strcat(filename, "/");
//		strcat(filename, argv[1]);
//		printf("%s\n", filename);
//	}
	// open file
	if( fVerbose )
		printf("Opening  file..\n");
	FILE * ptr = fopen(filename, "rb");
	if (ptr == NULL) {
		printf("Error opening file\n");
		exit(1);
	}
 
	int read = 0;

	// read header parts

	read = fread(header.riff, sizeof(header.riff), 1, ptr);
	if( fVerbose )
		printf("(1-4): %s \n", header.riff); 

	read = fread(buffer4, sizeof(buffer4), 1, ptr);
	if( fVerbose )
		printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);
 
	// convert little endian to big endian 4 byte int
	header.overall_size  = buffer4[0] | 
						(buffer4[1]<<8) | 
						(buffer4[2]<<16) | 
						(buffer4[3]<<24);

	if( fVerbose )
		printf("(5-8) Overall size: bytes:%u, Kb:%u \n", header.overall_size, header.overall_size/1024);

	read = fread(header.wave, sizeof(header.wave), 1, ptr);
	if( fVerbose )
		printf("(9-12) Wave marker: %s\n", header.wave);

	read = fread(header.fmt_chunk_marker, sizeof(header.fmt_chunk_marker), 1, ptr);
	if( fVerbose )
		printf("(13-16) Fmt marker: %s\n", header.fmt_chunk_marker);

	read = fread(buffer4, sizeof(buffer4), 1, ptr);
	if( fVerbose )
		printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);

	// convert little endian to big endian 4 byte integer
	header.length_of_fmt = buffer4[0] |
							(buffer4[1] << 8) |
							(buffer4[2] << 16) |
							(buffer4[3] << 24);
	if( fVerbose )
		printf("(17-20) Length of Fmt header: %u \n", header.length_of_fmt);

	read = fread(buffer2, sizeof(buffer2), 1, ptr); 
	if( fVerbose )
		printf("%u %u \n", buffer2[0], buffer2[1]);

	header.format_type = buffer2[0] | (buffer2[1] << 8);
	char format_name[10] = "";
	if (header.format_type == 1)
		strcpy(format_name,"PCM"); 
	else if (header.format_type == 6)
		strcpy(format_name, "A-law");
	else if (header.format_type == 7)
		strcpy(format_name, "Mu-law");

	if( fVerbose )
		printf("(21-22) Format type: %u %s \n", header.format_type, format_name);

	read = fread(buffer2, sizeof(buffer2), 1, ptr);
	if( fVerbose )
		printf("%u %u \n", buffer2[0], buffer2[1]);

	header.channels = buffer2[0] | (buffer2[1] << 8);
	if( fVerbose )
		printf("(23-24) Channels: %u \n", header.channels);

	read = fread(buffer4, sizeof(buffer4), 1, ptr);
	if( fVerbose )
		printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);

	header.sample_rate = buffer4[0] |
						(buffer4[1] << 8) |
						(buffer4[2] << 16) |
						(buffer4[3] << 24);

	if( fVerbose )
		printf("(25-28) Sample rate: %u\n", header.sample_rate);

	read = fread(buffer4, sizeof(buffer4), 1, ptr);
	if( fVerbose )
		printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);

	header.byterate  = buffer4[0] |
						(buffer4[1] << 8) |
						(buffer4[2] << 16) |
						(buffer4[3] << 24);
	if( fVerbose )
		printf("(29-32) Byte Rate: %u , Bit Rate:%u\n", header.byterate, header.byterate*8);

	read = fread(buffer2, sizeof(buffer2), 1, ptr);
	if( fVerbose )
		printf("%u %u \n", buffer2[0], buffer2[1]);

	header.block_align = buffer2[0] | (buffer2[1] << 8);
	if( fVerbose )
		printf("(33-34) Block Alignment: %u \n", header.block_align);

	read = fread(buffer2, sizeof(buffer2), 1, ptr);
	if( fVerbose )
		printf("%u %u \n", buffer2[0], buffer2[1]);

	header.bits_per_channel = buffer2[0] |
					(buffer2[1] << 8);
	if( fVerbose )
		printf("(35-36) Bits per channel: %u \n", header.bits_per_channel);

	read = fread(header.data_chunk_header, sizeof(header.data_chunk_header), 1, ptr);
	if( fVerbose )
		printf("(37-40) Data Marker: %s \n", header.data_chunk_header);

	read = fread(buffer4, sizeof(buffer4), 1, ptr);
	if( fVerbose )
		printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);

	header.data_size = buffer4[0] |
				(buffer4[1] << 8) |
				(buffer4[2] << 16) | 
				(buffer4[3] << 24 );
	if( fVerbose )
		printf("(41-44) Size of data chunk: %u \n", header.data_size);

	// calculate no.of samples
	long num_samples = (8 * header.data_size) / (header.channels * header.bits_per_channel);
	if( fVerbose )
		printf("Number of samples:%lu \n", num_samples);

	long size_of_each_sample = (header.channels * header.bits_per_channel) / 8;
	if( fVerbose )
		printf("Size of each sample:%ld bytes\n", size_of_each_sample);

	// calculate duration of file
	float duration_in_seconds = (float) header.overall_size / header.byterate;
	if( fVerbose ){
		printf("Approx.Duration in seconds=%f\n", duration_in_seconds);
	}
	fclose(ptr);
	return 0;

}



#endif