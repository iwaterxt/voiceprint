/*
 * wave-reader.c
 *
 *  Created on: Sep 19, 2017
 *      Author: tao
 */


# include "common.h"

static uint16 read_uint16(FILE* f , int32 swap)
{
	uint16 uinteger16 ;

	if(sizeof(uint16) != fread(&uinteger16 , 1 , sizeof(uint16) , f))  error("error on reading file" , __FILE__ , __LINE__)  ;

	if(swap)  SWAP2(&uinteger16) ;

	return uinteger16 ;
}

static uint32 read_uint32(FILE* f , int32 swap)
{
	uint32 uinteger32 ;

	if(sizeof(uint32) != fread(&uinteger32 , 1 , sizeof(uint32) , f))  error("error on reading file" , __FILE__ , __LINE__)  ;

	if(swap)  SWAP4(&uinteger32) ;

	return uinteger32 ;
}

static int read_quad_bytes(FILE* f , char* quad_bytes)
{
	assert(NULL != quad_bytes) ;

	if(4 != fread(quad_bytes , 1 , 4 , f))  error("error on reading file" , __FILE__ , __LINE__)  ;

	return 1 ;
}

static int expect_quad_bytes(FILE* f , const char* expected_quad_bytes)
{
	char read_bytes[5] ;
	read_bytes[4] = '\0' ;

	read_quad_bytes(f , &read_bytes[0]) ;

	if(0 != strcmp(read_bytes , expected_quad_bytes))
	{
		char err[100] ;
		sprintf(err , "%s got , but %s expected" , read_bytes , expected_quad_bytes) ;
		error(err , __FILE__ , __LINE__)  ;
	}

	return 1 ;
}

/*
 * refer to : https://ccrma.stanford.edu/courses/422/projects/WaveFormat/
 * read wave , and store the the actual sound data
 */
int wav_read(const char* wave , BaseFloat** data , uint32* length)
{
	FILE* f = fopen(wave , "rb") ;
	if(NULL == f)
	{
		char err[100] ;
		sprintf(err , "error : failed to open wave file : %s" , wave) ;
		error(err , __FILE__ , __LINE__) ;
	}

	char read_bytes[5] ;
	read_bytes[4] = '\0' ;

    	/*
	 * the canonical WAVE format starts with the RIFF header :
    	 * 4   ChunkID
     	 * 4   ChunkSize
     	 * 4   Format
        */

	read_quad_bytes(f , &read_bytes[0]) ;

	int32 riff = 0 ;
	if(0 == strcmp(read_bytes , "RIFF"))  riff = 1 ;
	else if(0 == strcmp(read_bytes , "RIFX"))  riff = 0 ;
	else {
		char err[100] ;
		sprintf(err , "chunk ID %s got , but RIFF or RIFX expected" , read_bytes) ;
		error(err , __FILE__ , __LINE__) ;
	}

	/*
	 * the default byte ordering assumed for WAVE data files is little-endian
	 * files written using the big-endian byte ordering scheme have the identifier RIFX instead of RIFF
	 */
	int32 swap = (riff == 1) ? 0 : 1 ;

	uint32 riff_chunk_size = read_uint32(f , swap) ;

	expect_quad_bytes(f , "WAVE") ;

	uint32 riff_chunk_read = 0 ;

	 /* WAVE included in riff_chunk_size */
	riff_chunk_read += 4 ;

	/*
	 * the "fmt " subchunk describes the sound data's format :
     	 * 4   Subchunk1ID
	 * 4   Subchunk1Size
	 * 2   AudioFormat
	 * 2   NumChannels
	 * 4   SampleRate
	 * 4   ByteRate
	 * 2   BlockAlign
	 * 2   BitsPerSample
	 * 2   ExtraParamSize
	 * X   ExtraParams
	 */

	expect_quad_bytes(f , "fmt ") ;

	uint32 subchunk1_size = read_uint32(f , swap) ;

	uint16 audio_format = read_uint16(f , swap) ,
           num_channels = read_uint16(f , swap) ;

	uint32 sample_rate = read_uint32(f , swap) ,
		   byte_rate = read_uint32(f , swap) ;

	uint16 block_align = read_uint16(f , swap) ,
		   bits_per_sample = read_uint16(f , swap) ;

	if(1 != audio_format)
	{
		char err[100] ;
		sprintf(err , "audio format %u got , but only support PCM data with audio format of 1" , audio_format) ;
		error(err , __FILE__ , __LINE__)  ;
	}

	if(subchunk1_size < 16)  error("expect PCM format data to have fmt chunk of size of at least 16" , __FILE__ , __LINE__)  ;
	else  for(uint32 i = 16 ; i < subchunk1_size ; ++i)   fgetc(f) ;

	if(num_channels <= 0)  error("no channels present" , __FILE__ , __LINE__)  ;

	if(bits_per_sample != 8 && bits_per_sample != 16 && bits_per_sample != 32)
	{
		char err[100] ;
		sprintf(err , "bits per sample %u got , but only support 8 / 16 / 32 bits per sample" , bits_per_sample) ;
		error(err , __FILE__ , __LINE__)  ;
	}

	if(byte_rate != sample_rate * num_channels * (bits_per_sample / 8))  error("unexpected byte rate" , __FILE__ , __LINE__)  ;

	if(block_align != num_channels * (bits_per_sample / 8))  error("unexpected block align" , __FILE__ , __LINE__)  ;

	/* size of what we just read , 4 bytes for "fmt " + 4 for subchunk1_size + subchunk1_size itself */
	riff_chunk_read += (8 + subchunk1_size) ;

	/*
	 * just ignore the "fact" chunk , for non-compressed data(which do not support anyway) , it doesn't contain useful information
	 */

	/* the "data" subchunk contains the size of the data and the actual sound :
	 * 4   Subchunk2ID
	 * 4   Subchunk2Size
	 * *   Data
	 */

	expect_quad_bytes(f , "data") ;
	riff_chunk_read += 4 ;

	uint32 data_chunk_size = read_uint32(f , swap) ;
	riff_chunk_read += 4 ;

	char* sound_data = (char*)malloc(data_chunk_size * sizeof(char)) ;
	if(NULL == sound_data)  error("no memory to be allocated" , __FILE__ , __LINE__) ;

	if(data_chunk_size != fread(sound_data , 1 , data_chunk_size , f))  error("error on reading file for actual sound data" , __FILE__ , __LINE__)  ;
	riff_chunk_read += data_chunk_size ;
	char* data_ptr = sound_data ;

	if(riff_chunk_read != riff_chunk_size)
	{
		char err[100] ;
		sprintf(err , "expect %u bytes in RIFF chunk , but get %u bytes(do not support reading multiple data chunks)" , riff_chunk_read , riff_chunk_size) ;
		error(err , __FILE__ , __LINE__) ;
	}

	if(0 != (data_chunk_size % block_align))
	{
		char err[100] ;
		sprintf(err , "data chunk size has unexpected length %u , block align is %u" , data_chunk_size , block_align) ;
		error(err , __FILE__ , __LINE__) ;
	}

	if(0 == data_chunk_size)  error("empty file (no data)" , __FILE__ , __LINE__) ;

	uint32 num_samples = data_chunk_size / block_align ;

	if(NULL == data)  error("empty pointer" , __FILE__ , __LINE__) ;
	if(NULL == *data || (*length) != num_channels * num_samples )
	{
		*length = num_channels * num_samples ;
		*data = (BaseFloat*)malloc((*length) * sizeof(BaseFloat)) ;
		if(NULL == *data)  error("no memory to be allocated" , __FILE__ , __LINE__) ;
	}

	for(uint32 i = 0 ; i < num_samples ; i++)
	{
		for(uint16 j = 0 ; j < num_channels ; j++)
		{
			switch(bits_per_sample)
			{
				case 8 :
					(*data)[j * num_samples + i] = (BaseFloat)(*data_ptr) ;
					data_ptr++ ;
					break ;
				case 16 :
				{
					int16 integer16 = *((int16*)(data_ptr)) ;
					if(swap)  SWAP2(&integer16) ;
					(*data)[j * num_samples + i]  = (BaseFloat)(integer16) ;
					data_ptr += 2 ;
					break ;
				}
				case 32 :
				{
					int32 integer32 = *((int32*)(data_ptr)) ;
					if(swap)  SWAP4(&integer32) ;
					(*data)[j * num_samples + i]  = (BaseFloat)(integer32) ;
					data_ptr += 4 ;
					break ;
				}
				default :
					error("not supported bits per sample" , __FILE__ , __LINE__) ;    // already checked
			}
		}
	}

	free(sound_data) ;

	fclose(f) ;

	return 1 ;
}



