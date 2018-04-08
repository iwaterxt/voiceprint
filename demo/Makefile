 GXX = gcc

#TOOLCHAIN = /opt/poky/1.5/sysroots/i686-pokysdk-linux/usr/bin/arm-poky-linux-gnueabi

#GXX = $(TOOLCHAIN)/arm-poky-linux-gnueabi-gcc

GXX_FLAGS = -g -msse -msse2 -Wall -fPIC  -Wno-sign-compare -Wno-unused-function  -std=c99
#GXX_FLAGS = -I/opt/poky/1.5/sysroots/cortexa9-vfp-neon-poky-linux-gnueabi/usr/include -Wall -O3 -std=c99 -march=armv7-a -mfpu=neon -ftree-vectorize -mfloat-abi=softfp -ffast-math --sysroot=/opt/poky/1.5/sysroots/cortexa9-vfp-neon-poky-linux-gnueabi
LIBS = -L../asrplus-engine -lvoiceprint -lm -lpthread -lstdc++ -L/opt/OpenBLAS/lib -lopenblas 
#LIBS = -L/opt/poky/1.5/sysroots/cortexa9-vfp-neon-poky-linux-gnueabi/usr/lib -L../asr-engine -lyoungtone -lm -lpthread -lstdc++ -L/opt/poky/1.5/sysroots/cortexa9-vfp-neon-poky-linux-gnueabi -L/opt/poky/1.5/sysroots/cortexa9-vfp-neon-poky-linux-gnueabi/usr/lib/alsa-lib

voiceprint : wave-reader.o voiceprint.o
	$(GXX) $(GXX_FLAGS) wave-reader.o voiceprint.o -o voiceprint $(LIBS)
	rm *.o -rf

voiceprint.o : voiceprint.c
	$(GXX) $(GXX_FLAGS) -o voiceprint.o -c voiceprint.c

wave-reader.o : wave-reader.c
	$(GXX) $(GXX_FLAGS) -o wave-reader.o -c wave-reader.c

clean :
	rm -rf *.o voiceprint
