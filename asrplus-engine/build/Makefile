CC = gcc


CC_OPTIONS = -I../ -O3 -g3 -Wall -c -std=c++11 -fmessage-length=0 

INCLUDE = -I/opt/OpenBLAS/include

LIBS = -lstdc++ -lpthread -lm -L/opt/OpenBLAS/lib -lopenblas

OBJ = common.o io.o posterior.o mult-thread.o data-interface.o math.o queue.o cmvn.o feature-common.o feature-function.o mfcc.o \
diag-gmm.o ivector-extractor.o Matrix.o SpMatrix.o TpMatrix.o Vector.o srfft.o asr-plusplus.o engine.o vpcontainer.o plda.o vad.o


voiceprint.a :  $(OBJ)
	ar cru libvoiceprint.a $(OBJ)
	#$(RANLIB) voiceprint.a
	chmod 775 libvoiceprint.a
	mv libvoiceprint.a ..
	rm *.o

# base
common.o : ../base/common.cc
	$(CC) $(CC_OPTIONS) ../base/common.cc
io.o : ../base/io.cc
	$(CC) $(CC_OPTIONS) ../base/io.cc
posterior.o : ../base/posterior.cc
	$(CC) $(CC_OPTIONS) $(INCLUDE) ../base/posterior.cc

# util
mult-thread.o : ../util/mult-thread.cc
	$(CC) $(CC_OPTIONS) ../util/mult-thread.cc

data-interface.o : ../util/data-interface.cc
	$(CC) $(CC_OPTIONS) ../util/data-interface.cc

math.o : ../util/math.cc
	$(CC) $(CC_OPTIONS) $(INCLUDE) ../util/math.cc

queue.o : ../util/queue.cc
	$(CC) $(CC_OPTIONS) ../util/queue.cc


# feature
cmvn.o : ../feature/cmvn.cc
	$(CC) $(CC_OPTIONS) ../feature/cmvn.cc

feature-common.o : ../feature/feature-common.cc
	$(CC) $(CC_OPTIONS) $(INCLUDE) ../feature/feature-common.cc

feature-function.o : ../feature/feature-function.cc
	$(CC) $(CC_OPTIONS) $(INCLUDE) ../feature/feature-function.cc

mfcc.o : ../feature/mfcc.cc
	$(CC) $(CC_OPTIONS) $(INCLUDE) ../feature/mfcc.cc

# gmm
diag-gmm.o : ../gmm/diag-gmm.cc
	$(CC) $(CC_OPTIONS) $(INCLUDE) ../gmm/diag-gmm.cc

# ivector
ivector-extractor.o : ../ivector/ivector-extractor.cc
	$(CC) $(CC_OPTIONS) $(INCLUDE) ../ivector/ivector-extractor.cc
plda.o : ../ivector/plda.cc
	$(CC) $(CC_OPTIONS) $(INCLUDE) ../ivector/plda.cc
vad.o : ../ivector/vad.cc
	$(CC) $(CC_OPTIONS) $(INCLUDE) ../ivector/vad.cc

# matrix
Matrix.o : ../matrix/Matrix.cc
	$(CC) $(CC_OPTIONS) $(INCLUDE) ../matrix/Matrix.cc

SpMatrix.o : ../matrix/SpMatrix.cc
	$(CC) $(CC_OPTIONS) $(INCLUDE) ../matrix/SpMatrix.cc

TpMatrix.o : ../matrix/TpMatrix.cc
	$(CC) $(CC_OPTIONS) $(INCLUDE) ../matrix/TpMatrix.cc

Vector.o : ../matrix/Vector.cc
	$(CC) $(CC_OPTIONS) $(INCLUDE) ../matrix/Vector.cc

srfft.o : ../matrix/srfft.cc
	$(CC) $(CC_OPTIONS) $(INCLUDE) ../matrix/srfft.cc
#engine

asr-plusplus.o : ../engine/asr-plusplus.cc
	$(CC) $(CC_OPTIONS) $(INCLUDE) ../engine/asr-plusplus.cc

engine.o : ../engine/engine.cc
	$(CC) $(CC_OPTIONS) $(INCLUDE) ../engine/engine.cc

vpcontainer.o : ../engine/vpcontainer.cc
	$(CC) $(CC_OPTIONS) $(INCLUDE) ../engine/vpcontainer.cc


clean :
	rm -rf *.o ../*.a
