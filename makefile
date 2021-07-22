all: ArchaicSeeker2
	

ArchaicSeeker2: ArchaicSeeker.cpp data.o gzfstream.o matching.o seeker.o
	g++ ArchaicSeeker.cpp data.o gzfstream.o matching.o seeker.o -o ArchaicSeeker2 ${LDFLAG} ${STATIC}

data.o: data.cpp data.hpp
	g++ -c data.cpp ${CFLAG} ${LDFLAG}

gzfstream.o: gzfstream.cpp gzfstream.hpp
	g++ -c gzfstream.cpp ${CFLAG} ${LDFLAG}

matching.o: matching.cpp matching.hpp
	g++ -c matching.cpp ${CFLAG} ${LDFLAG}

seeker.o: seeker.cpp seeker.hpp
	g++ -c seeker.cpp ${CFLAG} ${LDFLAG}

clean:
	rm -rf *.o ArchaicSeeker2

CFLAG=-I/home/unix/kyuan/link/102.InstalledSoftwares/include/ -L/home/unix/kyuan/link/102.InstalledSoftwares/lib/
LDFLAG=-lnlopt -lz
STATIC=-Wl,-Bstatic -lnlopt -lz -Wl,-Bdynamic -lm
