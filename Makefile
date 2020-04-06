DEBUG = 0
DUMPDATA = 1

CC = gcc
EXEC = ms2fft0
COMMON = -I./libmseed/ -I./fftw_module/ -I.
CFLAGS =  -Wall
LDFLAGS = -L./libmseed -Wl,-rpath,./libmseed
LDLIBS = -Wl,-Bstatic -lmseed -Wl,-Bdynamic -lm -lfftw3

OBJS = main.o standard_deviation.o fftw_module/fft.o

ifeq ($(DEBUG), 1)
CFLAGS += -O0 -g -DDEBUG=1
endif

ifeq ($(DUMPDATA), 1)
CFLAGS += -DDUMPDATA=1
endif

.PHONY: all clean

all: $(EXEC)

$(EXEC): $(OBJS)
	$(MAKE) -C libmseed/ static
	$(CC) $(COMMON) $(CFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

%.o: %.c
	$(CC) $(COMMON) $(CFLAGS) -c $< -o $@

clean:
	$(MAKE) -C libmseed/ clean
	rm -rf $(OBJS) $(EXEC)
