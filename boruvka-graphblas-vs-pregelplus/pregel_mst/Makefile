CCOMPILE = mpic++
CPPFLAGS = -I$(HADOOP_HOME)/include -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux -I. -Wno-deprecated -O2
LIB = -L$(HADOOP_HOME)/lib/native
LDFLAGS = -lhdfs

BUILD_TYPE ?= regular

ifeq ($(BUILD_TYPE),mpip)
    LDFLAGS += -lmpiP -lbfd
endif

.PHONY: all clean

all: run

run: run.cpp
	$(CCOMPILE) run.cpp $(CPPFLAGS) $(LIB) $(LDFLAGS) -o run

clean:
	-rm -f run