CXXFLAGS=-std=c++11 -g -Wall -Wextra -O3 -march=native
TARGET=cqed
OBJS=mediancut_util.o \
	 mediancut.o \
	 main.o

.PHONY: all
all: $(TARGET)

.PHONY: clean
clean:
	$(RM) $(TARGET) *.o

$(TARGET): $(OBJS)
	$(CXX) -g -o $@ $^
