CXXFLAGS=-std=c++11 -g -Wall -Wextra -D NDEBUG -O3 -march=native
LDFLAGS=-g
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
	$(CXX) $(LDFLAGS) -o $@ $^
