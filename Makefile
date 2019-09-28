CXXFLAGS=-g -Wall -Wextra
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
