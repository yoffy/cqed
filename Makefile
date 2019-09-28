TARGET=cqed
OBJS=main_6_34.o

.PHONY: all
all: $(TARGET)

.PHONY: clean
clean:
	$(RM) $(TARGET) *.o

$(TARGET): $(OBJS)
	$(CXX) -o $@ $^
