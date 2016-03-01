# ./src: source code path
# ./obj: object file path
# ./bin: binary file path

NAME = GIFT
GCC = g++
CFLAGS = -lboost_thread -lgsl -lgslcblas -lstdc++ -Wall -O3 
OBJDIR = obj
SRCDIR = src
BIN = bin
SRCS = $(notdir $(wildcard $(SRCDIR)/*.cpp))
OBJS = $(patsubst %.cpp, $(OBJDIR)/%.o, $(SRCS))
TARGET = $(BIN)/$(NAME)

$(TARGET):$(OBJS)
	$(GCC) -g $(CFLAGS) -o $@ $^
$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(GCC) -c -o $@ $<

.PHONY: clean
clean:
	-rm $(TARGET) $(OBJDIR)/*.o
									  
