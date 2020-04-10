# ./src: source code path
# ./obj: object file path
# ./bin: binary file path

NAME = GIFT
GCC = g++
CFLAGS = -lboost_thread -lgsl -lgslcblas -lstdc++ -Wall -O3 
OBJDIR = obj
SRCDIR = src
BINDIR = bin
SRCS = $(notdir $(wildcard $(SRCDIR)/*.cpp))
OBJS = $(patsubst %.cpp, $(OBJDIR)/%.o, $(SRCS))
TARGET = $(BINDIR)/$(NAME)

$(TARGET):$(OBJS)
	-mkdir -p ${BINDIR}
	$(GCC) -g $(CFLAGS) -o $@ $^
$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	-mkdir -p ${OBJDIR}
	$(GCC) -c -o $@ $<

.PHONY: install
install: ${TARGET}

.PHONY: clean
clean:
	-rm $(TARGET) $(OBJDIR)/*.o
