CXX        = g++
CXXFLAGS   = -std=c++17 -Wall -pthread

SRCDIR     = src
OBJDIR     = obj
TARGET     = main
SRCS       = $(wildcard $(SRCDIR)/*.cpp)
OBJS       = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCS))

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJDIR) $(TARGET)
