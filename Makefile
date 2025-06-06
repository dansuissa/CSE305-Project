CXX        = g++
CXXFLAGS   = -std=c++17 -Wall -pthread

SRCDIR     = src
OBJDIR     = obj
TARGET     = main

# Automatically find all .cpp files in src/
SRCS       = $(wildcard $(SRCDIR)/*.cpp)
# Map each src/foo.cpp to obj/foo.o
OBJS       = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCS))

all: $(TARGET)

# Link all object files into the final executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile rule for each .cpp in src/ to .o in obj/
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJDIR) $(TARGET)
