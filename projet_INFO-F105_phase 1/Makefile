CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -Wpedantic -D_GNU_SOURCE -Werror=all -O2

TARGET = compare_strings

SRCS = compare_strings.cpp utils_fasta.cpp utils_csv.cpp

OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)