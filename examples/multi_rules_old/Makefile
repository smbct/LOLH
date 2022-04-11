EXEC=prog

SRC_FILES_LOLH := $(wildcard ../../c++/src/*.cpp)
SRC_FILES_LOLH := $(patsubst ../../c++/src/main.cpp,,$(SRC_FILES_LOLH))
OBJ_FILES_LOLH := $(patsubst ../../c++/src/%.cpp, ../../c++/obj/%.o,$(SRC_FILES_LOLH))

SRC_FILES := $(wildcard *.cpp)
OBJ_FILES := $(patsubst %.cpp,%.o,$(SRC_FILES))

.PHONY = clean


all: main

run: main
	./main

../../c++/obj/%.o: ../../c++/src/%.cpp
	g++ -std=c++17 -o3 $< -c -o $@

./%.o: ./%.cpp
	g++ -std=c++17 -o3 -DGHOST_TRACE $< -c -o $@

main: $(OBJ_FILES) $(OBJ_FILES_LOLH)
	g++ -std=c++17 -o3 $(OBJ_FILES) $(OBJ_FILES_LOLH) -lghost -DGHOST_TRACE -pthread -o main

clean:
	rm *.o ./main
