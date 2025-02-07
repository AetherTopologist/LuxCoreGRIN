# Credits to Blender Cycles

ifeq ($(OS),Windows_NT)
	$(error On Windows, use "cmd //c make.bat" instead of "make")
endif

OS:=$(shell uname -s)

ifndef BUILD_CMAKE_ARGS
	BUILD_CMAKE_ARGS:=
endif

ifndef BUILD_DIR
	BUILD_DIR:=./build
endif

SOURCE_DIR:=$(shell pwd)

ifndef PYTHON
	PYTHON:=python3
endif

all: release

release:
	mkdir -p $(BUILD_DIR)
	cd $(BUILD_DIR) && cmake $(BUILD_CMAKE_ARGS) -DCMAKE_BUILD_TYPE=Release --preset conan-release $(SOURCE_DIR) && cmake --build --preset conan-release

debug:
	mkdir -p $(BUILD_DIR)
	cd $(BUILD_DIR) && cmake $(BUILD_CMAKE_ARGS) -DCMAKE_BUILD_TYPE=Release --preset conan-release $(SOURCE_DIR) && cmake --build --preset conan-debug

clean:
	rm -rf $(BUILD_DIR)

.PHONY: deps
deps:
	$(PYTHON) cmake/make_deps.py

