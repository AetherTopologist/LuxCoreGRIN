# Credits to Blender Cycles

# Environment variables that control the build:
#   LUX_SOURCE_DIR - The path to the top level of the source tree (default: .)
#   LUX_BINARY_DIR - The path to the top level of the build tree (default: ./out)
#   LUX_BUILD_TYPE - Specifies the build type (default: Release)


ifndef PYTHON
	PYTHON:=python3
endif

LUX-CMAKE = $(PYTHON) build-helpers/make/cmake.py

build-targets = pyluxcore luxcoreui luxcoreconsole luxcore doc

.PHONY: deps list-presets config luxcore pyluxcore luxcoreui luxcoreconsole install clean clear doc

all: luxcore pyluxcore luxcoreui luxcoreconsole

deps:
	$(LUX-CMAKE) deps

list-presets:
	$(LUX-CMAKE) list-presets

config:
	$(LUX-CMAKE) config

$(build-targets): %: config
	$(LUX-CMAKE) build-and-install $*

install:
	$(LUX-CMAKE) install all

package:
	$(LUX-CMAKE) build-and-install package

clean:
	$(LUX-CMAKE) clean

clear:
	$(LUX-CMAKE) clear
