#/bin/bash
if [ ! -f build/src/symmetry ]; then
	echo "Run build.sh first";
else
	./build/src/symmetry
fi

