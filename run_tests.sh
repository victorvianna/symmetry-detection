#/bin/bash
if [ ! -f build/test/symmetry_tests ]; then
	echo "Run build.sh first";
else
	./build/test/symmetry_tests
fi

