vanilla_dft: src/vanilla_dft.c
	gcc -Og src/vanilla_dft.c -o bin/vanilla_dft.exe -Wall

fft_v1: src/fft_v1.c
	gcc -Og src/fft_v1.c -o bin/fft_v1.exe -Wall

test: src/test.c
	gcc -Og src/test.c -o bin/test.exe -Wall
	bin/test.exe