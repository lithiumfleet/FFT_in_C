# FFT_in_C

pure c implementation of FFT and a concurrent version.(parallel_dft.c)

### results

~~Not fully testedm but have a glance:~~
~~test N=seq_len=4096*4~~
~~vanilla_dft: 3.734s~~
~~fft_v1: 3.968s~~
~~matlab: 3.018s~~
~~> emmm, sth wrong?~~


### todo

+ [x] fix stack overflow
+ [x] accelerate the programme


### reference

[introduction of fft](https://medium.com/swlh/the-fast-fourier-transform-fft-5e96cf637c38)

[measure time in win-cli](https://stackoverflow.com/questions/673523/how-do-i-measure-execution-time-of-a-command-on-the-windows-command-line)

> not a good chice...

[measure time in matlab](https://ww2.mathworks.cn/help/matlab/matlab_prog/measure-performance-of-your-program.html)
[measure time in c](https://www.geeksforgeeks.org/how-to-measure-time-taken-by-a-program-in-c/)

[ppt of fft](https://gr.xjtu.edu.cn/c/document_library/get_file?folderId=2468402&name=DLFE-138303.pdf)

[concurrent](https://www.analog.com/media/cn/technical-documentation/application-notes/390081591ee_263.pdf)

[mixed-base dft](https://www.bilibili.com/video/BV1Wo4y167rb/)