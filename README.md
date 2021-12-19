# cequal: C-based equalizer for educational purpose 

This equalizer demonstrates a typical conversion between time and frequency -domain, often applied to audio processing. 

* The current implementation uses smooth sine wave for sampling, in order to show numerical concept of DFT. 
* The square wave and thus FFT may be supported later. 

This document is exported from: [https://github.com/superpi15/cequal/](https://github.com/superpi15/cequal/)

## system prerequisite 

* Open MP
* g++ 

## compile
```shell
sh compile.sh
```

## run 

```shell
./cequal <YOUR_WAV_FILE>
```
Output file is named *tmp.wav*.

--- --- ---

## implementation notes 

![](https://latex.codecogs.com/svg.latex?{e^{i\phi}=\cos{\phi}+i\sin{\phi}})

![](https://latex.codecogs.com/svg.latex?{2\pi}) radians 
![](https://latex.codecogs.com/svg.latex?{=360\degree})

1Hz ![](https://latex.codecogs.com/svg.latex?{\approx}) 6.28 rad/sec

In C math library, sine function takes radians as its input. 


A sample is a data of a time step. 
Given sample rate ![](https://latex.codecogs.com/svg.latex?{\alpha}) 
and sample index ![](https://latex.codecogs.com/svg.latex?{k}), 
the elapsed time is <br/>
<p align="center">
<img src="https://latex.codecogs.com/svg.latex?{\frac{k}{\alpha}}"/>
</p>

Consider a frequency from a selected set of band ![](https://latex.codecogs.com/svg.latex?{f_i\in\\{f_1,f_2,...,f_n\\}}), the amplitude of sine function at the time is <br/>
<p align="center">
<img src="https://latex.codecogs.com/svg.latex?{{\sin}2\pi{f_i}\frac{k}{\alpha}={\sin}\omega_i\frac{k}{\alpha}}"/>
</p>

Shown above indicates that band frequency can be first converted into angular freqency 
![](https://latex.codecogs.com/svg.latex?{\omega})
to reduce computation overhead. 
