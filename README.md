# Cequal - Terminal DJ

<img src="img/dj.png" width=80%>

Cequal stands for: C-based equalizer for educational purpose. 

This equalizer demonstrates a typical conversion between time and frequency -domain, often applied to audio processing. 

It can amplify or decrease volumn in the selected bands. 

Better performance and higher spectrum resolution may be supported later. 

```diff
- Warning: gently adjust volumn especially in low frequncy 
```

## system prerequisite 

* pthread
* g++ 
* Open MP
* Open AL
* ncurses

## compile
```shell
sh compile.sh
```

## run 

```shell
./cequal <YOUR_WAV_FILE> [OUTPUT.WAV (optional)]
```
Output file by default is named *tmp.wav*.

## to do 

* better performance for higher spectrum resolution (to solve frequnecy leakage)
* refine ncurses refreshing synchronization 

This document is exported from: [https://github.com/superpi15/cequal/](https://github.com/superpi15/cequal/)


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


### FFT 

Define 
![](https://latex.codecogs.com/svg.latex?{W_N=e^{-j(2\pi/N)}}), DFT has the form 

<p align="center">
<img src="https://latex.codecogs.com/svg.latex?{X(k)=\sum^{N-1}_{n=0}x(n)W^{kn}_{N}}"/>
</p>

Partition into even and odd part 

<p align="center">
<img src="https://latex.codecogs.com/svg.latex?{X(k)=\sum^{N/2-1}_{r=0,n=2r}x(n)W^{kn}_{N}+\sum^{N/2-1}_{r=0,n=2r+1}x(n)W^{kn}_{N}}"/>
</p>

Using the identiry 
![](https://latex.codecogs.com/svg.latex?{W^2_N=(e^{-j(2\pi/N)})^2=e^{-j(2\pi/(N/2))}=W_{N/2}})

<p align="center">
<img src="https://latex.codecogs.com/svg.latex?{X(k)=\sum^{N/2-1}_{r=0}x(2r)W^{kr}_{N/2}+W^k_{N}\sum^{N/2-1}_{r=0}x(2r+1)W^{kr}_{N/2}}"/>
</p>

Observe the elements crossing half window size

<p align="center">
<img src="https://latex.codecogs.com/svg.latex?{X(k+N/2)=\sum^{N/2-1}_{r=0}x(2r)W^{kr}_{N/2}-W^k_{N}\sum^{N/2-1}_{r=0}x(2r+1)W^{kr}_{N/2}}"/>
</p>

Observe K+1-th term 

<p align="center">
<img src="https://latex.codecogs.com/svg.latex?{X(k+1)=\sum^{N/2-1}_{r=0}x(2r)W^{(k+1)r}_{N/2}+W^{(k+1)}_{N}\sum^{N/2-1}_{r=0}x(2r+1)W^{(k+1)r}_{N/2}}"/>
</p>

Observe K+(N/4)-th term 

<p align="center">
<img src="https://latex.codecogs.com/svg.latex?{X(k+N/4)=\sum^{N/2-1}_{r=0}x(2r)W^{(k+N/4)r}_{N/2}+W^{(k+N/4)}_{N}\sum^{N/2-1}_{r=0}x(2r+1)W^{(k+N/4)r}_{N/2}}"/>
</p>


