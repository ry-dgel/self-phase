# Self-Phase Modulation
## Code for simulation to be translated from matlab -> julia
Done by Rigel and Joseph as part of undergraduate thesis.

```
julia -p <number of sims to run at once/cores to use> multi-thread.jl <Parameter File> <Number of times to save data>
```

view params.yaml for example of parameter file
## Todo
* ~~Translate Main Loop~~
* ~~Translate "PPT.m"~~
* ~~Turn chunks of main loop into functions~~
* Optimize
* Test
* ~~Use plan fft/ifft to speed up all the transforms~~
* ~~Multi-Thread drifting~~
