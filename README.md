# Convected Dubins C++ Library

Convected Dubins paths are time-optimal paths planned in a known, steady, uniform wind. The algorithm is based on the paper by L. Techy and C. Woolsey referenced below. 

<p align="center"> 
<img src="http://arturwolek.com/img/dubins.png" width="300">
</p>

## Dependencies:
- Ubuntu 16.04
- cmake
- (Optional) Octave/MATLAB
  for generating .m files for plotting

## Build instructions:
- Run `./build.sh` to create a library in `./lib` and executables in `./bin`

## Usage instructions:
- See test programs in `./programs` for examples of how to incorporate into your own code
- `clean.sh` removes all compiled files, build files, and temporary files 

## References:

[1] Techy, L., & Woolsey, C. A. (2009). Minimum-time path planning for unmanned aerial vehicles in steady uniform winds. Journal of Guidance, Control, and Dynamics, 32(6), 1736–1746. 
https://doi.org/10.2514/1.44580

## Contact:
Artur Wolek, wolek@umd.edu
