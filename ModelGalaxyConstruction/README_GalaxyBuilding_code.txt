This folder contains various Metropolis Algorithm scripts for mock catalogue creation. Run times scale as the product of the number of stars and the number of Monte Carlo steps taken for each star in the Monte Carlo (MC) warming process. The Mersenne Twister is used for pseudorandom number generation, which must be downloaded as well, and included in the directory with the script (MersenneTwister.h). 

Approximately speaking, there needs to be enough MC steps such that:

N*step sizedistribution's characteristic spread>1

-The gaussian galaxy program is a toy model for testing  purposes.

-The erkal et al program builds a galaxy with a distorted halo, as in erkal et al.  Parameters from erkal et al are hard-coded into the program, and can  be changed in the script if desired.

-The chebyshev program takes an additional input to build the vertical distribution of stars.

Regardless of the particular model used, all scripts have the same general structure:
int main(int argC, char** argv)
{
 
  char* SEED = argv[1];
  string fileName = argv[2];

One can read off the arguments by looking at the argv[k] assignments at the top of main(). 


TO RUN:
g++ -O3 3d_args_BuildGalaxy.cpp
./a.out [pick an integer seed here, e.g. 433] [filename to write mock galaxy too]


FILE INPUTS:
Just the empty file to write a galaxy to.  If a file does not already exist, it will be created as named: [filename to write mock galaxy too]


OTHER INPUTS:
Many CRITICALLY IMPORTANT inputs reside in the code itself. All of the distribution function parameters are hard-coded at the top of the script. Monte Carlo step sizes, number of steps, and number of stars are all customizable within the script, right below the distribution function parameters.

Further, the distributionFunction() function has some easy to make cuts imposed within it. These cuts are hardcoded and will need to be changed by hand.

Next, in the main() function, the first for loop contains the initial conditions for the stars, before the stepping process. This initial point HAS to be within the volume you are simulating.

Finally, the phi coordinate for the Mock galaxy generation is drawn from a uniform distribution instead of a MC warming process.  I should have done this differently, but for now the phi cuts are hard-coded in:

float Phitemp = 179.5 + 0.5 * (2*ran1.rand() - 1); //phi 179-180  


OUTPUT:
The file specified in the command line to catch the data will be the output.


__________________


Not all cuts will be implemented in the MC process, as it is not practical.  As such, a script is also included to apply further cuts. 

To run:
g++ cutMockGalaxy.cpp
./a.out 

(the input and output files are specified within the code, should probably be updated to accept them as args)
