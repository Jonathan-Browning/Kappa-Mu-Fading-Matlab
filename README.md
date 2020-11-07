# Kappa-Mu-Fading-Matlab
The \kappa-\mu fading model implemented in Matlab and was created using Matlab 2018a.
Plots the theoretical and simulated, envelope porbability density functions (PDF).

Further details of this model can be found in Yacoub's paper: 
"The \kappa-\mu distribution and the \eta-\mu distribution".

Run main.m to start the GUI if Matlab is already installed.
Alternatively if Matlab isn't installed, can run the installer from the build folder, which requires an internet connection to download the required files.

The input \kappa accepts values in the range 0 to 50.
The input \mu accepts integer values in the range 1 to 10.
The input \hat{r} accepts values in the range 0.5 to 2.5.

When running the program the intial window appears:

![ScreenShot](https://raw.github.com/Jonathan-Browning/Kappa-Mu-Fading-Matlab/main/docs/window.png)

Entering values for the \kappa, \mu and the root mean sqaure of the signal \hat{r}:

![ScreenShot](https://raw.github.com/Jonathan-Browning/Kappa-Mu-Fading-Matlab/main/docs/inputs.png)

The theoretical evenlope PDF is plotted to compare with the simulation and gives the execution time for the theoretical calculation and simulation together:

![ScreenShot](https://raw.github.com/Jonathan-Browning/Kappa-Mu-Fading-Matlab/main/docs/results.png)