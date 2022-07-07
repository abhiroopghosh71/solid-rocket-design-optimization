# Solid rocket motor design using multi-objective optimization
Multi-objective optimization of a solid rocket motor to fit a pre-specified thrust profile and minimize insulation.

<h2>Note to macOS and Windows users</h4>
The code currently may not work properly on MAC-OS. One known issue is when you run ```make``` it might say:

```xcrun: error: invalid active developer path (/Library/Developer/CommandLineTools), missing xcrun at: /Library/Developer/CommandLineTools/usr/bin/xcrun```

Try to follow the instructions given in https://apple.stackexchange.com/questions/254380/why-am-i-getting-an-invalid-active-developer-path-when-attempting-to-use-git-a but its not guaranteed to fix the issue.

The code should run properly on Windows but it has not been extensively tested. Please contact the developer for any issues faced.

<h2>Running unit tests to ensure everything is in order:</h4>
1. Clone the repo:

    ```git clone https://github.com/abhiroopghosh71/DARPA-TRADES-CP3.git```
    
2. Change into the pymoo_optimization directory ```cd DARPA-TRADES-CP3/pymoo_optimization```

3. To install the dependencies use ```pip install requirements.txt```. If you are using Anaconda, it is recommended to create a new virtual environment using ```conda create --name <envname> --file requirements.txt```

4. Run ```make clean```

5. Run ```make```

6. Run ```python -m unittest tests/unit/test_rocket_evaluator.py tests/unit/test_pymoo_optimization.py```

<h2>A sample case for the standalone rocket burn simulator</h4>
From the ```pymoo_optimization``` directory, run ```python rocket_eval_demo.py```. The file has a predefined input in the ```__main__``` block.

Expected output:

```
Thrust Reward = 321114.5209178659, Simultaeinity Reward = 428436.2400948988
Stop Code 0
Thrust Profile
[7238.34365998 6355.11604302 6428.12892106 5944.83192139 5765.10659394
 5835.07203787 5593.34972486 5182.88571753 4967.62473903 4864.53466283
 4902.6610781  4636.60217346 4130.13557409 4075.20247541 3855.60918429
 3647.39638416 3461.04689244 3364.87654633 3114.91404598 3134.99595683
 3101.27984383 5437.7673139  5671.36289391    0.            0.
    0.            0.            0.            0.            0.
    0.            0.            0.            0.            0.
    0.            0.            0.            0.            0.        ]
Pressure Profile
[3289372.96273456 2888013.80011734 2921192.54931581 2701571.15686828
 2619899.78978036 2651693.7140755  2541849.47542816 2355325.08793923
 2257505.49487392 2210658.96758795 2227984.49705745 2107081.14539022
 1876930.94818628 1851968.07105024 1752179.77441134 1657563.03745279
 1572881.49612103 1529179.45506106 1415590.683716   1424716.37087569
 1409394.9852884  2471149.2108745  2577300.47311134       0.
       0.               0.               0.               0.
       0.               0.               0.               0.
       0.               0.               0.               0.
       0.               0.               0.               0.        ]
Residual fuel thickness for each segment
[2.94463669e-04 7.80797052e-05 2.30986229e-04 2.31323498e-03
 2.16117296e-04 3.08452273e-04 6.36774075e-04 1.25373874e-03
 6.28214145e-05 1.60950331e-04 1.22376606e-04 2.88977000e-04
 1.05326938e-04]
Ray burn depths for each star segment
[ 2  2  3  4 11 15  1  2  8 15 15 15  1  5 10 15 10 11]
Segment star burning or not
[0 0 0 0 0 0 0 0 0 0 0 0 0]
Time of circularization for each segment
[0.    0.    3.475 4.275 5.1   0.    0.    0.    0.    0.    0.    0.
 0.   ]
```

<h2>Running the optimization</h4>
The code uses NSGA-II [1] optimization algorithm to find the solid fuel configuration to match a given thrust profile.

1. Change into the pymoo_optimization directory.

    ```cd DARPA-TRADES-CP3/pymoo_optimization```

2. Run ```make clean```

3. Run ```make```

4. Run ```python optimize.py```. Add command line parameters if necessary.

The optimization will run in the default setting for 200 generations with a 100 population size and print out the Pareto Front in the end.

Expected output:

```
[[3.08867363e+07 1.23672841e-03]
 [2.97021486e+07 1.24345285e-03]
 [2.77866892e+07 1.36024597e-03]
 [3.09107749e+07 8.80772066e-04]
 [2.54478826e+07 1.45007993e-03]
 [2.25408356e+07 3.55183206e-03]
 [2.33277370e+07 2.67428477e-03]
 [2.33702714e+07 1.96810693e-03]
 [4.39457197e+07 7.59284727e-04]]
```

<h2>Important command line optimization parameters</h4>

1. ```--seed```: Sets the seed for the random number generator.

2. ```--target-thrust```: Thrust profile to be matched by the optimizer.

3. ```--ngen```: Number of generations of NSGA-II.

4. ```--popsize```: Population size of NSGA-II.

For example,
```python optimize.py --ngen 100 --popsize 200 --target-thrust baseline --seed 1234```

This runs the optimization for 100 generations and a population size of 200 with a seed of 1234.

# README for the Solid Fuel Rocket Simulator
<h4>Authored by Erik Goodman, Prof. at Michigan State University</h4>

This simulator is based on the rocket burn model described in the Challenge Problem 4 statement from
DARPA entitled, "TRADES Solid and Hybrid motor study statement", which originated with Tim Kibbey
at NASA's Jacobs Space Exploration Group.  

Goodman's rocket burn simulator, was written to simulate the burning of a model rocket of 
radius 7.62cm and length of cylindrical part about 53cm. The simulator is written in 'C' and can be 
called from a routine in another language if desired. It was written to provide an optimizer with a 
model that has a great deal of design flexibility while being organized such that it does not require 
a great deal of computer time to simulate, enabling many designs to be evaluated.

<h2>Segments and layers</h2>
It simulates a propellant load divided into SEGMENTS, each of which includes fuel distributed in LAYERS.
Layers are organized into a user-selectable number of NZONES zones, each of which represents 1/NZONES of
the web thickness. Within each zone, the thickness of each of the NLAYERSPERZONE layers can be optimized, 
within certain limits based on manufacturability. A special integer encoding of the layer thicknesses 
within any zone is used to reduce the size of the design space and ease the work of the optimizer. The
propellant in each layer is one of the eleven choices (0-10) as specifiec in the Kibbey statement of the
challenge problem.

<h2>STARGRAIN</h2>
IF the propellant in the cylindrical part of the rocket is distributed only in cylindrical layers--i.e., 
pragma "STARGRAIN" is undefined--all parts of the initial burn surface are equidistant from the 
nearest portion of the shell, as illustrated in the accompanying PowerPoint slides. The cylindrical 
segments can each be burned by tracking only a single radius (current depth of burn). This greatly 
reduces the simulation workload, allowing many more rockets to be simulated than would an arbitrary
distribution of propellants, so allowing a better job of optimization. 

<h2>Elliptical dome and corner segments</h2>
The "central dome" segment is similarly layered to the cylinder segments, and the initial burn surface 
is offset by a uniform amount (along surface normals) from the elliptical dome interior surface. Thus 
it burns with an INCREASING surface area, growing with the burn depth as illustrated in the PowerPoints.

The "corners" of the dome are irregular, as some of their circumference is defined by the ellipse, and 
some by the cylinder wall. They are modeled with a user-definable number of "pie slices" (each a solid 
of revolution in 3-D, of course). Each corner segment, in turn, is irregular, as its bounding rays are
of different lengths. In order to enable simultaneous burnout to be possible with propellant that is
distributed in a relatively simple and consistent way, a short final portion of each corner segment is
burned in a special distribution of two particular propellants, one faster burning than the other,
with a line defining the distribution boundary, so that the longer side of the corner slice can "catch
up" and burn out at the same time as the shorter side.

<h2>Nozzle segment</h2>
The nozzle segment, bounded by one cylindrical segment and the given insulation geometry, burns similarly
to the corner segments, Because it faces an irregular boundary distance, it is burned in a (fixed by me)
number of rays (7), when in the final ZONE of the segment, allowing an adequate representation of the burn. 
(See PowerPoints.) In all but the last zone, it burns as a single radius with a growing length, similarly 
to the corner segments. In three dimensions, the burn surface is a portion of a torus until the ray
burning commences.

<h2>Stargrain evolution and circularization</h2>
Because it is desirable to provide high initial thrust and sustained burning, the model allows simulation
of star or star-like cross-sections in the initial grain profile. That is allowed only in the cylindrical
segments of the rocket. It is not allowed in the short cylindrical segments adjacent to the nozzle and
the dome, in order that the burn geometry of the star portions not interact with the different geometries
of the nozzle and dome. But the remainder of the cylinder may be divided into a user-selectable number
of segments, and each segment may have a user-defined or an evolved (optimized) initial profile. In cross-
section, the segment is seen as having K-way symmetry and anti-symmetry, as illustrated in the 
PowerPoints. If straight lines are used, then choosing NSYMS = 12, for example, produces a 6-pointed star
in the cross-section, with the star commencing at what would otherwise be the first layer of the 
cylindrical burning segment. That is, the star shapes are all interior to the cylindrical NLAYERS
layers of propellant that is cylindrically distributed. The propellant in the stars is divided into
two further radial layers, with the propellant types evolvable. The depth of the stars or the profiles
of the initial burn surfaces, if not pure stars, is evolvable, and the initial profile is represented
as depths along rays emanating from the rocket's centerline. Only a 1/NSYMS portion of the burn need
be tracked by these rays, and all burns are tracked along these rays until all of them have burned into
the "normal" cylinder layer 0. At that time, a "circularization" process is begun, in order to enable
simultaneous burnout of all rays. That process will place faster-burning propellant in the region where
the initially furthest-from-shell rays are burning, and slower-burning propellant in the area where the
initially closest-to-shell ray is burning. As each subsequent ray catches up with the deeper, slow-burning
rays, it begins burning slow propellant, as well. This allows the entire segment eventually to begin 
burning in a circular arc (or cylinder, in 3-D), making simultaneous burnout possible and greatly speeding
the simulation speed, as the segment may again be represented by a single burn depth. (Please see the
PowerPoints for visualization of this process.)

<h2>Calling the simulator</h2>
For use with NSGA-II or other direct-calling optimizers, the simulator is designed to be called with the 
optimizable design variables defining a single rocket, and to return the resulting burn outcome and a 
code for burn termination. Outputs include the thrust profile produced at 0.5-second intervals, the 
pressures at those times, and the thickness of propellant unburned in each segment when the rocket's burn 
is terminated. These outputs allow the calling (optimizer) program to determine how well the rocket is 
meeting the objectives (closeness to desired thrust profile and simultaneity of burnout) and constraints 
(on pressures).  When the rocket is called from HEEDS, the parameters are instead passed via reading and 
writing of text files, provided for in the simulator by uncommenting the necessary write and read 
statements.

Burn terminates for the following reasons, represented by "StopCodes":
0) One of the segments burned to the shell of the rocket, halting the burn calculations. Remaining
   propellant depths in all segments are reported.
1) The burn terminated because the pressure fell below the minimum allowable pressure, PMIN.
2) The burn terminated because the pressure rose above the maximum allowable pressure, PMAX.
3) "Normal termination" -- the rocket ran out of the allotted burn time, currently set at 15 seconds.
4) Strange termination -- got a NaN from newton solver, for example, from a negative pressure call.

<h2>Configuring the simulator</h2>
The simulator must first be configured by filling in or commenting out some #define pragmas. They
tell the simulator how the user wants to configure the rocket--how many segments of each type are to
be used, how many zones, how many layers per zone, whether or not star segments will be used, and if
so, whether straight star edges or evolved shapes are to be used, etc. All of those #defines are
in file approcket.h. Here are the elements to be defined and their meanings (typical values used 
during development are shown, then any range information). Note that OTHER #defines MAY NOT be changed 
by the user as their values are "hard-coded" into the logic of the simulator. The user-settable 
#defines in file approcket.h are:

- ```#define NZONES 4```   (may be 2-6), the number of zones into which the web is divided EQUALLY (except for
                                  the irregular last zones of nozzle and corner segments)

- ```#define NLAYERSPERZONE  5```  (must be 2 or greater), the number of layers in each zone, each of which may
                                  have a different thickness and propellant type
- ```#define NCYLSEGS 6```  (must be 2 or more), the number of cylindrical segments INCLUDING one for the
                                  central dome segment, which, although not cylindrical, uses the same
                                  data structures
- ```#define NCORNERSEGS 6```  (1 or more, but should be 4 or more for any accuracy), the number of segments
                                  into which the region between the central dome and the cylindrical 
                                  segments is divided
- ```#define INITBURNRADIUS  0.03```  (may be 0.024 to 0.035 without requiring changing of the code calculating
                               the initial position of the corner between dome and cylinder segments)
                                  Defines the radius at which the web begins in pure cylinder segments,
                                  and implicitly, the dome thickness and the radius at which the deepest
                                  portions of any "star" shapes are located. 
- ```#define STARGRAIN```       (defined or commented out) asserts that the user wants the cylinder segments, 
                                  except those adjacent to nozzle and corner segments, to have an interior
                                  "star-shaped" portion of additional propellant, in each of NSYMS 
                                  fractions of the 2 pi radians of the cross-sections of the cylinder
                                  segments (see PowerPoint for illustration).
- ```#define NSYMS```		(MUST be defined if STARGRAIN is, and otherwise, will not be), tells the number of
                                  symmetrical (or antisymmetrical) regions into which the cylinder grain
                                  cross-section is to be divided, until such time as the rays have been
                                  circularized for the segment (again, see PowerPoints).
- ```#define INNERSTARPOINT  0.003``` (0.001 to INITBURNRADIUS), the depth (from centerline) of the "points" of a
                                  straight-sided star shape if not EVOLVINGSTARSHAPE, or the minimum
                                  allowable depth of any of the initial points if EVOLVINGSTARSHAPE.
- ```#define NRAYS	5```    (MUST be defined if STARGRAIN is, and otherwise, will not be), the number of "central" 
                                  rays along which the star burns are to be tracked. Two additional rays 
                                  (boundaries of the symmetry region) are ALWAYS also present, and are 
                                  NOT included in the count NRAYS. That is, NRAYS 5 results in tracking 
                                  along 7 total rays. Once the automatic differential choice of 
                                  propellants in early segments has allowed all ray burn depths to match, 
                                  these rays are not used for the remainder of the cylinder burn, which is 
                                  then tracked along a single burn depth radius.
- ```#define EVOLVINGSTARSHAPE 1```  (can be defined only if STARGRAIN is, but may be commented out even then), 
                                  flag that determines whether the "star shape" is truly a straight-sided
                                  star with NSYMS/2 points, or a much freer shape, determined by 
                                  optimization of the initial burn (propellant)  depths of each of the 
                                  NRAYS+2 rays. In order to reduce the size of the search space and to
                                  assure that unreasonable profiles are not introduced, the initial
                                  profiles are defined by a series of integer codes for each STARGRAIN
                                  segment, and the simulator automatically calculates the corresponding 
                                  initial depths along the rays. 
- ```#define NRAYDEPTHS  16```   (Can be 2 or greater), the number of depths into which the region from the 
                                  INNERSTARPOINT to the INITBURNRADIUS will be divided as possible initial
                                  depths along any of the rays defining the star shape.


In file approcket.c (near the top):

- ```#define NPROPELLTYPES 11``` (fixed by NASA for this problem, but reconfigurable with changing of the
                           array PropellChoices[]), number of reference burn rates (RDotRefs) available
                           to the optimizer.
- ```#define MAXTHRUST 7532.6``` (Must be set to the maximum thrust called for in the TargetThrusts array.)

- ```double MaxPressure = 3447000;``` If PHASE2 is not defined, it is set to 3447000Pa.
                           If optimizing throat area, value is calculate from the MAXTHRUST
                              specified above (from the thrust profile), using the value of
                              Isp determined from the throat area, etc. 

- ```#define MINPRESSURE 1379000``` Specified by NASA. Is automatically relaxed when burn time exceeds 10 
                               seconds, as required in order to continue meeting 25g thrust targets for
                               phase 1.

- ```double DeltaTime = 0.5```   (seconds, user may change to a number that is a multiple of 0.5 seconds, 
                             so 1, 1.5, 2, etc., are all fine, but 0.2, 0.3, for example, are not)

<h2>Calling the simulator</h2>
The 'C' language call from the optimizer that passes a rocket to approcket.c to be simulated is provided 
below, and can be seen in the accompanying test program, test_rocket.c:
```
objectivefunction(&tReasonStopped, &tThrustReward, &tSimultaneityReward, tTimes, tThrustProfile, tPressureProfile,
                          tTimesFine, tThrustProfileFine, tPressureProfileFine,
                          tBurnDistRemaining, tPropellantToEval, tLayerStartRadius, tRayBurnDepthInput, tPreCylStarRDotRefInput,
                          tRDotRefCatchupInput, tThroatAreaInput, &tDeltaV, tMaxThrust, tRayBurnDepths, tSegStarBurnStatus,
                          tTimestarBurnFinish, tSegResidualPerTimestepCoarse, tSegResidualPerTimestepFine, tTargetThrusts,
                          rayDepthFlag, tSegBurnAreaPerTimestepCoarse, tSegBurnAreaPerTimestepFine,
                          tBurnLayerPerTimestepCoarse, tBurnLayerPerTimestepFine);
```

The data structures being passed are shown in file test_rocket.c, for the particular case shown here:
- ```#define NZONES 4           ```
- ```#define NLAYERSPERZONE 5   ```
- ```#define NCYLSEGS 6         ```
- ```#define NCORNERSEGS 6      ```
- ```#define STARGRAIN          ```
- ```#define INNERSTARPOINT 0.03```
- ```#define EVOLVINGSTARSHAPE  ```
- ```#define NSYMS 12           ```
- ```#define NRAYS 5            ```
- ```#define NPROPELLTYPES 11   ```
- ```#define MAXPRESSURE 3447010```
- ```#define MINPRESSURE 1379000```
- ```double  DeltaTimeCoarse = 0.5;```


<h2>Compiling and running the code</h2>
Passing the Python unit tests at the beginning of the README ensures that the rocket 
simulator is also correctly running. To test the simulator separately as a C program, 
compile the test program that calls the simulator with a single rocket specification. 
On a Linux system,

```gcc -g ./test_rocket.c -lm```

should be all that's needed to compile it. It produces an executable a.out that,
when run, prints some information to standard output about the rocket.

Please report issues to me, ghoshab1@msu.edu

<h2>References:</h4>
1. A. Ghosh, E. Goodman, K. Deb, R. Averill, A. Diaz, "A Large-scale Bi-objective Optimization of Solid Rocket Motors Using Innovization," _2020 IEEE Congress on Evolutionary Computation (CEC)_, 1–8. https://doi.org/10.1109/CEC48606.2020.9185861
2. K. Deb, A. Pratap, S. Agarwal and T. Meyarivan, "A fast and elitist multiobjective genetic algorithm: NSGA-II," in _IEEE Transactions on Evolutionary Computation_, vol. 6, no. 2, pp. 182-197, April 2002. https://doi.org/10.1109/4235.996017

