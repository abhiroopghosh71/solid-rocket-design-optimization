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
propellant in each layer is one of eleven choices (0-10).

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

- ```double DeltaTimeCoarse = 0.5```   (seconds, user may change to a number that is a multiple of 0.5 seconds, 
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