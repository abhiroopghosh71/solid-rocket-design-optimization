/* PUT GLOBAL DECLARATIONS AND DEFINITIONS TO BE USED THROUGHOUT THIS FILE HERE. */
#define NZONES 4
// NZONES must be between 2 and 6 (more than 6 breaks the nozzle segment burn, although it could be 
//    fixed not to if necessary; just easier this way)!!!
#define NLAYERSPERZONE 5
// NLAYERSPERZONE must be 2 or greater!!!
#define NLAYERS (NZONES * NLAYERSPERZONE)
// So, NLAYERS will be set to 4 or greater.

// NOW inserting NOZZLE segment AFTER all other segs (i.e., after CORNER segs. 
/* NOTE:  NCYLSEGS includes all cylindrical segments and the central dome cap, NOT the corner segments in the dome. Insulation is not a seg. */
#define NCYLSEGS 6
// Must be 2 or more:  [0] is for central dome, remaining are for true cylinders, but include the seg next to dome/corners
// and the seg next to nozzle (but those could be the same segment).
#define NCORNERSEGS 6
// Can be any number, but should certainly be 4 or greater for any accuracy.
// NOZZLE is the final segment
// AUTOMATIC nozzle rays are two bdaries, shortest ray perp to straight edge, and ray to center of external arc.
// Six other rays are interpolated, half of them betw vert bdary and cent arc ray, others betw.
// central arc and shortest ray.
#define NOZZLESEG (NCYLSEGS+NCORNERSEGS)
// The insulation burn has been added, but IT IS NOT A SEGMENT! 
#define NTOTALSEGS (NCYLSEGS+NCORNERSEGS +1)
// So NTOTALSEGS is usually going to be at least 8.
#define NNOZRAYS 9 // This cannot be changed without doing so in the code, for now. */ 
                   // Specific rays are placed in key places relative to nozzle seg geometry.

#define NSEGINTERFACES (NCYLSEGS+NCORNERSEGS)
// Should be one for each cylinder/cylinder boundary, one for cyl/corner boundary, and among 
// all dome segments, including central (cyl) segment, and between nozzle seg and CylSeg 1 

#define CHAMBERRADIUS 0.0762  // Outer radius of the rocket motor

#define INITBURNRADIUS 0.0340
            /* You can vary if desired. Reasonable values maybe .030 to .038 */

// Uncomment if want to calculate DeltaV as another objective
//#define PHASE2

#define STARGRAIN

#ifdef STARGRAIN
#define INNERSTARPOINT 0.003  /* was original value for evolving and straight runs */
/* Number of symmetrical or antisymmetrical pie slices into which the circle is divided. */
#define NSYMS 12
/* This NRAYS is in addition (and interpolating between) 2 rays automatically defined along the 2 edges */
#define NRAYS 5
#define EVOLVINGSTARSHAPE 1
#define NRAYDEPTHS 16
#endif  /* of STARGRAIN */

#define NPROPELLTYPES 11  /* Coded from 0 to 10 */

#define  DeltaTime 0.025 /* seconds */
#define  DeltaTimeCoarse 0.5 /* seconds */

#define NumTimePointsCoarse 40 // > (int) (MaxBurnTime / DeltaTimeCoarse) + 1  = 31
#define NumTimePointsFine 1000 // > (int) (MaxBurnTime / DeltaTime) + 1 = 601

void objectivefunction(int *ReasonStopped, double *lThrustReward, double *lSimultaneityReward
		  , double *lTimes, double *lThrustProfile, double *lmPressureProfile
		  , double *lTimesFine, double *lThrustProfileFine, double *lmPressureProfileFine
		  , double *lBurnDistRemaining, int *lPropellantToEval
		  , double *lLayerStartRadius, int *lRayBurnDepthInput
		  , int *lPreCylStarRDotRefInput, int *lRDotRefCatchupInput
		  , double lThroatAreaInput, double *lDeltaV, double lMaxThrust, int *lRayBurnDepths
		  , int *lStarBurningThisSegNow, double *lTimeStarBurnFinish
		  , double lSegResidualPerTimeStep[][NumTimePointsCoarse]
		  , double lSegResidualPerTimeStepFine[][NumTimePointsFine], double lTargetThrusts[NumTimePointsFine], int lRayDepthFlag
		  , double lSegBurnAreaPerTimeStepCoarse[][NumTimePointsCoarse]
		  , double lSegBurnAreaPerTimeStepFine[][NumTimePointsFine]
		  , double lBurnLayerPerTimestepCoarse[][NumTimePointsCoarse]
		  , double lBurnLayerPerTimestepFine[][NumTimePointsFine]);
