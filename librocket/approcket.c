/*
    TITLE: Solid-fuel rocket simulator.
    AUTHOR: Erik Goodman, Michigan State University.
    VERSION: Alpha Release 1.08, Sept. 10, 2019.
             (Refactored by Abhiroop Ghosh on July 7, 2022)
    FILES:
        approcket.c
        approcket.h
    DESCRIPTION:
        Erik Goodman's 'C'-language simulator function of a solid-fuel rocket, for
        calling from an optimizer such as NSGA-II or HEEDS, begun 10/9/18.
        See accompanying README for code description and release notes.

        Simulator is designed to work with optimization problems where the thrust
        constraints are imposed only on half-second boundaries, as NASA originally
        specified, rather than at every timestep, as would be needed for a
        real rocket. This simulator is an early version which does not include the
        special propellant regions at segment boundaries needed to prevent premature
        burnout, which renders that code's burn calculations inaccurate at the
        segment boundaries. Future versions (under development since 2019) will
        include special "compatibilization regions" to satisfy the thrust constraints
        at every 25ms timestep.

    Based in part on:
        Plato_RocketDesignApp.cpp
        Created on: Sep 28, 2018
        Plato_SimpleRocket.hpp
        Created on: Aug 29, 2018
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "approcket.h"

FILE *InnovizDatafp;
FILE *Evaluationfp;
FILE *PropellantToEvalfp;
FILE *ThicknessesToEvalfp;
FILE *Graphfp;
FILE *Animatfp;
FILE *Restartfp;

int WantWriteRestartData;

FILE *outfp;
double Fitness = 0.;
int gen = 0;
long neval;
float pmutation;
int ncycles;
int cycle;
int ReasonStopped;
double ThrustReward, SimultaneityReward, StdDevReward;
double lThrustReward;
double lSimultaneityReward;
#define G0 (9.80665)
#define PI (3.141592653)

#ifdef STARGRAIN
/* Int array into which to decode the ints comprising the chromosome for calling
 * chromtointarray() */
int ChromValues[NLAYERS * 2 * NTOTALSEGS + (NCYLSEGS - 3) * (NRAYS + 3)];
/* Copy of ChromValues, almost:  "jammed" values will be inserted, ChromValues
 * properly distributed around them */
int PhenoValues[NLAYERS * 2 * NTOTALSEGS + (NCYLSEGS - 3) * (NRAYS + 3) + 100];
#else
/* Int array into which to decode the ints comprising the chromosome for calling
 * chromtointarray() */
int ChromValues[NLAYERS * 2 * NTOTALSEGS +
                10]; /* Guard against index getting off by 1 */
/* Copy of ChromValues, almost:  "jammed" values will be inserted, ChromValues
 * properly distributed around them */
int PhenoValues[NLAYERS * 2 * NTOTALSEGS + 100];
#endif
/* Then two arrays into which to store the values for manipulation */
int PropellType[NLAYERS + 1][NTOTALSEGS];
int LayerThicknessInts[NLAYERS][NTOTALSEGS];
double LayerThicknesses[NLAYERS][NTOTALSEGS];
double InitBurnRadius[NTOTALSEGS]; /* Can vary if desired. For cylinder only, I
                                      set to 3.077 below, to */
/* match initial thrust with as BIG a propellant load as possible for long burn.
 */
/* Must be set to same value for ALL segments except any cylinders NOT adjacent
 * to */
/* elliptical dome. InitBurnRadius from all nozzle segments must match that of
 */
/* cylinder seg [NCYLSEGS-1], with which they share an initial burn point. */
double RocketCylRadius = 0.0762;
double FinalBurnRadius[NTOTALSEGS];
double WebThickness[NTOTALSEGS];
double BdaryLen[NCORNERSEGS +
                1]; /* Only for corner segments, calculated internally */
double theta[NCORNERSEGS];

double ZoneThickness[NZONES][NTOTALSEGS];
double DistanceLeftToBurn[NTOTALSEGS];
double LayerStartRadius[(NLAYERS + 10)]
                       [NTOTALSEGS]; /* Later will deal with in detail */
double RDotRef[(NLAYERS + 9)][NTOTALSEGS];

double PropellChoices[NPROPELLTYPES + 1] = {0.00254, 0.00305, 0.00363, 0.00434, 0.00521,
                                            0.00622, 0.00744, 0.00892, .01064,  .01275,
                                            .01524,  0.}; /* Added so if still burning at shell
                                                             boundary, can give no fuel! */

/* Extended beyond 10 seconds in case can burn longer, same geometric rate of
 * decrease */
double MaxThrust = 7532.6;
double TargetThrusts[NumTimePointsFine];

double MaxPressure = 3447000.;

//#define MINPRESSURE 1379000.
#define MINPRESSURE 1000.

/* DECLARATIONS and initializations for parameters */
double RDot[NTOTALSEGS]; /* burn rate in seg, AFTER adjustment for press not
                            being Pref. */
unsigned MaxNumNewtonItr = 100;

double SegmentLength[NTOTALSEGS];
double RefBurnRate[NTOTALSEGS];
int irecord;
int irecordFine;

/* SegmentRadius is for each segment to update as the burn progresses,
 * incrementing each step by RDot * DeltaTime */
double SegmentRadius[NTOTALSEGS];
double NextSegRadius[NTOTALSEGS];
double NextNozRayRadius[NNOZRAYS];
double RefPressure = 3.447e6; /* Pa */
double Alpha = 0.38;
double Alphaminus1;
double ThroatDiameter = 0.04064;       /* m */
double CharacteristicVelocity = 1555.; /* m/sec */
double PropellantDensity = 1740.;      /* kg/m^3 */
double AmbientPressure = 101.325;      /* Pa */

// Added by Abhiroop 24th Sept 2019
double SegResidualPerTimeStep[NTOTALSEGS][NumTimePointsCoarse];
double SegResidualPerTimeStepFine[NTOTALSEGS][NumTimePointsFine];

// Added by Abhiroop 27th Feb 2020
double SegBurnAreaPerTimeStepCoarse[NTOTALSEGS][NumTimePointsCoarse];
double SegBurnAreaPerTimeStepFine[NTOTALSEGS][NumTimePointsFine];

// Added by Abhiroop Mar 3, 2020
double BurnLayerPerTimestepCoarse[NTOTALSEGS][NumTimePointsCoarse];
double BurnLayerPerTimestepFine[NTOTALSEGS][NumTimePointsFine];

/*  NOW, extend TotalBurnTime so can burn longer than 10 seconds. Set arb.
 * at 12. for now */
double TotalBurnTime = 16.;     /* sec */
double NewtonTolerance = 1.e-8; /* Pa */
double ThroatArea = 0.00129717; /* m^2, Is overridden if being evolved. */
double Time = 0.;
double Thrust = 0.;
double TotalPressure;
double DeltaPressure; /* communicates with Pressure Violation() */
double InvPrefAlpha;  /* intermediate value for scaling pressure calculation */
_Bool Burning = 1;    /* True as long as burning time does not exceed
                         TotalBurnTime, now 10 seconds */

// Stores data at 0.5 second intervals
double Times[NumTimePointsCoarse];
double ThrustProfile[NumTimePointsCoarse];
double mPressureProfile[NumTimePointsCoarse];

double TimesFine[NumTimePointsFine];
double ThrustFine[NumTimePointsFine];
double mPressureFine[NumTimePointsFine];

double SegmentArea[NTOTALSEGS];
double BdaryAngle[NCORNERSEGS + 1];
double ArcLength[NCORNERSEGS];
double Px, Py, Qx, Qy;
double MajAx;
double SemiMajAx;
double Eccentricity = 3.;
double MinAx;
double SemiMinAx;
double rayangle;
double bigtheta;
double initialchord, finalchord, depth, chordlength;
double IntegralContribs;
double PropellVolume[NTOTALSEGS];
double PropellMass[NTOTALSEGS];
double PropellantTotalVolume;
double PropellantTotalMass;

int WantThisGuyOutputFile;
int WantToAnimate = 233223; /* Put the neval of any guy you want to animate */
int FirstAnimationCall;
int ActivePtr = 0;

#ifdef STARGRAIN
int StarFieldsActivated;
/* Allocating these vars for all CYLSEGS, to keep subscripts easy; some won't be
 * used */
double PtAx[NCYLSEGS], PtAy[NCYLSEGS], PtBx[NCYLSEGS], PtBy[NCYLSEGS],
    StarTheta[NCYLSEGS];
double RayAngle[NRAYS + 2][NCYLSEGS], BurnPtx[NRAYS + 2][NCYLSEGS];
double BurnPty[NRAYS + 2][NCYLSEGS], RayBurnDir[NRAYS + 2][NCYLSEGS],
    RayBurnDepth[NRAYS + 2][NCYLSEGS];
double PtBRadius[NCYLSEGS];
double RayInitDepthChoices[NRAYDEPTHS][NCYLSEGS];
double RayiToRayiminus1[NRAYS + 2][NCYLSEGS];
double SpecLayerOffset[NRAYS + 2][NCYLSEGS];
double RayRDotRef[NRAYS + 2][NCYLSEGS];
double RayRDot[NRAYS + 2][NCYLSEGS];
double RDotRefCatchup[NCYLSEGS];
double LatestCatchupDepthBurned[NRAYS + 2][NCYLSEGS];
double SumOfProducts[NCYLSEGS];
double PreCylStarRDotRef[2][NCYLSEGS];
double StarRayOuterInnerBdary[NCYLSEGS];
double FurthestRayDistTtoQ[NCYLSEGS];
int StarBurningThisSegNow[NTOTALSEGS];
double TimeStarBurnFinish[NTOTALSEGS];
double lSegResidualPerTimeStep[NTOTALSEGS][100];
int StarRayInitialPropellantType[2][NCYLSEGS];
int RayInLayer[NRAYS + 2][NCYLSEGS];
int FirstCompensatingLayerNo[NCYLSEGS];
int FarRayNo[NCYLSEGS];
int ShortestRayNo[NCYLSEGS];
int CompensatingNow[NCYLSEGS];
int NRAYSinPreCylLayer[NCYLSEGS];
double localoffsets[NRAYS + 2][NCYLSEGS];
int NCompensLayers[NCYLSEGS];
int NRaysToCatchUp[NCYLSEGS];
int StartCatchupLayer[NCYLSEGS];
int CatchingUpRay[NRAYS + 2][NCYLSEGS];
double StarRayArea[NRAYS + 2][NCYLSEGS];
double RayCylIntersectx[NRAYS + 2][NCYLSEGS];
double RayCylIntersecty[NRAYS + 2][NCYLSEGS];
#endif

double BurnDistRemaining[NTOTALSEGS];
double ThicknessesToEval[NLAYERS * NTOTALSEGS];
/* Inputs for propellants will come as int codes. */
int PropellantToEval[NLAYERS * NTOTALSEGS];
#ifdef STARGRAIN
int RayBurnDepthInput[(NCYLSEGS - 3) * (NRAYS + 2)];
int PreCylStarRDotRefInput[NCYLSEGS - 3];
int RDotRefCatchupInput[NCYLSEGS - 3];
int RayDepthHEEDS[(NCYLSEGS - 3) * (NRAYS + 1)];
#endif

double Reward;             /* to be used in tallying fitness */
double SimultaneityReward; /* so can watch separately */
double StdDevReward;
double StdDev; /* of remaining burn distances */

double localmax = 0.; /* Used in printing thrust file in calc_fitness() */
int innovcounter =
    0; /* counts new best guys printed to file for innovization use. */
int StopCode = 0;
#ifdef STARGRAIN
int FirstRayCallThisGuy[NCYLSEGS]; /* Flag to control initializing print */
#endif

void RecordThrustAndResidualData(int istep);
double newton(double TotalPressure, double ThroatArea);
double jacobian(double TotalPressure, double ThroatArea);
double residual(double TotalPressure, double ThroatArea);
double FindPropellantType(int SegNo);
double calc_fitness(void);
void UpdateInsulationBurn(void);
void CalcPropellantMass(void);

#ifdef STARGRAIN
void UpdateStarRayBurns(int SegNo);
void RayDebugPrint(int SegNo);
void SetupStarBurn(int SegNo);
void CalcRayBurnDirnsAndLengths(int SegNo);
void DefineInitStarProfile(int RayNo, int SegNo);
void CalcStarBurnIntegrContribs(int SegNo);
void CalcPerRayRDot(int SegNo);
void CalcStarRayBurnDirnsAndLengths(void);
void EndCatchUp(int j, int SegNo);
int FindShortestRayNo(int SegNo);
int FindFarRay(int SegNo);
void StarAreaCalc(void);
#endif

void decodechromosome(void);
void BurnedToShellFinish(void);
int BadPressure(void);
void PrintAllSegs(void);

double bisection(double (*fntocall)(double param1, double *param2),
                 double *secondreturn, double lowerguess, double upperguess);
double fnQy(double guess, double *Qy);
double fnEllipse(double guessd, double *c);
void MakeCornerSegBdaries(double *BdaryAngle, double *c, double *d);
double DistBetwPts(double ax, double ay, double bx, double by);
void SetUpGeometry(void);
double CalcCornerSpecialPropellantMix(int SegNo);
void RecordThrustProfile(int irecord);
void RecordThrustProfileFine(int irecord);
void SetupNozzleGeometry(void);
void SprayNozzleRays(void);
void MakeNozzleProfileRedefinitionLayer(void);
void FindNozzleRayPropellants(void);
double CalcNozzleRaysTotalIntegralContribs(void);
void CalcNozRayRDots(void);
void SetupInsulationGeometry(void);

void InitializeRemainingGlobals(void);
void InitializationsNeededForDecoding();
void CalcSegBurnAreas(void);
void CalcBurnIntegralContributions(void);
int CalcSegRefBurnRates(void);
int PressureViolation(void);
int ReadGenome(int *AlreadyJammedFlag);
int WriteAnimationFile(long neval, double Time);

int UpdateBurnRadiiAndReachedShell(void);
void WhenBurnHitsShell(int SegNo);

void WriteThisGuyOutputFile(void);
int DistributeHEEDSInputs(void);

double ThetaSlice; /* Now the same for all dome ring segments */
double ellipsec[NCORNERSEGS + 1], ellipsed[NCORNERSEGS + 1];
double B[NTOTALSEGS]; /* point along longer side of corner seg where special
                         fast */
                      /* propellant terminates. */
int CurrLayerIndex[NTOTALSEGS];
int ThrustInterval;
double RemainingDistToBurnInSeg[NTOTALSEGS];
static int FirstInitGlobalsCall = 1; /* For stuff to be done only once */
int AlreadyJammedFlag = 1;
double THIS_ANGLE; /* Kluge to communicate between MakeCornerSegBdaries() and
                      fnEllipse() */

/* Adding Phase 2 variables and objective, DeltaV: */
#ifdef PHASE2
void CalcThroatAreaDependentParams(void);
double mbo; /* mass at burnout */
double exitDiameter = 0.14605;
double amp = 113.3;
double bmp = 2755.;
double epsilon0 = 36.;
double epsilon;
double mprop;
#endif

/* Now vars used for nozzle segment */
int NozBurningRays = 0;
double NozRayAngle[NNOZRAYS];
double NozEndFastBurnPts[NNOZRAYS];
double NozRayFinalBurnRadius[NNOZRAYS];
double NozShortestRayAngle;
double NozRayBurnStartRadius;
double NozRayRadius[NNOZRAYS];
double NozRayBurnPty[NNOZRAYS];
double NozRayRefBurnRate[NNOZRAYS];
double NozRaySegBurnArea[NNOZRAYS];
double NozRayRDot[NNOZRAYS];

/* Adding vars used for insulation burn */
double InsulDensity = 1250.;        /* in kg/m^3 */
double InsulRefBurnRate = 0.000635; /* in m/sec */
double InsulBurnArea;
double InsulBottomBurnDepth;
double InsulPtx, InsulPty;
/* Next pt is movingintersection of nozzle propellant burn line with insulation
 * top */
double InsulPropBdaryPtx, InsulPropBdaryPty;
double InsulShellPtx, InsulShellPty;
double InsulTopPtAtShellx, InsulTopPtAtShelly;
double InsulTopBurningLength;
double InsulTopSlope;
double InsulTopBurnLineAngle;
double InsulTopBurnNormal;
double InsulHalfAngle;
double InsulPointBurnRadius;
double InsulPropBdaryBurnArea, InsulBottomBurnArea;
double InsulPropBdaryPerimeter;
double InsulDeltaBottomDepth;

/* Special debugging flag */
int RayDebug = 0;

int rayDepthFlag = 0;

void objectivefunction(
    int *ReasonStopped, double *lThrustReward, double *lSimultaneityReward,
    double *lTimes, double *lThrustProfile, double *lmPressureProfile,
    double *lTimesFine, double *lThrustProfileFine,
    double *lmPressureProfileFine, double *lBurnDistRemaining,
    int *lPropellantToEval, double *lLayerStartRadius, int *lRayBurnDepthInput,
    int *lPreCylStarRDotRefInput, int *lRDotRefCatchupInput,
    double lThroatAreaInput, double *lDeltaV, double lMaxThrust,
    int *lRayBurnDepths, int *lStarBurningThisSegNow,
    double *lTimeStarBurnFinish,
    double lSegResidualPerTimeStep[][NumTimePointsCoarse],
    double lSegResidualPerTimeStepFine[][NumTimePointsFine],
    double lTargetThrusts[NumTimePointsFine], int lRayDepthFlag,
    double lSegBurnAreaPerTimeStepCoarse[][NumTimePointsCoarse],
    double lSegBurnAreaPerTimeStepFine[][NumTimePointsFine],
    double lBurnLayerPerTimestepCoarse[][NumTimePointsCoarse],
    double lBurnLayerPerTimestepFine[][NumTimePointsFine]) {
  rayDepthFlag = lRayDepthFlag;
  for (int t = 0; t < NumTimePointsFine; t++) {
    TargetThrusts[t] = lTargetThrusts[t];
  }
  double CallingPressure;
  int SegNo, istep = 0, i, k;
  int Code;
  irecord = 0;
  irecordFine = 0;

  int Debug = 0;

  if (Debug)
    if (!(outfp = fopen("./OutputFile", "w"))) {
      printf("Output file %s cannot be written.\n", "OutputFile");
    }
    /* Unpack the calling parameters into global vars */
#ifdef PHASE2
  ThroatArea = lThroatAreaInput;
  DeltaV = *lDeltaV;
  MaxThrust = lMaxThrust;
#endif
  for (int SegNo = 0; SegNo < NTOTALSEGS; SegNo++) {
    LayerStartRadius[NLAYERS][SegNo] =
        10.; /* way outside rocket; shouldn't be crossed, so */
             /* never look past layer NLAYERS's propellant type */
    for (int i = 0; i < NLAYERS; i++) {
      PropellantToEval[i + NLAYERS * SegNo] =
          lPropellantToEval[i + NLAYERS * SegNo];
      LayerStartRadius[i][SegNo] = lLayerStartRadius[i + NLAYERS * SegNo];
    }
  }
#ifdef STARGRAIN
  for (int i = 0; i < ((NCYLSEGS - 3) * (NRAYS + 1)); i++) {
    RayBurnDepthInput[i] = *(lRayBurnDepthInput + i);
  }
  for (int i = 0; i < (NCYLSEGS - 3); i++) {
    PreCylStarRDotRefInput[i] = *(lPreCylStarRDotRefInput + i);
    RDotRefCatchupInput[i] = *(lRDotRefCatchupInput + i);
  }
#endif
  if (!DistributeHEEDSInputs()) {
    fprintf(outfp, "\n*** Inputs from NSGA-II or HEEDS could not be processed. "
                   " Quitting ***");
    exit(-222);
  }

  decodechromosome(); /* or here, decode the values passed in */

  InitializeRemainingGlobals();
  /* Now set up the geometry for the "star" burn special interior sections. */
  /* Must follow decoding of chromosome to work */
#ifdef STARGRAIN
  SetupStarBurn(SegNo);
  StarAreaCalc();
#endif /* of STARGRAIN */
  CalcPropellantMass();
#ifdef PHASE2
  CalcThroatAreaDependentParams();
#endif

  /* NOW READY TO DO BURN SIMULATION FOR THIS PARTICULAR GRAIN DESIGN */
  while (Burning == 1) {
    for (SegNo = 0; SegNo < NTOTALSEGS; SegNo++) {
      /* Initialize for this step to receive additions to calculate new value */
      NextSegRadius[SegNo] = SegmentRadius[SegNo];
    }
#ifdef STARGRAIN
    CalcStarRayBurnDirnsAndLengths();
#endif
    /* Calculate the reference burn rates for all segments */
    /* At same time, are checking PREVIOUS step to update CurrLayerIndex and
     * check for hit shell */
    Code = CalcSegRefBurnRates();
    if (Code == -1) {
      BurnedToShellFinish();
      goto FinishUp;
    }
    /* Calculate segment burn areas for all segments */
    CalcSegBurnAreas();
    CalcBurnIntegralContributions();

    /* NEWTON CALL HERE!!! Calls by VALUE with pressure, then resets pressure to
     * value returned. */
    CallingPressure = (TotalPressure + RefPressure) / 2.;
    TotalPressure = newton(CallingPressure, ThroatArea);
    if (BadPressure()) {
      Fitness = calc_fitness();
      goto FinishUp;
    }
    if (Debug)
      PrintAllSegs();
    /* Pressure okay, so reward burn time */
    Reward += 24000. * DeltaTime; /* Adds 600. per timestep of 0.025sec, less
                                     per step as gets shorter */
    /* Here is where changes in burn radius are done after current step's burn,
     * and */
    /* thrust is calculated from pressure in newton() */
    for (SegNo = 0; SegNo < NTOTALSEGS;
         SegNo++) { /* Equation (5) in Kibbey doc. */
      if (SegNo == NTOTALSEGS - 1 && NozBurningRays) {
        CalcNozRayRDots();
      }
#ifdef STARGRAIN
      else if (StarBurningThisSegNow[SegNo]) {
        CalcPerRayRDot(SegNo);
      }
#endif
      else {
// Modification by Abhiroop Ghosh. Sep 16, 2019
#ifdef STARGRAIN
        if (SegNo >= 2 && SegNo <= (NCYLSEGS - 2) &&
            TimeStarBurnFinish[SegNo] == 0) {
          TimeStarBurnFinish[SegNo] = Time;
          if (Debug) {
            fprintf(outfp, "\nTimeStarBurnFinish [%d] = %0.14f\n", SegNo, Time);
          }
        }
#endif
        RDot[SegNo] =
            RefBurnRate[SegNo] * pow(TotalPressure, Alpha) * InvPrefAlpha;
      }
    }
#ifdef PHASE2
    //    Thrust = Cf*At*Pc; Cf = Isp/C* = Isp/1555.;
    Thrust = Isp / 1555. * ThroatArea * (TotalPressure - AmbientPressure);
#else
    Thrust = 269.0 * G0 * ThroatArea * (TotalPressure - AmbientPressure) /
             CharacteristicVelocity;
#endif
    //      // Added by Abhiroop 09/24/2019. Records data on a finer time
    //      interval RecordThrustProfileFine(irecordFine);
    //
    //      // Record segment-wise residual per time step
    //      for (int s=0; s < NTOTALSEGS; s++)
    //      {
    //        SegResidualPerTimeStepFine[s][irecordFine] = FinalBurnRadius[s] -
    //        SegmentRadius[s];
    //      }
    //      irecordFine++;
    //
    //      if (istep % ThrustInterval == 0) {
    //         RecordThrustProfile(irecord);
    //
    //         // Added by Abhiroop 09/24/2019. Records data at 0.5 second
    //         interval
    //         // Record segment-wise residual per time step
    //          for (int s=0; s < NTOTALSEGS; s++)
    //          {
    //            SegResidualPerTimeStep[s][irecord] = FinalBurnRadius[s] -
    //            SegmentRadius[s];
    //          }
    //         irecord++;
    //      }
    RecordThrustAndResidualData(istep);
    if (UpdateBurnRadiiAndReachedShell()) {
      StopCode = 0;
      Time += DeltaTime;
      goto FinishUp;
    }
    /* At this point, NO LAYER HAS YET hit the shell, keep burning */
    /* Final actions in this timestep:  HERE MUST UPDATE BURN DEPTH ACC TO THIS
     * TIMESTEP'S BURN RATE */
    Time += DeltaTime;
    if (WantToAnimate == neval) {
      WriteAnimationFile(neval, Time);
    }
    if (Debug) {
      fprintf(outfp,
              "\nJust updated Time******** to %0.14f by DeltaTime %0.14f\n",
              Time, DeltaTime);
    }
    Burning = (Time + DeltaTime) < TotalBurnTime;
    for (SegNo = 0; SegNo < NTOTALSEGS; SegNo++) {
      SegmentRadius[SegNo] = NextSegRadius[SegNo];
    }
    if (NozBurningRays) {
      for (i = 0; i < NNOZRAYS; i++) {
        NozRayRadius[i] = NextNozRayRadius[i];
      }
    }
    istep++;

    if (Debug) {
      fflush(outfp);
    }
    /* Done with 1 timestep within Burning loop for this rocket */
  }
  /* If/when reach here, Burning loop is over, done with this rocket, burn
   * reached 10 seconds, reward. */
  Reward += 10000.;
  Fitness = calc_fitness();
  if (Debug) {
    fprintf(
        outfp,
        "\nExiting, 15 second burn completed for rocket neval %ld, gen %d.\n",
        neval, gen);
  }
  if (Debug) {
    fclose(outfp);
  }
  StopCode = 3;
FinishUp:
  // WriteThisGuyOutputFile();
#ifdef PHASE2
  *lDeltaV = DeltaV;
#endif
  RecordThrustAndResidualData(0);
  *ReasonStopped = StopCode;
  *lThrustReward = ThrustReward;
  *lSimultaneityReward = SimultaneityReward;
  // Added by Abhiroop 09/24/2019.
  // Record thrust profile in 0.5 second intervals
  for (i = 0; i < NumTimePointsCoarse; i++) {
    lTimes[i] = Times[i];
    lThrustProfile[i] = ThrustProfile[i];
    lmPressureProfile[i] = mPressureProfile[i];
  }
  // Added by Abhiroop 09/24/2019.
  // Record thrust profile in 0.025 second intervals
  for (i = 0; i < NumTimePointsFine; i++) {
    lTimesFine[i] = TimesFine[i];
    lThrustProfileFine[i] = ThrustFine[i];
    lmPressureProfileFine[i] = mPressureFine[i];
  }
  for (int i = 0; i < (NCYLSEGS - 3) * (NRAYS + 1); i++) {
    lRayBurnDepths[i] = RayDepthHEEDS[i];
  }
  // Added by Abhiroop Ghosh. Sep 16, 2019
  // Latest Update Mar 3, 2020 -- record which layer burning at each timestep
  for (int i = 0; i < NTOTALSEGS; i++) {
    lStarBurningThisSegNow[i] = StarBurningThisSegNow[i];
    if (StarBurningThisSegNow[i]) {
      lTimeStarBurnFinish[i] = -1;
    } else {
      lTimeStarBurnFinish[i] = TimeStarBurnFinish[i];
    }
    lBurnDistRemaining[i] = BurnDistRemaining[i];

    for (int k = 0; k < NumTimePointsCoarse; k++) {
      lSegResidualPerTimeStep[i][k] = SegResidualPerTimeStep[i][k];
      lSegBurnAreaPerTimeStepCoarse[i][k] = SegBurnAreaPerTimeStepCoarse[i][k];
      lBurnLayerPerTimestepCoarse[i][k] = BurnLayerPerTimestepCoarse[i][k];
    }

    for (int k = 0; k < NumTimePointsFine; k++) {
      lSegResidualPerTimeStepFine[i][k] = SegResidualPerTimeStepFine[i][k];
      lSegBurnAreaPerTimeStepFine[i][k] = SegBurnAreaPerTimeStepFine[i][k];
      lBurnLayerPerTimestepFine[i][k] = BurnLayerPerTimestepFine[i][k];
    }
  }

  return;
}

void RecordThrustAndResidualData(int istep) {
  // Added by Abhiroop 24th Sept. 2019. Records data on a finer time interval
  RecordThrustProfileFine(irecordFine);

  // Record segment-wise residual per time step
  for (int s = 0; s < NTOTALSEGS; s++) {
    SegResidualPerTimeStepFine[s][irecordFine] =
        FinalBurnRadius[s] - SegmentRadius[s];
    SegBurnAreaPerTimeStepFine[s][irecordFine] = SegmentArea[s];
    BurnLayerPerTimestepFine[s][irecordFine] = CurrLayerIndex[s];
  }
  irecordFine++;

  if (istep % ThrustInterval == 0) {
    RecordThrustProfile(irecord);

    // Added by Abhiroop 24th Sept. 2019. Records data at 0.5 second interval
    // Record segment-wise residual per time step
    for (int s = 0; s < NTOTALSEGS; s++) {
      SegResidualPerTimeStep[s][irecord] =
          FinalBurnRadius[s] - SegmentRadius[s];
      SegBurnAreaPerTimeStepCoarse[s][irecord] = SegmentArea[s];
      BurnLayerPerTimestepCoarse[s][irecord] = CurrLayerIndex[s];
    }
    irecord++;
  }
}

void PrintAllSegs(void) {
  int SegNo, k;

  /* DEBUG PRINTOUTS OF ALL SEGS */
  for (SegNo = 0; SegNo < NTOTALSEGS; SegNo++) {
    /* Note that CYLSEGS will not still be burning rays when any */
    /* seg burns out, due to fastest/slowest ref burn rate ratio */
    if (SegNo != NTOTALSEGS - 1) {
#ifdef STARGRAIN
      if (StarBurningThisSegNow[SegNo]) {
        fprintf(outfp, "\nStar Burning seg [%d]", SegNo);
        for (k = 0; k < NRAYS + 2; k++) {
          fprintf(outfp,
                  "\nFinalBurnRadius[%d] %0.14f, RayBurnDepth[%d][%d] %0.14f",
                  SegNo, FinalBurnRadius[k], k, SegNo, RayBurnDepth[k][SegNo]);
        }
      } else {
        fprintf(outfp,
                "\nIn objfunc(), FinalBurnRadius[%d] %0.14f, SegmentRadius[%d] "
                "%0.14f",
                SegNo, FinalBurnRadius[SegNo], SegNo, SegmentRadius[SegNo]);
        fprintf(outfp, "\nSeg %d contributes %0.14f", SegNo,
                fabs(FinalBurnRadius[SegNo] - SegmentRadius[SegNo]) * 100.);
      }
#else
      fprintf(outfp,
              "\nIn objfunc(), FinalBurnRadius[%d] %0.14f, SegmentRadius[%d] "
              "%0.14f",
              SegNo, FinalBurnRadius[SegNo], SegNo, SegmentRadius[SegNo]);
      fprintf(outfp, "\nSeg %d contributes %0.14f", SegNo,
              fabs(FinalBurnRadius[SegNo] - SegmentRadius[SegNo]) * 100.);
#endif
    }
    /* Nozzle seg could still be burning cylindrically, or as rays */
    else { /* Nozzle seg */
           /*Prints of nozzle seg or ray burns */
      if (NozBurningRays) {
        for (k = 0; k < NNOZRAYS; k++) {
          fprintf(outfp,
                  "\nIn objfunc(), NozBurningRays,  NozRayFinalBurnRadius[%d] "
                  "%0.14f, NozRayRadius[%d] %0.14f",
                  k, NozRayFinalBurnRadius[k], k, NozRayRadius[k]);
        }
      } else {
        fprintf(outfp,
                "\nIn objfunc(), FinalBurnRadius[%d] %0.14f, SegmentRadius[%d] "
                "%0.14f",
                SegNo, FinalBurnRadius[SegNo], SegNo, SegmentRadius[SegNo]);
        fprintf(outfp, "\nSeg %d contributes %0.14f", SegNo,
                fabs(FinalBurnRadius[SegNo] - SegmentRadius[SegNo]) * 100.);
      }
    }
  }
}

int BadPressure(void) {
  int Debug = 0;

  if (TotalPressure == -444.44444) {
    /* Jacobian generated a negative pressure, failed; abandon this rocket. */
    Fitness = calc_fitness();
    return 1;
  }
  if (TotalPressure == -9999.) {
    /* Jacobian generated NaN, failed; abandon this rocket. */
    Fitness = calc_fitness();
    Fitness = 1.; /* Overriding with a very low fitness. */
    return 1;
  }
  if (PressureViolation()) {
    if (TotalPressure == -555.55555) {
      /* PressureViolation() found bad pressure, calculate fitness, end the run
       */
      Fitness = calc_fitness();
      if (Debug)
        fprintf(outfp,
                "\n *** Pressure violation ended eval of rocket %ld, Gen %d, "
                "at time %0.14f",
                neval, gen, Time);
    }
    return 1;
  }
  return 0;
}

void CalcPropellantMass() {

  int Debug = 0;

  int SegNo;
  double cornerdomecylpartheight;
  double DomePropellVolume, DomeCylVolume, DomePropellMass;
  double nozcylvol, noztriangvol, begxoftriangle, endxoftriangle, triangtop,
      triangbot;
  double NozPropellVolume, NozPropellMass;
  double nozbiginsultrianglevol, bigtriangleheight;
  double smallinsultriangvolcorrn, upperradsmalltriangle;
  double topmissingtriangcorrn;

#ifdef STARGRAIN
  int RayNo;
#endif

  PropellantTotalVolume = 0.;
  PropellantTotalMass = 0.;

  /* See Excel spreadsheet for details of derivation of this model */
  /* First handle non-star cylinder segments: diff of cyl volumes. */
  for (SegNo = 1; SegNo < NCYLSEGS; SegNo++) {
    PropellVolume[SegNo] = PI *
                           (FinalBurnRadius[SegNo] * FinalBurnRadius[SegNo] -
                            SegmentRadius[SegNo] * SegmentRadius[SegNo]) *
                           SegmentLength[SegNo];
    /* Now add in the "star-shaped" portions of those segments currently burning
     * stars. */
    /* Treats arc where star joins cyl seg is just a line. */
#ifdef STARGRAIN
    if (StarBurningThisSegNow[SegNo]) {
      for (RayNo = 1; RayNo < NRAYS + 2; RayNo++) {
        PropellVolume[SegNo] +=
            NSYMS * StarRayArea[RayNo][SegNo] * SegmentLength[SegNo];
      }
    }
    PropellMass[SegNo] = PropellVolume[SegNo] * PropellantDensity;
    if (Debug) {
      fprintf(outfp, "\nCyl seg %d mass %0.14f", SegNo, PropellMass[SegNo]);
    }
#endif
    PropellantTotalVolume += PropellVolume[SegNo];
    PropellantTotalMass += PropellVolume[SegNo] * PropellantDensity;
  }
  if (Debug) {
    fprintf(outfp, "\nTotal of non-dome cyl segs mass, incl stars %0.14f",
            PropellantTotalMass);
  }

  /* Now dome (segment 0 and corner segments, together */
  DomePropellVolume = (4. / 3.) * PI * SemiMajAx * SemiMajAx * SemiMinAx;
  /* cyl part height = WT - semiminor axis + epsilon for curve avgd across width
   */
  cornerdomecylpartheight =
      (FinalBurnRadius[1] - InitBurnRadius[1]) - SemiMinAx + 0.00337;
  DomeCylVolume =
      PI * FinalBurnRadius[1] * FinalBurnRadius[1] * cornerdomecylpartheight;
  DomePropellMass = (DomePropellVolume + DomeCylVolume) * PropellantDensity;
  PropellantTotalVolume += DomePropellVolume + DomeCylVolume;
  PropellantTotalMass += DomePropellMass;
  if (Debug) {
    fprintf(outfp, "\nDome propellant mass %0.14f", DomePropellMass);
  }

  /* This completes all segments except the nozzle. */
  /* Nozzle segment is partly cylindrical minus a correction for small
   * insulation triangle, */
  begxoftriangle = 0.4 * (FinalBurnRadius[1] - InitBurnRadius[1]);
  upperradsmalltriangle =
      InitBurnRadius[1] + begxoftriangle * tan((17. / 180.) * PI);
  smallinsultriangvolcorrn = 0.5 * 0.4 * PI *
                             (FinalBurnRadius[1] - InitBurnRadius[1]) *
                             (upperradsmalltriangle * upperradsmalltriangle -
                              InitBurnRadius[1] * InitBurnRadius[1]);
  nozcylvol = PI * 0.8 * (FinalBurnRadius[1] - InitBurnRadius[1]) *
                  (FinalBurnRadius[1] * FinalBurnRadius[1] -
                   InitBurnRadius[1] * InitBurnRadius[1]) -
              smallinsultriangvolcorrn;
  /* Nozzle seg also includes a cut-off triangular section also reduced by a
   * truncated triangle at top */
  triangtop = 1.15 * FinalBurnRadius[1];
  triangbot = InitBurnRadius[1];
  endxoftriangle = 1.2 * (FinalBurnRadius[1] - InitBurnRadius[1]);
  noztriangvol = 0.5 * PI * (endxoftriangle - begxoftriangle) *
                 (triangtop * triangtop - triangbot * triangbot);
  topmissingtriangcorrn =
      0.5 * 0.15 * (FinalBurnRadius[1] - InitBurnRadius[1]) * PI *
      (triangtop * triangtop - FinalBurnRadius[1] * FinalBurnRadius[1]);
  noztriangvol -= topmissingtriangcorrn;
  if (Debug) {
    fprintf(outfp,
            "\ntriangtop %e, triangbot %e, endxoftriangle %e, noztriangvol %e, "
            "topmissingtriangcorrn %e, noztriangvol %e, PI %e",
            triangtop, triangbot, endxoftriangle, noztriangvol,
            topmissingtriangcorrn, noztriangvol, PI);
  }
  /* now correct noztriangvol for insulation loss by */
  /* subtracting (big insul triangle - small insul triangle) */
  bigtriangleheight =
      1.2 * (FinalBurnRadius[1] - InitBurnRadius[1]) * tan((17. / 180.) * PI);

  if (Debug) {
    fprintf(
        outfp,
        "\nFinalBurnRadius[1] %e, InitBurnRadius[1] %e, bigtriangleheight %e",
        FinalBurnRadius[1], InitBurnRadius[1], bigtriangleheight);
  }
  nozbiginsultrianglevol = 1.2 * (FinalBurnRadius[1] - InitBurnRadius[1]) *
                           0.5 * PI *
                           (bigtriangleheight * bigtriangleheight -
                            InitBurnRadius[1] * InitBurnRadius[1]);
  noztriangvol -= (nozbiginsultrianglevol - smallinsultriangvolcorrn);
  NozPropellVolume = nozcylvol + noztriangvol;
  NozPropellMass = NozPropellVolume * PropellantDensity;
  if (Debug) {
    fprintf(outfp,
            "\nbigtriangleheight %e, nozbiginsultrianglevol %e, noztriangvol "
            "%e, NozPropellVolume %e",
            bigtriangleheight, nozbiginsultrianglevol, noztriangvol,
            NozPropellVolume);
    fprintf(outfp, "\nNozzle propellant mass %e", NozPropellMass);
  }

  PropellantTotalVolume += NozPropellVolume;
  PropellantTotalMass += NozPropellMass;
  if (Debug) {
    fprintf(outfp, "\nTotal propellant volume %e, total propellant mass %e",
            PropellantTotalVolume, PropellantTotalMass);
  }
  return;
}

#ifdef PHASE2
void CalcThroatAreaDependentParams(void) {
  int Debug = 0;

  ThroatDiameter = 2. * sqrt(ThroatArea / PI);
  if (Debug)
    fprintf(outfp, "\nThroat Diameter %0.14f", ThroatDiameter);
  epsilon = (exitDiameter / ThroatDiameter) * (exitDiameter / ThroatDiameter);
  Isp = amp * log(epsilon / epsilon0) + bmp;
  /*    Thrust = Cf*At*Pc; Cf = Isp/C* = Isp/1555.; */
  /*    MaxPressure = MaxThrust/(Cf*At) + AmbientPressure */
  /*    MaxPressure = MaxThrust/((Isp/C*)*ThroatArea) + AmbientPressure */
  /*    MaxThrust was determined by the desired target thrust profile, if doing
   * PHASE2 */
  MaxPressure = (MaxThrust / ((Isp / 1555.) * ThroatArea) + AmbientPressure);
  mbo = 12.13 + 0.164 * (MaxPressure / 1000000. - 3.42);
  mprop = PropellantTotalMass;
  DeltaV = Isp * log(1. + mprop / mbo);
  if (Debug) {
    fprintf(outfp,
            "\nneval %ld, Isp %0.14f, MaxPressure %0.14f, mbo %0.14f, mprop "
            "%0.14f, DeltaV %0.14f",
            neval, Isp, MaxPressure, mbo, mprop, DeltaV);
  }
}
#endif

void StarAreaCalc(void) {
  int SegNo, RayNo;

  int Debug = 0;

  double a, b, c, d, p, q;
  /* From current data structure, calc area between adjacent rays */
  /* Uses area formula for gen quadrilateral area=1/4*sqrt(4p^2*q^2 -
   * (b^2+d*2-a^2-c^2)^2) */
  for (SegNo = 2; SegNo < NCYLSEGS - 1; SegNo++) {
    for (RayNo = 1; RayNo < NRAYS + 2; RayNo++) {
      /* No area for ray 0 */
      a = DistBetwPts(BurnPtx[RayNo][SegNo], BurnPty[RayNo][SegNo],
                      BurnPtx[RayNo - 1][SegNo], BurnPty[RayNo - 1][SegNo]);
      b = DistBetwPts(BurnPtx[RayNo - 1][SegNo], BurnPty[RayNo - 1][SegNo],
                      RayCylIntersectx[RayNo - 1][SegNo],
                      RayCylIntersecty[RayNo - 1][SegNo]);
      c = DistBetwPts(RayCylIntersectx[RayNo - 1][SegNo],
                      RayCylIntersecty[RayNo - 1][SegNo],
                      RayCylIntersectx[RayNo][SegNo],
                      RayCylIntersecty[RayNo][SegNo]);
      d = DistBetwPts(BurnPtx[RayNo][SegNo], BurnPty[RayNo][SegNo],
                      RayCylIntersectx[RayNo][SegNo],
                      RayCylIntersecty[RayNo][SegNo]);
      p = DistBetwPts(BurnPtx[RayNo - 1][SegNo], BurnPty[RayNo - 1][SegNo],
                      RayCylIntersectx[RayNo][SegNo],
                      RayCylIntersecty[RayNo][SegNo]);
      q = DistBetwPts(BurnPtx[RayNo][SegNo], BurnPty[RayNo][SegNo],
                      RayCylIntersectx[RayNo - 1][SegNo],
                      RayCylIntersecty[RayNo - 1][SegNo]);
      StarRayArea[RayNo][SegNo] =
          0.25 *
          pow(4. * p * p * q * q - pow((b * b + d * d - a * a - c * c), 2),
              0.5);
      if (Debug)
        fprintf(outfp,
                "\n a %0.14f, b %0.14f, c %0.14f, d %0.14f, p %0.14f, q "
                "%0.14f, StarRayArea[%d][%d] %0.14f",
                a, b, c, d, p, q, RayNo, SegNo, StarRayArea[RayNo][SegNo]);
    }
  }
  return;
}

void BurnedToShellFinish() {

  int Debug = 0;

  /* A ray burned through to shell; fitness not yet calculated, call it here */
  Fitness = calc_fitness();
  if (Debug) {
    fprintf(outfp,
            "\n *** A ray burned through to shell, stopping eval of rocket "
            "%ld, Gen %d, at time %0.14f\nThrustReward %0.14f, "
            "SimultaneityReward %0.14f, StdDevReward %0.14f, fitness %0.14f.\n",
            neval, gen, Time, ThrustReward, SimultaneityReward, StdDevReward,
            Reward);
  }
}

void InitializationsNeededForDecoding(void) {

  int SegNo;

  int Debug = 0;

  /* Things that must be set before chromosome is decoded */
  for (SegNo = 0; SegNo < NTOTALSEGS; SegNo++) {
    /* SET UP TO BURN THE CYLINDRICAL SEGMENTS HERE */
#ifndef EVOLVEINITBURNRADIUS
    /* If evolving it, has already been set. */
    InitBurnRadius[SegNo] = INITBURNRADIUS; /* Unless overridden if evolving it
                                               (in decodechromosome()). */
#endif
    FinalBurnRadius[SegNo] = RocketCylRadius;
#ifdef STARGRAIN
    if (SegNo > 1 && SegNo < NCYLSEGS - 1) {
      PtAx[SegNo] = INNERSTARPOINT;
      if (Debug)
        fprintf(outfp, "\nPtAx[%d] %0.14f", SegNo, PtAx[SegNo]);
    }
    /* Unless overridden if evolving it in decodechromosome() */
#endif
  }
  if (Debug)
    fprintf(outfp, "\n Nozzle FinalBurnRadius %0.14f",
            FinalBurnRadius[NTOTALSEGS - 1]);

  /* Now handle segment 0, the central elliptical dome segment. */
  /* Segment[0] length setting, the ellipse */
  /* SegmentLength[0] will change with depth of burn, must calculate below. */
  // InitBurnRadius[NCYLSEGS-1] is the radius for which I want Px and Qy
  // calculated.
  /* Calculate the Qy point where the burn area of the elliptical dome will
   * start, and Px, */
  /* the point where the outer edge of ellipse central portion burn stops, given
   */
  /* the value InitBurnRadius[NCYLSEGS], the initial burn radius of the segment
   * adjacent to the */
  /* elliptical dome segment. */
  MajAx = RocketCylRadius * 2.;
  MinAx = MajAx / Eccentricity;
  SemiMajAx = MajAx / 2.;
  SemiMinAx = MinAx / 2.;
  Qx = InitBurnRadius[NCYLSEGS - 1];
  InitBurnRadius[0] = InitBurnRadius[NCYLSEGS - 1];
  FinalBurnRadius[0] = RocketCylRadius;

  // original : Px = bisection(&fnQy, &Qy, 3., 4.); /* Px returned as a new
  // value, Qy set inside function */
  Px = bisection(&fnQy, &Qy, InitBurnRadius[0] * 100.,
                 5.0); /* Px returned as a new value, Qy set inside function */
  /* Convert the Px, Py, Qx, Qy to meters. */
  Px /= 100.;
  Qx /= 100.;
  Py /= 100.;
  Qy /= 100.;
  if (Debug)
    fprintf(outfp, "\nPx %0.14f, Py %0.14f, Qx %0.14f, Qy %0.14f\n", Px, Py, Qx,
            Qy);
  /* Now initialize the segment boundary points and angles for the corner "ring"
   * cap dome segments */
  MakeCornerSegBdaries(BdaryAngle, ellipsec, ellipsed);
  /* Set the outermost corner boundary lengths already known. */
  BdaryLen[0] = InitBurnRadius[NCYLSEGS - 1];
  BdaryLen[NCORNERSEGS] = DistBetwPts(Px, Py, Qx, Qy);
  /* NOTE: BdaryLen is length of a boundary; (ellipsec[],ellipsed[]) is the
   * point on outer ellipse at end of the boundary. */
  for (SegNo = 0; SegNo < NCORNERSEGS; SegNo++) {
    if (Debug)
      fprintf(outfp,
              "\nFor CornerSegs, BdaryAngle[%d] %0.14f, ellipsec[%d] %0.14f, "
              "ellipsed[%d] %0.14f\n",
              SegNo, BdaryAngle[SegNo], SegNo, ellipsec[SegNo], SegNo,
              ellipsed[SegNo]);
    /* BdaryAngle now has angles from "corner" of web to each point
     * (ellipsec,ellipsed) on outer surface of ellipse. */
    /* Will use these in burn calculations below. */
    /* Calculate lengths of seg boundaries: */
    BdaryLen[SegNo] = DistBetwPts(ellipsec[SegNo], ellipsed[SegNo],
                                  InitBurnRadius[NCYLSEGS], Qy);
    if (Debug)
      fprintf(outfp, "\nBdaryLen[%d] %0.14f", SegNo, BdaryLen[SegNo]);
  }
  if (Debug)
    fprintf(outfp, "\nBdaryLen[%d] %0.14f", NCORNERSEGS, BdaryLen[NCORNERSEGS]);
  /* Now set seg length to mean of lengths of its two sides, for now. */
  for (SegNo = NCYLSEGS; SegNo < NTOTALSEGS - 1; SegNo++) {
    InitBurnRadius[SegNo] = InitBurnRadius[NCYLSEGS - 1];
    FinalBurnRadius[SegNo] =
        (BdaryLen[SegNo - NCYLSEGS] + BdaryLen[SegNo - NCYLSEGS + 1]) / 2. +
        InitBurnRadius[SegNo];
    if (Debug)
      fprintf(outfp, "\nCorner Seg FinalBurnRadius[%d] %0.14f\n", SegNo,
              FinalBurnRadius[SegNo]);
  }
  /* Set up burn radii for nozzle segment. */
  InitBurnRadius[NTOTALSEGS - 1] =
      InitBurnRadius[1]; /* matches adjacent cyl seg */
  /* Final burn radii are each different, along rays, so set up separately */
}

int WriteGraphicsOutputFile(long neval) {
  int i, SegNo;

  if (!(Graphfp = fopen("./GraphicsFile", "w"))) {
    printf("Graphics file %s cannot be opened for write.\n", "GraphicsFile");
    //   exit(777);
  }
  fprintf(Graphfp, "Evaluation number in THIS population, neval\n %ld", neval);
  fprintf(Graphfp, "\nRocket Radius\n %0.14f", 0.0762);
  fprintf(Graphfp, "\nCylinder length\n %0.14f", 0.5388);
  fprintf(Graphfp, "\nDome ellipse major axis\n %0.14f", 0.1524);
  fprintf(Graphfp, "\nDome ellipse minor axis\n %0.14f", 0.508);
  fprintf(Graphfp, "\nDome lower ellipse major axis\n %0.14f", 0.1524);
  fprintf(Graphfp, "\nDome lower ellipse minor axis\n %0.14f", 0.065);
  fprintf(Graphfp, "\nInitBurnRadius for all segs\n %0.14f", InitBurnRadius[0]);
  fprintf(Graphfp, "\nOffset between ellipses\n %0.14f", InitBurnRadius[0]);
  fprintf(Graphfp, "\nBoundary of dome seg on upper ellipse, x=\n %0.14f", Px);
  fprintf(Graphfp, "\nBoundary of dome seg on upper ellipse, y=\n %0.14f", Py);
  fprintf(Graphfp, "\nDomeCornerPointx on lower ellipse\n %0.14f", Qx);
  fprintf(Graphfp, "\nDomeCornerPointy on lower ellipse\n %0.14f", Qy);
  fprintf(Graphfp, "\nNZONES\n %d", NZONES);
  fprintf(Graphfp, "\nNLAYERSPERZONE\n %d", NLAYERSPERZONE);
  fprintf(Graphfp, "\nNLAYERS\n %d", NLAYERS);
  fprintf(Graphfp, "\nNCYLSEGS\n %d", NCYLSEGS);
  fprintf(Graphfp, "\nNCORNERSEGS\n %d", NCORNERSEGS);
  fprintf(Graphfp, "\nNTOTALSEGS\n %d", NTOTALSEGS);
  fprintf(Graphfp, "\nNNOZRAYS\n %d", NNOZRAYS);
  fprintf(Graphfp, "\nSTARGRAIN\n %s", "STARGRAIN");
#ifdef NSYMS
  fprintf(Graphfp, "\nNSYMS\n %d", NSYMS);
#endif
  fprintf(Graphfp, "\nNPROPELLTYPES\n %d", NPROPELLTYPES);

  /* SEGMENT BOUNDARY DEFINITIONS, ROCKET DIMENSIONS, etc. */
  /* Now dims for corner segs */
  fprintf(Graphfp, "\nbigtheta\n %0.14f", bigtheta);
  fprintf(Graphfp, "\nThetaSlice\n %0.14f", ThetaSlice);
  for (i = 0; i < NCORNERSEGS + 1; i++) {
    fprintf(Graphfp, "\nBdaryAngle[%d]\n %0.14f", i, BdaryAngle[i]);
    fprintf(Graphfp, "\nCorner ray bdary intersection x coord: C[%d]\n %0.14f",
            i, ellipsec[i]);
    fprintf(Graphfp, "\nCorner ray bdary intersection y coord: D[%d]\n %0.14f",
            i, ellipsed[i]);
  }
  /* Then dims for cyl segs, excluding the different dome cyl seg already
   * specified  */
  for (i = 1; i < NCYLSEGS; i++) {
    fprintf(Graphfp, "\nCylindrical Segment Length[%d]\n %0.14f", i,
            SegmentLength[i]);
  }
#ifdef STARGRAIN
  fprintf(Graphfp,
          "\nNumber of INTERNAL NRAYS (not counting seg boundary rays):\n %d",
          NRAYS);
  fprintf(Graphfp, "\nNRAYDEPTHS\n %d", NRAYDEPTHS);
#ifdef EVOLVINGSTARSHAPE
  fprintf(Graphfp, "\nEVOLVINGSTARSHAPE\n %d", EVOLVINGSTARSHAPE);
#endif
  for (SegNo = 2; SegNo < NCYLSEGS - 1; SegNo++) {
    /* if NOT EVOLVINGSTARSHAPE, then they determine the line which is the star
     * edge. */
    /* But IF EVOLVINGSTARSHAPE, then they are simply the intersections of the
     * star edge with rays */
    fprintf(Graphfp, "\nRayBurnDepth[%d][%d] = PtAx[%d]\n %0.14f", 0, SegNo,
            SegNo, PtAx[SegNo]);
    fprintf(Graphfp, "\nRayBurnDepth[%d][%d]  PtAy[%d]\n %0.14f", 0, SegNo,
            SegNo, PtAy[SegNo]);
    fprintf(Graphfp, "\nRay PtBRadius[%d]\n %0.14f", SegNo, PtBRadius[SegNo]);
    fprintf(Graphfp, "\nRay PtBx[%d]\n %0.14f", SegNo, PtBx[SegNo]);
    fprintf(Graphfp, "\nRay PtBy[%d]\n %0.14f", SegNo, PtBy[SegNo]);
    for (i = 1; i < NRAYS + 1; i++) {
      fprintf(Graphfp, "\nSeg %d Ray %d Init Burn Depth %0.14f", SegNo, i,
              RayBurnDepth[i][SegNo]);
      fprintf(Graphfp, "\nSeg %d Ray %d Init Burn Pt x %0.14f", SegNo, i,
              BurnPtx[i][SegNo]);
      fprintf(Graphfp, "\nSeg %d Ray %d Init Burn Pt y %0.14f", SegNo, i,
              BurnPty[i][SegNo]);
    }
    for (i = 0; i < NRAYS + 2; i++) {
      fprintf(
          Graphfp,
          "\nSpecial faster propellant catchup bdary depth[%d][%d]\n %0.14f", i,
          SegNo, SpecLayerOffset[i][SegNo]);
    }
  }
#endif /* of STARGRAIN */
  /* Then dims for nozzle */
  for (i = 0; i < NNOZRAYS; i++) {
    fprintf(Graphfp, "\nNozRayAngle[%d]\n %0.14f", i, NozRayAngle[i]);
  }
  for (i = 0; i < NNOZRAYS; i++) {
    fprintf(Graphfp, "\nNozEndFastBurnPts[%d]\n %0.14f", i,
            NozEndFastBurnPts[i]);
  }
  for (i = 0; i < NNOZRAYS; i++) {
    fprintf(Graphfp, "\nNozRayFinalBurnRadius[%d]\n %0.14f", i,
            NozRayFinalBurnRadius[i]);
  }
  /* Then dims for insulation */
  fprintf(Graphfp, "\nTip of insulation, x = \n%0.14f", InsulPtx);
  fprintf(Graphfp, "\nTip of insulation, y = \n%0.14f", InsulPty);
  fprintf(Graphfp, "\nShell end horizontal insulation line, x = \n%0.14f",
          InsulShellPtx);
  fprintf(Graphfp, "\nShell end horizontal insulation line, y = \n%0.14f",
          InsulShellPty);
  fprintf(Graphfp, "\nTop of insulation hits shell, x = \n%0.14f",
          InsulTopPtAtShellx);
  fprintf(Graphfp, "\nTop of insulation hits shell, y = \n%0.14f",
          InsulTopPtAtShelly);
  fprintf(Graphfp, "\nFurthest Burned Pt along top of insul, x = \n%0.14f",
          InsulPropBdaryPtx);
  fprintf(Graphfp, "\nFurthest Burned Pt along top of insul, y = \n%0.14f",
          InsulPropBdaryPty);
  fprintf(Graphfp, "\nInsulTopBurningLength\n %0.14f", InsulTopBurningLength);
  fprintf(
      Graphfp,
      "\nInsulHalfAngle along which top insulation burn is tracked:\n%0.14f",
      InsulHalfAngle);

  /* Now print layer start radii and propellanttypes for ALL segments, ALL
   * layers. */
  for (SegNo = 0; SegNo < NTOTALSEGS; SegNo++) {
    for (i = 0; i < NLAYERS; i++) {
      fprintf(Graphfp, "\nLayerStartRadius[%d][%d]\n %0.14f", i, SegNo,
              LayerStartRadius[i][SegNo]);
    }
    for (i = 0; i < NLAYERS; i++) {
      fprintf(Graphfp, "\nPropellType[%d][%d] %d", i, SegNo,
              PropellType[i][SegNo]);
    }
  }
  fprintf(Graphfp, "\nThat concludes the information defining the initial "
                   "configuration of the rocket!!!\n");

  fclose(Graphfp);
  Graphfp = NULL;
  return 1;
}

int WriteAnimationFile(long neval, double Time) {
  int k;
  int i, SegNo;
  double dummy;

  if (FirstAnimationCall) {
    FirstAnimationCall = 0;
    if (!(Animatfp = fopen("./AnimationFile", "w"))) {
      printf("Animation file %s cannot be opened for write.\n",
             "AnimationFile");
      //    exit(777);
    }
    fprintf(Animatfp, "Evaluation number in THIS population, neval\n %ld",
            neval);
    fprintf(Animatfp, "\nRocket Radius\n %0.14f", 0.0762);
    fprintf(Animatfp, "\nCylinder length\n %0.14f", 0.5388);
    fprintf(Animatfp, "\nDome ellipse major axis\n %0.14f", 0.1524);
    fprintf(Animatfp, "\nDome ellipse minor axis\n %0.14f", 0.508);
    fprintf(Animatfp, "\nDome lower ellipse major axis\n %0.14f", 0.1524);
    fprintf(Animatfp, "\nDome lower ellipse minor axis\n %0.14f", 0.065);
    fprintf(Animatfp, "\nInitBurnRadius for all segs\n %0.14f",
            InitBurnRadius[0]);
    fprintf(Animatfp, "\nOffset between ellipses\n %0.14f", InitBurnRadius[0]);
    fprintf(Animatfp, "\nBoundary of dome seg on upper ellipse, x=\n %0.14f",
            Px);
    fprintf(Animatfp, "\nBoundary of dome seg on upper ellipse, y=\n %0.14f",
            Py);
    fprintf(Animatfp, "\nDomeCornerPointx on lower ellipse\n %0.14f", Qx);
    fprintf(Animatfp, "\nDomeCornerPointy on lower ellipse\n %0.14f", Qy);
    fprintf(Animatfp, "\nNZONES\n %d", NZONES);
    fprintf(Animatfp, "\nNLAYERSPERZONE\n %d", NLAYERSPERZONE);
    fprintf(Animatfp, "\nNLAYERS\n %d", NLAYERS);
    fprintf(Animatfp, "\nNCYLSEGS\n %d", NCYLSEGS);
    fprintf(Animatfp, "\nNCORNERSEGS\n %d", NCORNERSEGS);
    fprintf(Animatfp, "\nNTOTALSEGS\n %d", NTOTALSEGS);
    fprintf(Animatfp, "\nNNOZRAYS\n %d", NNOZRAYS);
    fprintf(Animatfp, "\nSTARGRAIN\n %s", "STARGRAIN");
#ifdef NSYMS
    fprintf(Animatfp, "\nNSYMS\n %d", NSYMS);
#endif
    fprintf(Animatfp, "\nNPROPELLTYPES\n %d", NPROPELLTYPES);

    /* SEGMENT BOUNDARY DEFINITIONS, ROCKET DIMENSIONS, etc. */
    /* Now dims for corner segs */
    fprintf(Animatfp, "\nbigtheta\n %0.14f", bigtheta);
    fprintf(Animatfp, "\nThetaSlice\n %0.14f", ThetaSlice);
    for (SegNo = NCYLSEGS; SegNo < NCYLSEGS + NCORNERSEGS; SegNo++) {
      /* Need to force calculation of B[SegNo] here */
      dummy = CalcCornerSpecialPropellantMix(SegNo);

      fprintf(Animatfp, "\nBdaryAngle[%d]\n %0.14f", SegNo - NCYLSEGS,
              BdaryAngle[SegNo - NCYLSEGS]);
      fprintf(Animatfp,
              "\nCorner ray bdary intersection x coord: C[%d]\n %0.14f",
              SegNo - NCYLSEGS, ellipsec[SegNo - NCYLSEGS]);
      fprintf(Animatfp,
              "\nCorner ray bdary intersection y coord: D[%d]\n %0.14f",
              SegNo - NCYLSEGS, ellipsed[SegNo - NCYLSEGS]);
      if (SegNo < NCYLSEGS + NCORNERSEGS) {
        fprintf(Animatfp,
                "\nPoint B, depth (measured from start of layer %d) of faster "
                "propellant on longer side of corner segment[%d], B = %0.14f",
                NLAYERS - 2, SegNo, B[SegNo]);
      }
    }
    fprintf(Animatfp, "\nBdaryAngle[%d]\n %0.14f", SegNo - NCYLSEGS,
            BdaryAngle[SegNo - NCYLSEGS]);
    fprintf(Animatfp, "\nCorner ray bdary intersection x coord: C[%d]\n %0.14f",
            NCORNERSEGS, ellipsec[NCORNERSEGS]);
    fprintf(Animatfp, "\nCorner ray bdary intersection y coord: D[%d]\n %0.14f",
            NCORNERSEGS, ellipsed[NCORNERSEGS]);

    fprintf(Animatfp, "\nCorner segs slow special layer propellant type: %d",
            0);
    fprintf(Animatfp, "\nCorner segs fast special layer propellant type: %d",
            5);

    /* Then dims for cyl segs, excluding the different dome cyl seg already
     * specified  */
    for (SegNo = 1; SegNo < NCYLSEGS; SegNo++) {
      fprintf(Animatfp, "\nCylindrical Segment Length[%d]\n %0.14f", SegNo,
              SegmentLength[SegNo]);
    }
#ifdef STARGRAIN
    fprintf(Animatfp,
            "\nNumber of INTERNAL NRAYS (not counting seg boundary rays):\n %d",
            NRAYS);
    fprintf(Animatfp, "\nNRAYDEPTHS\n %d", NRAYDEPTHS);
    for (SegNo = 2; SegNo < NCYLSEGS - 1; SegNo++) {
      fprintf(Animatfp,
              "\nLayer number where compensation burning in Seg %d begins, %d",
              SegNo, StartCatchupLayer[SegNo]);
      fprintf(Animatfp,
              "\nStarting depth of layer where compensation burning begins in "
              "Seg %d, %0.14f",
              SegNo, LayerStartRadius[StartCatchupLayer[SegNo]][SegNo]);
    }
#ifdef EVOLVINGSTARSHAPE
    fprintf(Animatfp, "\nEVOLVINGSTARSHAPE\n %d", EVOLVINGSTARSHAPE);
#endif
    for (SegNo = 2; SegNo < NCYLSEGS - 1; SegNo++) {
      /* if NOT EVOLVINGSTARSHAPE, then they determine the line which is the
       * star edge. */
      /* But IF EVOLVINGSTARSHAPE, then they are simply the intersections of the
       * star edge with rays */
      fprintf(Animatfp, "\nRay PtAx[%d][%d]\n %0.14f", 0, SegNo, PtAx[SegNo]);
      fprintf(Animatfp, "\nRay PtAy[%d][%d]\n %0.14f", 0, SegNo, PtAy[SegNo]);
      fprintf(Animatfp, "\nRay PtBRadius[%d][%d]\n %0.14f", NRAYS + 1, SegNo,
              PtBRadius[SegNo]);
      fprintf(Animatfp, "\nRay PtBx[%d][%d]\n %0.14f", NRAYS + 1, SegNo,
              PtBx[SegNo]);
      fprintf(Animatfp, "\nRay PtBy[%d][%d]\n %0.14f", NRAYS + 1, SegNo,
              PtBy[SegNo]);
      for (i = 0; i < NRAYS + 2; i++) {
        fprintf(Animatfp, "\nSeg %d Ray %d Init Burn Depth %0.14f", SegNo, i,
                RayBurnDepth[i][SegNo]);
        fprintf(Animatfp, "\nSeg %d Ray %d Init Burn Pt x %0.14f", SegNo, i,
                BurnPtx[i][SegNo]);
        fprintf(Animatfp, "\nSeg %d Ray %d Init Burn Pt y %0.14f", SegNo, i,
                BurnPty[i][SegNo]);
      }
      for (i = 0; i < NRAYS + 2; i++) {
        if (SpecLayerOffset[i][SegNo] != 0.0) {
          fprintf(Animatfp,
                  "\nDepth at which Ray [%d] in seg [%d] STOPS burning fastest "
                  "'Catchup' propellant\n %0.14f",
                  i, SegNo, LatestCatchupDepthBurned[i][SegNo]);
        }
      }
    }
#endif /* of STARGRAIN */
    /* Then dims for nozzle */
    for (i = 0; i < NNOZRAYS; i++) {
      fprintf(Animatfp, "\nNozRayAngle[%d]\n %0.14f", i, NozRayAngle[i]);
    }
    for (i = 0; i < NNOZRAYS; i++) {
      fprintf(Animatfp, "\nNozRay special fast propellant type: 6");
      fprintf(Animatfp, "\nNoz End of Fast Burn Pts[%d]\n %0.14f", i,
              NozEndFastBurnPts[i]);
      fprintf(Animatfp, "\nNozRay special slow propellant type: 0");
    }
    for (i = 0; i < NNOZRAYS; i++) {
      fprintf(Animatfp, "\nNozRayFinalBurnRadius[%d]\n %0.14f", i,
              NozRayFinalBurnRadius[i]);
    }
    /* Then dims for insulation */
    fprintf(Animatfp, "\nTip of insulation, x = \n%0.14f", InsulPtx);
    fprintf(Animatfp, "\nTip of insulation, y = \n%0.14f", InsulPty);
    fprintf(Animatfp, "\nShell end horizontal insulation line, x = \n%0.14f",
            InsulShellPtx);
    fprintf(Animatfp, "\nShell end horizontal insulation line, y = \n%0.14f",
            InsulShellPty);
    fprintf(Animatfp, "\nTop of Insulation meets shell at x = %0.14f",
            InsulTopPtAtShellx);
    fprintf(Animatfp, "\nTop of Insulation meets shell at y = %0.14f",
            InsulTopPtAtShelly);
    fprintf(Animatfp, "\nFurthest Burned Pt along top of insul, x = \n%0.14f",
            InsulPropBdaryPtx);
    fprintf(Animatfp, "\nFurthest Burned Pt along top of insul, y = \n%0.14f",
            InsulPropBdaryPty);
    fprintf(Animatfp, "\nInsulTopBurningLength\n %0.14f",
            InsulTopBurningLength);
    fprintf(
        Animatfp,
        "\nInsulHalfAngle along which top insulation burn is tracked:\n%0.14f",
        InsulHalfAngle);

    /* Now print layer start radii and propellanttypes for ALL segments, ALL
     * layers. */
    for (SegNo = 0; SegNo < NTOTALSEGS; SegNo++) {
      for (i = 0; i < NLAYERS; i++) {
        fprintf(Animatfp, "\nLayerStartRadius[%d][%d]\n %0.14f", i, SegNo,
                LayerStartRadius[i][SegNo]);
      }
      for (i = 0; i < NLAYERS; i++) {
        fprintf(Animatfp, "\nPropellType[%d][%d] %d", i, SegNo,
                PropellType[i][SegNo]);
      }
    }
#ifdef STARGRAIN
    for (SegNo = 2; SegNo < NCYLSEGS - 1; SegNo++) {
      fprintf(Animatfp, "\nPreCylStarRDotRef[0] (outer) Segment[%d] %0.14f",
              SegNo, PreCylStarRDotRef[0][SegNo]);
      fprintf(Animatfp, "\nPreCylStarRDotRef[1] (inner) Segment[%d] %0.14f",
              SegNo, PreCylStarRDotRef[1][SegNo]);
      fprintf(Animatfp, "\nRDotRefCatchupInput[%d] %d", SegNo,
              RDotRefCatchupInput[SegNo - 2]);
    }
#endif /* of STARGRAIN */
    fprintf(Animatfp, "\nThat concludes the information defining the initial "
                      "configuration of the rocket!!!\n");
    fprintf(Animatfp, "\nThe information below is the state after each "
                      "timestep of burning for this rocket.");
  }
  fprintf(Animatfp, "\n\n### ### Time = %0.14f\n", Time);
  for (SegNo = 0; SegNo < NTOTALSEGS; SegNo++) {
    /* Note that CYLSEGS will not still be burning rays when any */
    /* seg burns out, due to fastest/slowest ref burn rate ratio */
    if (SegNo != NTOTALSEGS - 1) {
#ifdef STARGRAIN
      if (StarBurningThisSegNow[SegNo]) {
        fprintf(Animatfp, "\nStar Burning seg [%d]", SegNo);
        for (k = 0; k < NRAYS + 2; k++) {
          fprintf(Animatfp,
                  "\nFinalBurnRadius[%d] %0.14f, RayBurnDepth[%d][%d] %0.14f",
                  SegNo, FinalBurnRadius[SegNo], k, SegNo,
                  RayBurnDepth[k][SegNo]);
        }
      } else {
        fprintf(Animatfp,
                "\nFinalBurnRadius[%d] %0.14f, SegmentRadius[%d] %0.14f", SegNo,
                FinalBurnRadius[SegNo], SegNo, SegmentRadius[SegNo]);
      }
#else
      fprintf(Animatfp,
              "\nFinalBurnRadius[%d] %0.14f, SegmentRadius[%d] %0.14f", SegNo,
              FinalBurnRadius[SegNo], SegNo, SegmentRadius[SegNo]);
#endif
    }
    /* Nozzle seg could still be burning cylindrically, or as rays */
    else { /* Nozzle seg */
           // Put in prints of nozzle seg or ray burns
      if (NozBurningRays) {
        for (k = 0; k < NNOZRAYS; k++) {
          fprintf(Animatfp,
                  "\nNozBurningRays,  NozRayFinalBurnRadius[%d] %0.14f, "
                  "NozRayRadius[%d] %0.14f",
                  k, NozRayFinalBurnRadius[k], k, NozRayRadius[k]);
        }
      } else {
        fprintf(Animatfp,
                "\nFinalBurnRadius[%d] %0.14f, SegmentRadius[%d] %0.14f", SegNo,
                FinalBurnRadius[SegNo], SegNo, SegmentRadius[SegNo]);
      }
    }
  }
  return 1;
}

void FindNozzleRayPropellants(void) {
  int i;
  int CurrNozLayer;
  int Debug = 0;

  if (Debug)
    fprintf(outfp,
            "\n*** *** *** gen %d, neval %ld, Nozzle is Burning Rays!!!!!", gen,
            neval);
  /* In Nozzle segment after conversion to burning rays at boundary of last
   * zone. */
  if (NozRayRadius[8] >= NozEndFastBurnPts[8]) {
    /* NozRay 8 is the longest, last to equalize, so when this is true, we are
     */
    /* no longer equalizing rays, they are all same dist from outer shell. */
    /* So impose the same evolved propellant type on them all, so they burn out
     * together. */
    /* Find out what layer ray[0] is in, then give that ray's evolved propellant
     * type */
    /* to ALL the nozzle rays. */
    CurrNozLayer = NLAYERS - NLAYERSPERZONE;
    while (NozRayRadius[0] < LayerStartRadius[CurrNozLayer][NTOTALSEGS - 1]) {
      CurrNozLayer++;
    }
    for (i = 0; i < NNOZRAYS; i++) {
      NozRayRefBurnRate[i] =
          PropellChoices[PropellType[CurrNozLayer][NTOTALSEGS - 1]];
    }
    if (Debug) {
      fprintf(outfp,
              "\n*** Nozzle rays now equal, burning evolved propell type in "
              "layer %d, prop type %d",
              CurrNozLayer, PropellType[CurrNozLayer][NTOTALSEGS - 1]);
      fprintf(outfp, "\nRefBurnRate = %0.14f",
              PropellChoices[PropellType[CurrNozLayer][NTOTALSEGS - 1]]);
    }
  }
  for (i = 0; i < NNOZRAYS; i++) {
    if (NozRayRadius[i] >= NozRayFinalBurnRadius[i]) {
      NozRayRefBurnRate[i] = PropellChoices[11]; /* Namely, 0. */
    } else if (NozRayRadius[i] < NozEndFastBurnPts[i]) {
      NozRayRefBurnRate[i] = PropellChoices[6];
    } else {
      NozRayRefBurnRate[i] =
          PropellChoices[0]; /* i.e., distances not yet equalized */
    }
    if (Debug) {
      fprintf(outfp, "\nNozEndFastBurnPts[%d] %0.14f", i, NozEndFastBurnPts[i]);
      fprintf(outfp, "\nNozRayRefBurnRate[%d] = %0.14f", i,
              NozRayRefBurnRate[i]);
      fprintf(outfp, "\nBEFORE update, NozRayRadius[%d] = %0.14f", i,
              NozRayRadius[i]);
    }
  }
}

void SetupNozzleGeometry(void) {
  int Debug = 0;

  /* NOZZLE segment is segment number NTOTALSEGS-1, just after corner segs */
  NozRayBurnStartRadius =
      InitBurnRadius[1] +
      (FinalBurnRadius[1] - InitBurnRadius[1]) * (NZONES - 1) / (double)NZONES;
  if (Debug) {
    fprintf(outfp,
            "\nIn SetupNozzle Geom, FinalBurnRadius[1] %0.14f, "
            "InitBurnRadius[1] %0.14f",
            FinalBurnRadius[1], InitBurnRadius[1]);
    fprintf(outfp, "\nNozRayBurnStartRadius %0.14f", NozRayBurnStartRadius);
  }
  FinalBurnRadius[NTOTALSEGS - 1] =
      FinalBurnRadius[1]; /* set to match adjacent seg 1, */
  /* but not actually this radius; however, IS the radius to use to divide the
   */
  /* nozzle segment into zones,  as matches adjacent segment[1]. */
  SprayNozzleRays();
  MakeNozzleProfileRedefinitionLayer();
}

void MakeNozzleProfileRedefinitionLayer(void) {

  // NOTE:  THIS is for NOZZLE segment only!!!
  double BN, JQ, AK, AL, QS, LK, TQ, violetprofile[NNOZRAYS], burnratio;
  int i;
  int Debug = 0;

  /* The profile redefinition layer will consist of type 4 propellant, burning
   * at */
  /* 0.00521 m/sec, or about 2.05 times as fast as type 0, which is 0.00254 m/s.
   */
  /* The thickness of the faster-burning propellant on each ray will be
   * determined */
  /* from the offset profile of the outer shell, translated so equidistant from
   * that */
  /* shell everywhere, but with its distance on ray A (ray[0]) equal to the
   * thickness */
  /* of the final zone on that ray. */
  burnratio = PropellChoices[6] / PropellChoices[0];
  /* by equality of burn times: */
  BN = NozRayFinalBurnRadius[8] -
       (InitBurnRadius[1] +
        (FinalBurnRadius[1] - InitBurnRadius[1]) * (NZONES - 1) / NZONES);
  AK = (FinalBurnRadius[1] - InitBurnRadius[1]) /
       NZONES; /* one zone thickness */
  AL = (burnratio * AK - BN) / (burnratio - 1.);
  LK = AK - AL;
  if (Debug) {
    fprintf(outfp, "\nNozRayFinalBurnRadius[8] %0.14f",
            NozRayFinalBurnRadius[8]);
    fprintf(outfp,
            "\nburnratio %0.14f, BN %0.14f, AL %0.14f, AK %0.14f, LK  %0.14f, "
            "AL %0.14f\n",
            burnratio, BN, AL, AK, LK, AL);
  }
  NozEndFastBurnPts[0] = NozRayBurnStartRadius;
  for (i = 1; i < NNOZRAYS; i++) {
    violetprofile[i] = NozRayFinalBurnRadius[i] - AK;
    JQ = NozRayFinalBurnRadius[i] - NozRayBurnStartRadius;
    TQ = JQ - AL;
    QS = burnratio * (TQ - LK) / (burnratio - 1.);
    NozEndFastBurnPts[i] = violetprofile[i] + QS;
    if (Debug) {
      fprintf(outfp,
              "\nJQ %0.14f, TQ %0.14f, QS %0.14f, violetprofile[%d] %0.14f, "
              "NozEndFastBurnPts[%d] %0.14f",
              JQ, TQ, QS, i, violetprofile[i], i, NozEndFastBurnPts[i]);
    }
  }
}

void SprayNozzleRays(void) {
  int i;
  int Debug = 0;

  /* Nozzle will burn in rays after passing NozRayBurnStartRadius, for all of
   * last zone */
  /* Refer to PowerPoint diagram for ray definitions. */
  /* Total of 9 rays in following order: A, F, E, D, G, H, C, J, B (left to
   * right) */
  NozRayAngle[0] = PI / 2.;
  NozRayAngle[8] = 11. * PI / 180.;
  NozRayAngle[6] = +.58468; /* from spreadsheet */
  NozRayAngle[3] = 57. / 180 * PI;
  NozRayAngle[7] = (NozRayAngle[6] + NozRayAngle[8]) / 2.;
  NozRayAngle[5] = NozRayAngle[6] + (NozRayAngle[3] - NozRayAngle[6]) * 1. / 3.;
  NozRayAngle[4] = NozRayAngle[6] + (NozRayAngle[3] - NozRayAngle[6]) * 2. / 3.;
  NozRayAngle[1] = NozRayAngle[3] + (NozRayAngle[0] - NozRayAngle[3]) * 2. / 3.;
  NozRayAngle[2] = NozRayAngle[3] + (NozRayAngle[0] - NozRayAngle[3]) * 1. / 3.;
  NozShortestRayAngle = NozRayAngle[0];

  /* Now calculate Final Burn Radius for each ray, in same order as before. */
  NozRayFinalBurnRadius[0] = FinalBurnRadius[1];
  // Corrected 6/28/19
  NozRayFinalBurnRadius[8] = FinalBurnRadius[1] * 10. / 8.6;
  NozRayFinalBurnRadius[6] = FinalBurnRadius[1] * 9.1 / 8.6;
  NozRayFinalBurnRadius[3] = FinalBurnRadius[1] * 9.8 / 8.6;
  NozRayFinalBurnRadius[7] =
      (NozRayFinalBurnRadius[8] + NozRayFinalBurnRadius[6]) / 2.;
  /* last 4 have negligible error due to curvature near D, but not worth working
   * out. */
  NozRayFinalBurnRadius[5] =
      (2. * NozRayFinalBurnRadius[6] + NozRayFinalBurnRadius[3]) / 3.;
  NozRayFinalBurnRadius[4] =
      (NozRayFinalBurnRadius[6] + 2. * NozRayFinalBurnRadius[3]) / 3.;
  NozRayFinalBurnRadius[2] =
      (NozRayFinalBurnRadius[0] + 2. * NozRayFinalBurnRadius[3]) / 3.;
  NozRayFinalBurnRadius[1] =
      (2. * NozRayFinalBurnRadius[0] + NozRayFinalBurnRadius[3]) / 3.;
  if (Debug) {
    for (i = 0; i < NNOZRAYS; i++) {
      fprintf(outfp,
              "\nIn SprayNozzleRays, NozRayAngle[%d] %0.14f, "
              "NozRayFinalBurnRadius[%d] %0.14f",
              i, NozRayAngle[i], i, NozRayFinalBurnRadius[i]);
    }
  }
}

#ifdef STARGRAIN
void CalcStarRayBurnDirnsAndLengths(void) {
  int SegNo;

  /* Called during burn, not initialization */
  for (SegNo = 2; SegNo < NCYLSEGS - 1; SegNo++) {
    if (StarBurningThisSegNow[SegNo])
      CalcRayBurnDirnsAndLengths(SegNo);
  }
}

void CalcPerRayRDot(int SegNo) {
  int i;
  int Debug = 0;

  /* Here, adjust ref burn rates for actual pressure for the RAYS being burnt in
   * STAR burn cylinders */
  for (i = 0; i < NRAYS + 2; i++) {
    RayRDot[i][SegNo] =
        RayRDotRef[i][SegNo] * pow(TotalPressure, Alpha) * InvPrefAlpha;
    // for normal segs, is: RDot[SegNo] = RefBurnRate[SegNo] *
    // pow(TotalPressure, Alpha) * InvPrefAlpha;
    if (Debug)
      fprintf(outfp, "\nRayRDotRef[%d][%d] = %0.14e", i, SegNo,
              RayRDotRef[i][SegNo]);
    if (Debug)
      fprintf(outfp, "\nRayRDot[%d][%d] = %0.14e", i, SegNo, RayRDot[i][SegNo]);
  }
  return;
}
#endif

void CalcNozRayRDots(void) {
  int i;

  for (i = 0; i < NNOZRAYS; i++) {
    NozRayRDot[i] =
        NozRayRefBurnRate[i] * pow(TotalPressure, Alpha) * InvPrefAlpha;
  }
  return;
}

#ifdef STARGRAIN
void CalcStarBurnIntegrContribs(int SegNo) {
  int i;
  double LocRDotRef;

  int Debug = 0;

  /* Need to calculate the integral contribution of a segment full of rays (star
   * burn segment, first 3 layers). */
  /* Sum up the product of each ray's arc length (half of edge lengths,
   * center-to-center for the rest, */
  /* with its burn rate number. Add up the products, divide by total length to
   * get equivalent RDot. */
  /* NOTE: will update burn depths with trig functions later */
  /* First handle the two edge rays separately. Give each HALF of the width to
   * next ray. */
  SumOfProducts[SegNo] = 0.;
  if (Debug) {
    fprintf(outfp,
            "\nRayiToRayiminus1[0][%d] is not used... half length is of [1].",
            SegNo);
    for (i = 1; i < NRAYS + 2; i++) {
      fprintf(
          outfp,
          "\nIn CalcStarBurnIntegrContribs, RayiToRayiminus1[%d][%d] = %0.14f",
          i, SegNo, RayiToRayiminus1[i][SegNo]);
    }
  }
  for (i = 0; i < NRAYS + 2; i++) {
    if (Debug) {
      fprintf(outfp, "\nRayInLayer[%d][%d] = %d", i, SegNo,
              RayInLayer[i][SegNo]);
      if (RayInLayer[i][SegNo] == -1) {
        fprintf(outfp, "\nPreCylStarRDotRef[0][%d] = %0.14f", SegNo,
                PreCylStarRDotRef[0][SegNo]);
        fprintf(outfp, "\nPreCylStarRDotRef[1][%d] = %0.14f", SegNo,
                PreCylStarRDotRef[1][SegNo]);
      } else {
        fprintf(
            outfp,
            "\nRayBurnDepth[%d][%d] = %0.14f, LayerStartRadius[%d][%d] %0.14f",
            i, SegNo, RayBurnDepth[i][SegNo], RayInLayer[i][SegNo], SegNo,
            LayerStartRadius[i][SegNo]);
        fprintf(outfp, "\nNextLayer is %d, its LayerStartRadius is %0.14f",
                RayInLayer[i][SegNo] + 1,
                LayerStartRadius[RayInLayer[i][SegNo]][SegNo]);
      }
    }
    if (RayInLayer[i][SegNo] == -1) {
      /* Modification to allow for 2 layers in PreCylStar regions */
      if (RayBurnDepth[i][SegNo] < StarRayOuterInnerBdary[SegNo]) {
        LocRDotRef = PreCylStarRDotRef[0][SegNo];
      } else {
        LocRDotRef = PreCylStarRDotRef[1][SegNo];
      }
      /* -1 is a flag saying are burning in star part, not a cylinder layer */
      if (i == 0) {
        SumOfProducts[SegNo] +=
            ((RayiToRayiminus1[1][SegNo]) / 2.) * LocRDotRef;
        if (Debug)
          fprintf(outfp,
                  "\n Segment %d's Ray %d is in Star point propellant rate "
                  "%0.14f\n",
                  SegNo, i, LocRDotRef);
      } else if (i == NRAYS + 1) {
        SumOfProducts[SegNo] +=
            ((RayiToRayiminus1[NRAYS + 1][SegNo]) / 2.) * LocRDotRef;
        if (Debug)
          fprintf(outfp,
                  "\n Segment %d's Ray %d is in Star point propellant rate "
                  "%0.14f\n",
                  SegNo, i, LocRDotRef);
      } else {
        SumOfProducts[SegNo] +=
            ((RayiToRayiminus1[i][SegNo] + RayiToRayiminus1[i - 1][SegNo]) /
             2.) *
            LocRDotRef;
        if (Debug)
          fprintf(outfp,
                  "\n Segment %d's Ray %d is in Star point propellant rate "
                  "%0.14f\n",
                  SegNo, i, LocRDotRef);
      }
      if (Debug) {
        fprintf(
            outfp,
            "\nIn CalcStarBurnIntegrContribs Layer -1, SumOfProducts = %0.14f",
            SumOfProducts[SegNo]);
      }
    } else {
      // HERE correcting bug found by Abhiroop, 9/17/19
      if (i == 0) {
        SumOfProducts[SegNo] +=
            RayiToRayiminus1[i][SegNo] / 2. * RayRDotRef[i][SegNo];
        if (Debug) {
          fprintf(outfp,
                  "\nIn CalcStarBurnIntegrContribs Seg %d, Ray 0, "
                  "SumOfProducts = %0.14f",
                  SegNo, SumOfProducts[SegNo]);
        }
      } else if (i == NRAYS + 1) {
        SumOfProducts[SegNo] +=
            RayiToRayiminus1[NRAYS + 1][SegNo] / 2. * RayRDotRef[i][SegNo];
        if (Debug)
          fprintf(outfp,
                  "\nIn CalcStarBurnIntegrContribs Last Ray, SumOfProducts = "
                  "%0.14f",
                  SumOfProducts[SegNo]);
      } else {
        SumOfProducts[SegNo] +=
            (RayiToRayiminus1[i][SegNo] + RayiToRayiminus1[i - 1][SegNo]) / 2. *
            RayRDotRef[i][SegNo];
        if (Debug) {
          fprintf(
              outfp,
              "\nIn CalcStarBurnIntegrContribs Ray %d, SumOfProducts = %0.14f",
              i, SumOfProducts[SegNo]);
          fprintf(outfp,
                  "\n Segment %d's Ray %d is burning RayRDotRef = %0.14f\n",
                  SegNo, i, RayRDotRef[i][SegNo]);
        }
      }
    }
    if (Debug) {
      fprintf(outfp,
              "\n in Star seg within cylinder, after Ray[%d], "
              "SumOfProducts[%d] = %0.14f",
              i, SegNo, SumOfProducts[SegNo]);
    }
  }
  SumOfProducts[SegNo] *=
      SegmentLength[SegNo]; /* Length is along vertical axis of rocket. */
  SumOfProducts[SegNo] *=
      NSYMS; /* To adjust for the fraction of circle simulated. */
  if (Debug) {
    fprintf(
        outfp,
        "\n*** After mult by seg length, Final SumOfProducts[%d] = %0.14f\n",
        SegNo, SumOfProducts[SegNo]);
  }
}
#endif

void UpdateInsulationBurn(void) {
  int Debug = 0;

  if (!NozBurningRays) {
    InsulPropBdaryPtx =
        cos(11. / 180. * PI) * (SegmentRadius[12] - InitBurnRadius[12]);
    InsulPropBdaryPty =
        InitBurnRadius[12] +
        sin(11. / 180. * PI) * (SegmentRadius[12] - InitBurnRadius[12]);
  } else {
    // 6/27/19 corrected to what was originally intended
    InsulPropBdaryPtx =
        cos(11. / 180. * PI) * (NozRayRadius[8] - InitBurnRadius[12]);
    InsulPropBdaryPty =
        InitBurnRadius[12] +
        sin(11. / 180. * PI) * (NozRayRadius[8] - InitBurnRadius[12]);
    if (Debug)
      fprintf(outfp,
              "\nNozRayFinalBurnRadius[8] %0.14f, InsulPropBdaryPtx %0.14f, "
              "InsulPropBdaryPty %0.14f",
              NozRayRadius[8], InsulPropBdaryPtx, InsulPropBdaryPty);
  }
  InsulDeltaBottomDepth =
      InsulRefBurnRate * pow(TotalPressure / RefPressure, 0.8) * DeltaTime;
  InsulBottomBurnDepth += InsulDeltaBottomDepth;
  if (Debug) {
    fprintf(outfp,
            "\nUpdating BurnRadii: InsulRefBurnRate %0.14e, TotalPressure "
            "%0.14e, RefPressure %0.14e",
            InsulRefBurnRate, TotalPressure, RefPressure);
    fprintf(outfp, "\npow(TotalPressure/RefPressure, 0.8) %0.14e",
            pow(TotalPressure / RefPressure, 0.8));
    fprintf(outfp,
            "\nInsulDeltaBottomDepth %0.14e, InsulPropBdaryPtx %0.14e, "
            "InsulPropBdaryPty %0.14e",
            InsulDeltaBottomDepth, InsulPropBdaryPtx, InsulPropBdaryPty);
    fprintf(outfp, "\nInsulBottomBurnDepth %0.14e", InsulBottomBurnDepth);
  }
  InsulPointBurnRadius = InsulBottomBurnDepth / sin(InsulHalfAngle);
  InsulPtx = InsulPointBurnRadius * cos(InsulHalfAngle);
  InsulPty = InitBurnRadius[1] + InsulPointBurnRadius * sin(InsulHalfAngle);
  if (Debug) {
    fprintf(outfp,
            "\nIn UpdateInsulationBurn, InsulPointBurnRadius %0.14f, InsulPtx "
            "%0.14f, InsulPty %0.14f",
            InsulPointBurnRadius, InsulPtx, InsulPty);
  }
  if (InsulBottomBurnDepth > 0.0) {
    // 6/26/19 Corrected below, was using wrong variable.
    InsulShellPtx -= InsulDeltaBottomDepth / tan((56.5 / 180.) * PI);
  }
  InsulShellPty = InitBurnRadius[1] + InsulBottomBurnDepth;
  if (Debug) {
    fprintf(outfp, "\nInsulBottomBurnDepth %0.14e, tan((56.5/180.)*PI) %0.14e",
            InsulBottomBurnDepth, tan((56.5 / 180.) * PI));
  }
  /* Here updating for new bottom burn depth */
  InsulTopBurningLength =
      DistBetwPts(InsulPropBdaryPtx, InsulPropBdaryPty, InsulPtx, InsulPty);
  /* As long as driving the InsulBurnRadius from the bottom burn, we COULD
   * calculate these below, */
  /* but we will NOT use them to calculate the new burn front... we are using
   * the new */
  /* InsulPointBurnRadius coords and the SegmentRadius[1], since the seg 1
   * propellant */
  /* will burn much faster and determine where the distant burn point is. */
  /* BUT WE MUST capture the top insulation burned in the residue() calculation,
   */
  /* by including the insul top area in the total insul area being burned. */
  /* Burn top portion of insulation that is exposed */
  /* Calculate normal burn direction */

  /* I don't need to update anything now with this info... it will be used to
   * calc burn area contrib, */
  /* but does not need updating here. */
  InsulTopSlope =
      (InsulPropBdaryPty - InsulPty) / (InsulPropBdaryPtx - InsulPtx);
  InsulTopBurnLineAngle = atan(InsulTopSlope);
  InsulTopBurnNormal = InsulTopBurnLineAngle - PI / 2.;
}

#ifdef STARGRAIN
void UpdateStarRayBurns(int SegNo) {
  double Angle, Angle1, Angle2;
  double DepthIncrease, DepthIncrease1, DepthIncrease2;
  double AvgBurnDir;
  int i;

  int Debug = 0;

  /* Doing RAY burns in a STAR segment, still in star part. */
  /* Must update burn depth on each radial RAY, instead of setting NextSegRadius
   * below */
  /* Ray[0] and Ray[NRAYS+1], because of the symmetry boundaries, are handled
   * separately . */
  /* For [0], use burn line between depths on 0 and 1 to calc burn dirn, but
   * calc the offset */
  /* into adjacent seg, then back that off to a new depth on Ray[0]. */
  /* If a concave corner (i.e., RayAngle[0] is positive), just advances along
   * ray by */
  /* amount of burn depth change. */

  if (Time == 0. && FirstRayCallThisGuy[SegNo]) {
    RayDebugPrint(SegNo);
    FirstRayCallThisGuy[SegNo] = 0;
  }

  /* HERE IS THE GUTS OF THE RAY BURN MODEL!!! */
  /* First handle Ray 0 */
  Angle = RayAngle[0][SegNo] - RayBurnDir[0][SegNo];
  DepthIncrease = (RayRDot[0][SegNo] * DeltaTime) / cos(Angle);
  if (Debug)
    fprintf(outfp,
            "\n\n\nfor Ray 0, Angle %0.14f, RayRDot*DT %0.14f, DepthIncrease "
            "%0.14f, RayRDot[0][%d] %0.14f, cos(Angle) %0.14f",
            Angle, RayRDot[0][SegNo] * DeltaTime, DepthIncrease, SegNo,
            RayRDot[0][SegNo], cos(Angle));
  if (RayBurnDir[0][SegNo] > 0. &&
      DepthIncrease > RayRDot[0][SegNo] * DeltaTime) {
    /* This constrains rate of advance of concave corner to be the burn rate
     * along the */
    /* bisector of the angle (representing normal to a small arc there), which
     * is the */
    /* ray itself, which reflects the actual upper limit on that rate of
     * advance. */
    DepthIncrease = RayRDot[0][SegNo] * DeltaTime;
  }
  if (Debug) {
    fprintf(outfp, "\n!!! Ray [0][%d] Angle %0.14e, DepthIncrease %0.14e",
            SegNo, Angle, DepthIncrease);
  }
  RayBurnDepth[0][SegNo] += DepthIncrease;
  BurnPtx[0][SegNo] = RayBurnDepth[0][SegNo] * cos(RayAngle[0][SegNo]);
  BurnPty[0][SegNo] = RayBurnDepth[0][SegNo] * sin(RayAngle[0][SegNo]);

  for (i = 1; i < NRAYS + 2; i++) {
    if (i == NRAYS + 1) {
      /* First handle final ray, which actually reflects back onto itself, in
       * what CANNOT be a convex corner */
      /* since no other ray can be deeper than Point B, Ray[NRAYS+1]. So this
       * ray WILL always be burning only */
      /* along the ray itself, since that is also the bisector, in this case. */
      RayBurnDepth[NRAYS + 1][SegNo] += RayRDot[NRAYS + 1][SegNo] * DeltaTime;
      BurnPtx[NRAYS + 1][SegNo] =
          RayBurnDepth[NRAYS + 1][SegNo] * cos(RayAngle[NRAYS + 1][SegNo]);
      BurnPty[NRAYS + 1][SegNo] =
          RayBurnDepth[NRAYS + 1][SegNo] * sin(RayAngle[NRAYS + 1][SegNo]);
    } else {
      /* Handling INTERNAL rays */

      // EDG 031619 REVISED RayBurn model.
      //
      Angle1 = RayAngle[i][SegNo] - RayBurnDir[i][SegNo];
      if (Angle1 > 0.) {
        DepthIncrease1 = RayRDot[i][SegNo] * DeltaTime;
      } else {
        DepthIncrease1 = (RayRDot[i][SegNo] * DeltaTime) / cos(Angle1);
      }
      Angle2 = RayAngle[i][SegNo] - RayBurnDir[i + 1][SegNo];
      if (Angle2 <= 0.) {
        DepthIncrease2 = RayRDot[i][SegNo] * DeltaTime;
      } else {
        DepthIncrease2 = (RayRDot[i][SegNo] * DeltaTime) / cos(Angle2);
      }
      if (DepthIncrease1 > DepthIncrease2) {
        DepthIncrease = DepthIncrease1;
      } else {
        DepthIncrease = DepthIncrease2;
      }

      RayBurnDepth[i][SegNo] += DepthIncrease;
      BurnPtx[i][SegNo] = RayBurnDepth[i][SegNo] * cos(RayAngle[i][SegNo]);
      BurnPty[i][SegNo] = RayBurnDepth[i][SegNo] * sin(RayAngle[i][SegNo]);
      if (Debug) {
        fprintf(outfp,
                "\n!!! RayRDot[%d][%d] %0.14f, RayAngle[%d][%d] %0.14f, "
                "RayBurnDir[%d][%d] %0.14f",
                i, SegNo, RayRDot[i][SegNo], i, SegNo, RayAngle[i][SegNo], i,
                SegNo, RayBurnDir[i][SegNo]);
        fprintf(outfp,
                "\n!!! HERE WE ARE: Angle betw ray and incoming burn dir [%d] "
                "%0.14f",
                i, Angle1);
        fprintf(outfp,
                "\n!!! HERE WE ARE: Angle betw ray and outgoing burn dir [%d] "
                "%0.14f",
                i + 1, Angle2);
        fprintf(outfp,
                "\n!!! DepthIncrease1 %0.14f, DepthIncrease2, %0.14f, "
                "DepthIncrease %0.14f",
                DepthIncrease1, DepthIncrease2, DepthIncrease);
        fprintf(outfp,
                "\nIn UpdateStarRayBurns, RayBurnDepth[%d][%d] = %0.14f\n", i,
                SegNo, RayBurnDepth[i][SegNo]);
        fprintf(outfp,
                "\nIn UpdateStarRayBurns, RayAngle %0.14f, Burnptx[%d][%d] = "
                "%0.14f, BurnPty[%d][%d] = %0.14f",
                RayAngle[i][SegNo], i, SegNo, BurnPtx[i][SegNo], i, SegNo,
                BurnPty[i][SegNo]);
      }
    }
  }
  RayDebugPrint(SegNo);
}

void RayDebugPrint(int SegNo) {
  int i;

  /* REWORKING all the random, in situ, RAY debug prints into one comprehensive
   * print to allow comparison of */
  /* all relevant quantities side by side, and from timestep to timestep. */
  if (RayDebug) {
    if (FirstRayCallThisGuy[SegNo]) {
      fprintf(outfp, "\n Time: INITIALIZED at Time 0.0, Segment [%d]", SegNo);
    } else {
      fprintf(outfp, "\nneval %ld, AFTER burn at Time %0.14f, Segment [%d]",
              neval, Time, SegNo);
    }
    fprintf(outfp, "\nNo, Ray Burn Depth,    BurnPtx, BurnPty, RayAngle, "
                   "BurnDir,   Layer, NextBdary");

    fprintf(outfp, "\nStarRayOuterInnerBdary[%d] %0.14f", SegNo,
            StarRayOuterInnerBdary[SegNo]);

    for (i = 0; i < NRAYS + 2; i++) {
      if (RayInLayer[i][SegNo] == -1) {
        fprintf(outfp,
                "\n%d %0.14e %0.14f %0.14f %0.14f  %0.14f   %d    %0.14f", i,
                RayBurnDepth[i][SegNo], BurnPtx[i][SegNo], BurnPty[i][SegNo],
                RayAngle[i][SegNo], RayBurnDir[i][SegNo], RayInLayer[i][SegNo],
                LayerStartRadius[0][SegNo]);
        if (RayBurnDepth[i][SegNo] < StarRayOuterInnerBdary[SegNo])
          fprintf(outfp, "\nPreCylStarRDotRef[0][%d] %0.14f", SegNo,
                  PreCylStarRDotRef[0][SegNo]);
        else
          fprintf(outfp, "\nPreCylStarRDotRef[1][%d] %0.14f", SegNo,
                  PreCylStarRDotRef[1][SegNo]);

      } else if (!CompensatingNow[SegNo]) {
        fprintf(outfp,
                "\n%d %0.14e %0.14f %0.14f %0.14f  %0.14f   %d    %0.14f   "
                "Prop in layer %0.14f",
                i, RayBurnDepth[i][SegNo], BurnPtx[i][SegNo], BurnPty[i][SegNo],
                RayAngle[i][SegNo], RayBurnDir[i][SegNo], RayInLayer[i][SegNo],
                LayerStartRadius[RayInLayer[i][SegNo] + 1][SegNo],
                RDotRef[RayInLayer[i][SegNo]][SegNo]);
      } else {
        /* Compensating, SOME rays burning catchup fuel; print it, too */
        fprintf(outfp,
                "\n%d %0.14e %0.14f %0.14f %0.14f  %0.14f   %d    %0.14f    "
                "(see below)",
                i, RayBurnDepth[i][SegNo], BurnPtx[i][SegNo], BurnPty[i][SegNo],
                RayAngle[i][SegNo], RayBurnDir[i][SegNo], RayInLayer[i][SegNo],
                LayerStartRadius[RayInLayer[i][SegNo] + 1][SegNo]);
        if (RayBurnDepth[i][SegNo] <
                RayBurnDepth[ShortestRayNo[SegNo]][SegNo] &&
            RayInLayer[i][SegNo] > StartCatchupLayer[SegNo] &&
            CatchingUpRay[i][SegNo] < 2) {
          fprintf(
              outfp,
              "\nRay %d burning Catchup propellant RDotRefCatchup[%d] = %0.14f",
              i, SegNo, RayRDotRef[i][SegNo]);
        }
        fprintf(outfp,
                "\nRegardless, RayRDotRef being burned now in Ray %d is %0.14f",
                i, RayRDotRef[i][SegNo]);
      }
    }
  }
}
#endif

int UpdateBurnRadiiAndReachedShell(void) {
  int SegNo, i;

  int Debug = 0;

  /* At this point, we KNOW the actual pressure, so can scale the RDotRefs using
   * that to get RDots */
  /* The key update to the burn radius for this segment */
  /* Because of burn overlaps at seg boundaries, we do not update SegmentRadius
   * until all burns for */
  /* this step have been calculated */
  UpdateInsulationBurn();

  for (SegNo = 0; SegNo < NTOTALSEGS; SegNo++) {
#ifdef STARGRAIN
    if (SegNo > 1 && SegNo < NCYLSEGS - 1 && StarBurningThisSegNow[SegNo]) {
      /* Doing ray burns in a star segment */
      UpdateStarRayBurns(SegNo);
    } else if (SegNo == NTOTALSEGS - 1 && NozBurningRays) {
#else
    if (SegNo == NTOTALSEGS - 1 && NozBurningRays) {
#endif
      CalcNozRayRDots();
      if (Debug)
        fprintf(outfp, "\nBurning Nozzle Rays.");
      for (i = 0; i < NNOZRAYS; i++) {
        NextNozRayRadius[i] = NozRayRadius[i] + NozRayRDot[i] * DeltaTime;
        // WORKING HERE
        if (NextNozRayRadius[i] > NozRayFinalBurnRadius[i]) {
          WhenBurnHitsShell(SegNo);
          return (1);
        }
      }
    } else {
      NextSegRadius[SegNo] += RDot[SegNo] * DeltaTime;
      if (NextSegRadius[SegNo] > FinalBurnRadius[SegNo]) {
        WhenBurnHitsShell(SegNo);
        return (1);
      }
    }
  }
  // SPECIAL CORNER SEG BURN INITIALIZATION
  /* NOW for very first timestep of CORNER segments, their area started at 0, so
   * they will not burn at all!  But */
  /* they would burn exactly the same amounts as cyl segs they are adjacent to,
   * so just copy that amount in as their */
  /* change in NextSegRadius[SegNo] for the FIRST TIMESTEP ONLY. */
  if (Time == 0.0) {
    for (SegNo = 0; SegNo < NCORNERSEGS; SegNo++) {
      NextSegRadius[SegNo + NCYLSEGS] = NextSegRadius[NCYLSEGS - 1];
    }
  }
  /* EXACTLY the same for the initial NOZZLE burn, initial area is 0 at the
   * point, so do first step */
  /* manually here, to match adjacent segment, which is correct, same propellant
   * type. */
  if (Time == 0.0) {
    NextSegRadius[NCYLSEGS - 1] = NextSegRadius[1];
  }
  return (0);
}

void WhenBurnHitsShell(int SegNo) {

  int Debug = 0;

  /* Take these actions when burn hits shell */
  if (Debug) {
    fprintf(outfp,
            "\nSegment %d Burn reached shell, RDot[%d] = %0.14f, Thrust "
            "%0.14f, NextSegRadius[%d] %0.14f, Time %0.14f\n",
            SegNo, SegNo, RDot[SegNo], Thrust, SegNo, NextSegRadius[SegNo],
            Time);
  }
  /* ASSIGN FITNESS, RETURN THIS INDIVIDUAL */
  Burning = (Time + DeltaTime) < TotalBurnTime;
  /* Reward for completing burn without blowing up */
  Reward += 20000.;
  Fitness = calc_fitness();
  if (Debug) {
    fprintf(outfp,
            "\n *** Outer wall reached, neval %ld, gen %d. Fitness = %0.14f\n",
            neval, gen, Fitness);
  }
  StopCode = 0;
}

void RecordThrustProfile(int irecord) {
  int Debug = 0;

  Times[irecord] = Time;
  ThrustProfile[irecord] = Thrust;
  mPressureProfile[irecord] = TotalPressure;
  if (Debug) {
    fprintf(outfp, "\nTime %0.14f, ThrustProfile %0.14f, Pressure %0.14f\n",
            Times[irecord], ThrustProfile[irecord], mPressureProfile[irecord]);
  }
}

// Added by Abhiroop 09/24/2019. Used to record thrust profile on a much finer
// scale
void RecordThrustProfileFine(int irecord) {
  int Debug = 0;

  TimesFine[irecord] = Time;
  ThrustFine[irecord] = Thrust;
  mPressureFine[irecord] = TotalPressure;
  if (Debug) {
    fprintf(outfp, "\nTime %0.14f, ThrustProfile %0.14f, Pressure %0.14f\n",
            TimesFine[irecord], ThrustFine[irecord], mPressureFine[irecord]);
  }
}

void MakeCornerSegBdaries(double *BdaryAngle, double *c, double *d) {
  /* calculate thetabig, angle between boundaries of cylindrical segment
   * [NCYLSEGS-1] and central dome segment. */
  /* Then divide it to form boundaries of dome segments. */
  int SegNo;
  double n, adjacent;
  int Debug = 0;

  if (Debug)
    fprintf(outfp,
            "\nPy %0.14f, Qy %0.14f, Px %0.14f, Qx %0.14f, "
            "InitBurnRadius[NCYLSEGS-1] %0.14f\n",
            Py, Qy, Px, Qx, InitBurnRadius[NCYLSEGS - 1]);
  /* NOTE:  index [NCYLSEGS-1] refers to the LAST CYLINDRICAL LAYER. [0] is the
   * central dome */

  bigtheta = asin((Py - Qy) / DistBetwPts(Px, Py, Qx, Qy));
  // pow((Px-Qx)*(Px-Qx) + (Py - Qy) * (Py - Qy), 0.5));
  ThetaSlice = bigtheta / (NCORNERSEGS);
  BdaryAngle[0] = 0.;
  BdaryAngle[NCORNERSEGS] = bigtheta;
  for (SegNo = 1; SegNo < NCORNERSEGS + 1; SegNo++) {
    BdaryAngle[SegNo] = BdaryAngle[SegNo - 1] + ThetaSlice;
  }

  /* Now for each angle above the initial 0 angle, calculate point where it
   * intercepts either */
  /* the cylinder (lower points) or ellipse (upper points) as coordinates (c,d),
   * with line */
  /* emanating from point (m,n), the "corner" point where ellipse web hits
   * cylinder web. */
  /* The cylinder intersections are simple, y = InitBurnRadius * tan(theta), but
   * the higher points */
  /* require an elliptical integral, so must iterate to find it.  Guess a d;
   * then */
  /* try to minimize the elliptical eqn -1 (so it goes over and under 0.) */
  adjacent = (RocketCylRadius - InitBurnRadius[0]);
  n = Qy;
  ellipsec[NCORNERSEGS] = Px;
  ellipsed[NCORNERSEGS] = Py;
  for (SegNo = 0; SegNo < NCORNERSEGS; SegNo++) {
    if (adjacent * tan(BdaryAngle[SegNo]) <= fabs(n) + 0.0000001) {
      ellipsed[SegNo] = n + adjacent * tan(BdaryAngle[SegNo]);
      ellipsec[SegNo] = RocketCylRadius;
    } else {
      /* ellipse intersection case */
      /* Will get ellipsec[i] returned, not sending anything */
      THIS_ANGLE = BdaryAngle[SegNo];
      ellipsed[SegNo] = bisection(&fnEllipse, &(ellipsec[SegNo]), 0., .0253);
      if (Debug) {
        fprintf(outfp,
                "\n**** InMakeCornerSegBdaries back from bisection, "
                "ellipsec[%d] %0.14f, ellipsed[%d] %0.14f\n",
                SegNo, ellipsec[SegNo], SegNo, ellipsed[SegNo]);
      }
    }
  }
}

double fnEllipse(double guessd, double *this_c) {
  /* Calculate ellipse equation - 1 so goes above and below 0. */
  /* d passed in is a guess at a d that will satisfy the (c,d) point
   * intersecting the ellipse */
  /* with a line from (m,n) ("corner" of web) at angle THIS_ANGLE. */
  double temptozero, sqrtterm, d_estim;
  double n = Qy;
  double m = InitBurnRadius[NCYLSEGS - 1];

  *this_c = RocketCylRadius *
            pow(1. - guessd * guessd / (SemiMinAx * SemiMinAx), 0.5);

  sqrtterm = pow(1. - ((guessd * guessd) / (SemiMinAx * SemiMinAx)), 0.5);
  d_estim = n + (RocketCylRadius * sqrtterm - m) * tan(THIS_ANGLE);
  temptozero = d_estim - guessd;
  return (temptozero); /* sending back the difference from ideal ellipse */
}

double bisection(double (*fntocall)(double param1, double *param2),
                 double *secondreturn, double lowerguess, double upperguess) {
  /* Should work for any one or two args passed in as pointers, so long as
   * define */
  /* a function that is positive at one bound and negative at other and doesn't
   * have */
  /* multiple internal zeroes. */

  /* Function to find zero of is defined in the function called. */
  double answer = 100.;
  double threshold = 0.00001; /* Should reach before quitting from count */
  double guess, lower, upper;
  double temp;
  int counter = 0;

  int Debug = 0;

  lower = lowerguess;
  upper = upperguess;
  /* Are iterating on values of Px as the guesses */
  /* Count of 30 guarantees accuracy of answer to a millionth of initial
   * interval. */
  while (fabs(answer) > threshold && counter < 30) {
    guess = (lower + upper) / 2.;
    /* Call function to zero */
    answer = fntocall(guess, secondreturn);
    counter++;
    if (Debug)
      fprintf(outfp,
              "\nIn bisection, lower %0.14f, upper %0.14f, answer %0.14f",
              lower, upper, answer);
    /* Second call to function to zero */
    temp = fntocall(lower, secondreturn);
    if (temp * answer < 0.) {
      upper = guess;
    } else {
      lower = guess;
    }
  }
  if (counter >= 30) {
    printf("\n*** OOPS *** Bisection failed to find a root, quitting.\n");
    exit(999);
  }
  return (guess);
}

double fnQy(double guess, double *Qy)
/* Must return difference between target length webthickness and distance from
   (Px,Py) to (Qx,Qy) */
{
  double webthickness, distbetwpts, answer;

  int Debug = 0;

  Qx = InitBurnRadius[0] *
       100; /* Here working in cm, so convert from meters stored in IBR */
  /* Here is the code for calculating Qy from the guess at Px sent in */
  if (Debug)
    fprintf(outfp, "\nInitBurnRadius[0] %0.14f", InitBurnRadius[0]);
  webthickness =
      7.62 - InitBurnRadius[0] *
                 100; /* InitBurnRadius passed in as meters, convert for here */
  Py = (SemiMinAx * 100.) *
       pow((1. - (guess * guess) / ((SemiMajAx * 100. * SemiMajAx * 100.))),
           0.5);
  if (Debug)
    fprintf(outfp, "\nIn fnQy, Py %0.14f", Py);
  *Qy = Py - ((pow(7.62, 2) * webthickness * Py) /
              pow(pow((SemiMinAx * 100), 4) * guess * guess +
                      pow((SemiMajAx * 100.), 4) * Py * Py,
                  0.5));
  /* Then calculate the quantity to be zeroed */
  distbetwpts = DistBetwPts(Qx, *Qy, guess, Py);
  answer = distbetwpts - webthickness;
  if (Debug)
    fprintf(outfp, "\nanswer to zero, %0.14f, Qy %0.14f", answer, *Qy);
  return (answer);
}

void SetUpGeometry(void) {
  int i, SegNo;
  int Debug = 0;

#ifdef STARGRAIN
  /* Partial initialization needed when EVOLVINGSTARSHAPE, else okay. */
  /* THESE initializations must PRECEDE decodechromosome(); remainder must
   * FOLLOW decodechromosome() */
  /*         nozzle segment into zones,  as matches adjacent segment[1]. */
  for (SegNo = 2; SegNo < NCYLSEGS - 1; SegNo++) {
    /* if NOT EVOLVINGSTARSHAPE, then they determine the line which is the star
     * edge. */
    /* But IF EVOLVINGSTARSHAPE, then they are simply the limits for the
     * positions being */
    /* evolved. */
#ifndef EVOLVEPTAXOFSTAR
    PtAx[SegNo] = INNERSTARPOINT;
#endif
    PtAy[SegNo] = 0.0;
    PtBRadius[SegNo] = InitBurnRadius[SegNo];
    // PtBRadius[SegNo] = InitBurnRadius[SegNo]; /* Use this if want PtB at
    // InitBurnRadius, else use smth LESS. */
    PtBx[SegNo] = PtBRadius[SegNo] * cos(2. * PI / NSYMS);
    PtBy[SegNo] = PtBRadius[SegNo] * sin(2. * PI / NSYMS);
    if (Debug) {
      fprintf(outfp,
              "\nIn SetUpGeometry, PtAx[%d] = %0.14f, PtAy[%d] = %0.14f, "
              "PtBx[%d] = %0.14f, PtBy[%d] = %0.14f\n",
              SegNo, PtAx[SegNo], SegNo, PtAy[SegNo], SegNo, PtBx[SegNo], SegNo,
              PtBy[SegNo]);
    }
  }
  /* Now determine the boundary between the two propellant types in pre-cyl
   * parts of star segments */
  /* Placed here so will be correct whether or not PtAx and/or InitBurnRadius
   * are being evolved. */
  for (SegNo = 2; SegNo < NCYLSEGS - 1; SegNo++) {
    StarRayOuterInnerBdary[SegNo] = (PtAx[SegNo] + InitBurnRadius[SegNo]) / 2.;
    if (Debug)
      fprintf(outfp, "\nPtAx[%d] %0.14f, StarRayOuterInnerBdary[%d] = %0.14f",
              SegNo, PtAx[SegNo], SegNo, StarRayOuterInnerBdary[SegNo]);
  }
#endif

  SetupNozzleGeometry();

  /* Initialize interval between recording thrusts */
  ThrustInterval = (int)(0.5 / DeltaTime + 0.5);

  /* Calculate lengths of initial and final burn chords for the central ellipse
   * burn segment */
  /* (this step allows InitialSegmentRadius to be open for optimization.) */
  initialchord = 2. * InitBurnRadius[0] * 100.; /* fitted equation in cm */
  finalchord = 2.3407279 * InitBurnRadius[0] * 100. + 0.774180243;
  FirstInitGlobalsCall = 0;
}

double DistBetwPts(double pt1x, double pt1y, double pt2x, double pt2y) {
  return (sqrt((pt1x - pt2x) * (pt1x - pt2x) + (pt1y - pt2y) * (pt1y - pt2y)));
}

void SetupInsulationGeometry(void) {

  /* Set up Insulation geometry here, must reset for each new individual */
  /* x coordinates will be tracked relative to corner being x = 0. */
  /* Therefore, must subtract IBR from Ray[8] to get proper length of top. */
  /* BURN depth will be tracked in y relative to initial location treated as 0.
   */
  /* OTHER y coord's will be relative to center of rocket */
  /* --i.e., have InitBurnRadius[0] ( == [1]) added in */

  InsulPtx = 0.;
  InsulPty = InitBurnRadius[1];
  InsulHalfAngle =
      (11. / 180. * PI) / 2.; /* Will track insulation burn along this ray */
  InsulBottomBurnDepth = 0.;
  InsulShellPtx = (0.043) * 4.543 / 3.38;
  InsulShellPty = InitBurnRadius[1];
  InsulTopPtAtShellx =
      NozRayFinalBurnRadius[8] * cos(11. / 180. * PI) - InitBurnRadius[1];
  InsulTopPtAtShelly =
      sin(11. / 180. * PI) * NozRayFinalBurnRadius[8] + InitBurnRadius[1];
  InsulPropBdaryPtx = InsulPtx; /* only initially */
  InsulPropBdaryPty = InsulPty; /* only initially */
  InsulTopBurningLength =
      DistBetwPts(InsulPropBdaryPtx, InsulPropBdaryPty, InsulPtx, InsulPty);
}

#ifdef STARGRAIN
void SetupStarBurn(int SegNo) {
  /* MUST BE SURE DOES NOT OVERWRITE things computed by decodechromosome(). */
  double SectorAngle;
  int i;

  int Debug = 0;

  for (SegNo = 2; SegNo < NCYLSEGS - 1; SegNo++) {
    /* CALLED ONLY DURING INITIALIZATION */
    NRAYSinPreCylLayer[SegNo] = NRAYS + 2;
    /* Now set up the geometry for the "star" burn special interior sections. */
    SectorAngle = 2. * PI / NSYMS;
    StarBurningThisSegNow[SegNo] = 1;
    /* For actual star (linear sides), set up points A and B defining one side
     */
    /* of a star's point. B is deep side (depth starts at 0, closest to shell */
    /* and A is x-axis side (ray [0]), extending toward center an evolvable
     * distance */

    for (i = 0; i < NRAYS + 2; i++) {
      RayInLayer[i][SegNo] = -1;
    }
    /* Note: calculated AFTER decodechromosome, so if PtAx being evolved, is
     * already loaded */
    StarTheta[SegNo] = atan(PtBy[SegNo] / (PtBx[SegNo] - PtAx[SegNo]));
    /* Set up for either straight-sided star or evolved star shape */
    /* Now find the angles of the rays. Consider x axis ray[0] */
    /* and last ray as [NRAYS+1] */
    /* Much of work for straightsided star already done in decodechromosome */
#ifndef EVOLVINGSTARSHAPE
    /* Starting depth for straight-sided star */
    RayBurnDepth[0][SegNo] = PtAx[SegNo];
    /* Deep end of star for straight-sided star */
    RayBurnDepth[NRAYS + 1][SegNo] = InitBurnRadius[SegNo];
#endif /* of NOT EVOLVINGSTARSHAPE */
    RayAngle[0][SegNo] = 0.;
    RayAngle[NRAYS + 1][SegNo] = 2. * PI / NSYMS;
    BurnPtx[0][SegNo] = RayBurnDepth[0][SegNo];
    BurnPty[0][SegNo] = 0.;
    BurnPtx[NRAYS + 1][SegNo] = PtBx[SegNo];
    BurnPty[NRAYS + 1][SegNo] = PtBy[SegNo];
    if (Debug) {
      fprintf(
          outfp,
          "\nIn         SetupStarBurn, RayBurnDepth[0][%d] %0.14f, RayAngle "
          "%0.14f, Burnptx[%d][%d] = %0.14f, BurnPty[%d][%d] = %0.14f",
          SegNo, RayBurnDepth[0][SegNo], RayAngle[0][SegNo], 0, SegNo,
          BurnPtx[0][SegNo], 0, SegNo, BurnPty[0][SegNo]);
    }
    for (i = 1; i < NRAYS + 1; i++) {
      RayAngle[i][SegNo] =
          RayAngle[i - 1][SegNo] + (2. * PI / NSYMS) / (NRAYS + 1);
      /* And find where they intersect the initial burn curve to set BurnPtx,
       * BurnPty for each i */
      /* This defines either straight star or evolved star shape */
      DefineInitStarProfile(i, SegNo);
    }
    if (Debug) {
      fprintf(
          outfp,
          "\nIn         SetupStarBurn, RayBurnDepth[%d][%d] %0.14f, RayAngle "
          "%0.14f, Burnptx[%d][%d] = %0.14f, BurnPty[%d][%d] = %0.14f",
          NRAYS + 1, SegNo, RayBurnDepth[NRAYS + 1][SegNo],
          RayAngle[NRAYS + 1][SegNo], NRAYS + 1, SegNo,
          BurnPtx[NRAYS + 1][SegNo], NRAYS + 1, SegNo,
          BurnPty[NRAYS + 1][SegNo]);
    }
    /* For volume, mass calculations, determine where each ray intersects layer
     * 0 of cylinder */
    for (i = 0; i < NRAYS + 2; i++) {
      RayCylIntersectx[i][SegNo] =
          InitBurnRadius[SegNo] * cos(RayAngle[i][SegNo]);
      RayCylIntersecty[i][SegNo] =
          InitBurnRadius[SegNo] * sin(RayAngle[i][SegNo]);
    }
    /* Now determine the normal vector at each point (direction of burn == of
     * advance of front */
    /* Now calculate on each radial line what the burn depth is. */
    CalcRayBurnDirnsAndLengths(SegNo);
  }
}

void CalcRayBurnDirnsAndLengths(int SegNo) {
  /* ONLY FOR STAR RAY CYL SEGMENTS!!! */
  /* ONLY called during burn, NOT during initialization */

  int i;
  double slope[NRAYS + 3]; /* local temp; slope of line segment between ray i
                              and ray i+1, */
  double denom;
  /* except uses normals as exact slopes at both ends. */
  int Debug = 0;

  if (Debug) {
    for (i = 0; i < NRAYS + 2; i++) {
      fprintf(outfp,
              "\nComing INTO CalcRayBurnDirnsAndLengths, RayBurnDepth[%d][%d] "
              "%0.14f",
              i, SegNo, RayBurnDepth[i][SegNo]);
      fprintf(outfp, "\nBurnPtx[%d][%d] %0.14f, BurnPty[%d][%d] %0.14f", i,
              SegNo, BurnPtx[i][SegNo], i, SegNo, BurnPty[i][SegNo]);
    }
  }
  /* First calculate slopes of all between-ray lines, all NRAYS + 3 of them
   * (counting slopes for end points) */
  /* But will end up with only NRAYS+2 RayBurnDirs. */
  /* If not evolving star shape, slope of initial line is By/(Bx-Ax). */
  /* If evolving, then each ray will have own initial burn direction */
  denom = (BurnPtx[1][SegNo] - BurnPtx[0][SegNo]);
  if (denom != 0.) {
    slope[0] = (BurnPty[1][SegNo] - BurnPty[0][SegNo]) / denom;
  } else {
    slope[0] = 1000.; /* arbitrary, but near enough to infinite to work fine. */
  }
  if (BurnPtx[1][SegNo] - BurnPtx[0][SegNo] > 0.) {
    RayBurnDir[0][SegNo] = atan(slope[0]) - PI / 2.;
  } else {
    RayBurnDir[0][SegNo] = atan(slope[0]) + PI / 2.;
  }
  if (Debug) {
    fprintf(outfp, "\nRayBurnDir[%d][%d] = %0.14f", 0, SegNo,
            RayBurnDir[0][SegNo]);
  }
  denom = (BurnPtx[NRAYS + 1][SegNo] - BurnPtx[NRAYS][SegNo]);
  slope[NRAYS + 1] =
      (BurnPty[NRAYS + 1][SegNo] - BurnPty[NRAYS][SegNo]) / denom;
  if (BurnPtx[NRAYS + 1][SegNo] - BurnPtx[NRAYS][SegNo] > 0.) {
    RayBurnDir[NRAYS + 1][SegNo] = atan(slope[NRAYS + 1]) - PI / 2.;
  } else {
    RayBurnDir[NRAYS + 1][SegNo] = atan(slope[NRAYS + 1]) + PI / 2.;
  }
  for (i = 1; i < NRAYS + 1; i++) {
    denom = (BurnPtx[i][SegNo] - BurnPtx[i - 1][SegNo]);
    slope[i] = (BurnPty[i][SegNo] - BurnPty[i - 1][SegNo]) / denom;
    if (BurnPtx[i][SegNo] - BurnPtx[i - 1][SegNo] > 0.) {
      RayBurnDir[i][SegNo] = atan(slope[i]) - PI / 2.;
    } else {
      RayBurnDir[i][SegNo] = atan(slope[i]) + PI / 2.;
    }
    if (Debug) {
      fprintf(outfp, "\natan(slope[%d]) = %0.14f", i, atan(slope[i]));
      fprintf(outfp, "\nRayBurnDir[%d][%d] = %0.14f", i, SegNo,
              RayBurnDir[i][SegNo]);
    }
  }
  if (Debug) {
    fprintf(outfp, "\nRayBurnDir[%d][%d] = %0.14f", NRAYS + 1, SegNo,
            RayBurnDir[NRAYS + 1][SegNo]);
  }
  /* RayBurnDirs are angles from x axis in radians, NOT slopes. */
  /* first and last rays don't average slopes, will use half lengths */
  for (i = 1; i < NRAYS + 2; i++) {
    if (Debug) {
      fprintf(outfp,
              "\nIn CalcRayBurnDir... BurnPtx[%d][%d] %0.14f, BurnPty[%d][%d] "
              "%0.14f, BurnPtx[%d][%d] %0.14f, BurnPty[%d][%d] %0.14f",
              i, SegNo, BurnPtx[i][SegNo], i, SegNo, BurnPty[i][SegNo], i - 1,
              SegNo, BurnPtx[i - 1][SegNo], i - 1, SegNo,
              BurnPty[i - 1][SegNo]);
    }
    RayiToRayiminus1[i][SegNo] =
        DistBetwPts(BurnPtx[i][SegNo], BurnPty[i][SegNo], BurnPtx[i - 1][SegNo],
                    BurnPty[i - 1][SegNo]);
  }
  /* NOTE:  special half-length calculations for first and last rays are handled
   * in CalcStarRayIntegralContribs(), not here */
  /* This code calculates all that is needed to support that. */
  if (Debug) {
    for (i = 0; i < NRAYS + 2; i++) {
      fprintf(outfp,
              "\nLeaving CalcRayBurnDirnsAndLengths, Ray[%d][%d] has slope "
              "%0.14f, angle %0.14f, distbetwpts %0.14f",
              i, SegNo, slope[i], atan(slope[i]), RayiToRayiminus1[i][SegNo]);
    }
  }
}

void DefineInitStarProfile(int RayNo, int SegNo) {
  /* Calculates initial star burn point (x,y) for one ray in one segment in a
   * cylindrical segment of rocket.  */
  /* NOTE:  This routine is easily extendable to one that will allow a broad
   * choice of initial burn */
  /* profiles, just by entering the intersections of the profile with the rays
   * perpendicular to the */
  /* rocket centerline in a 1/Nth (n-way symmetry or anti-symmetry) section of
   * the circular cross- */
  /* section. Must be defined so as to avoid slivers during the burn. IF the
   * "special layers" are */
  /* sufficient to allow circularization of the burn once it is internal to the
   * cylinder, then it can */
  /* also be designed so as to burn out simultaneously, avoiding cusps or
   * unburned slivers at the */
  /* rocket shell. Define new shape by replacing y=mx+b with the desired shape,
   * calculating the */
  /* intersection with ray formula. */

  double ABslopem, ABintb;
  int Debug = 0;

#ifndef EVOLVINGSTARSHAPE
  /* Calculate intersection point between a radial line and the initial line
   * defining the star edge */
  /* AB in y=mx+b form is: */
  ABslopem = PtBy[SegNo] / (PtBx[SegNo] - PtAx[SegNo]);
  ABintb = -ABslopem * PtAx[SegNo];
  /* Then intersection with ray RayNo is: */
  BurnPty[RayNo][SegNo] =
      ABintb / (1. - ABslopem / tan(RayAngle[RayNo][SegNo]));
  BurnPtx[RayNo][SegNo] = BurnPty[RayNo][SegNo] / tan(RayAngle[RayNo][SegNo]);
  RayBurnDepth[RayNo][SegNo] =
      DistBetwPts(BurnPtx[RayNo][SegNo], BurnPty[RayNo][SegNo], 0., 0.);
  if (Debug) {
    fprintf(outfp, "\nFrom DefineInitStarProfile, RayBurnDepth[%d][%d] %0.14f",
            RayNo, SegNo, RayBurnDepth[RayNo][SegNo]);
    fprintf(outfp, "\nFrom DefineInitStarProfile, BurnPtx[%d][%d] %0.14f",
            RayNo, SegNo, BurnPtx[RayNo][SegNo]);
    fprintf(outfp, "\nFrom DefineInitStarProfile, BurnPty[%d][%d] %0.14f",
            RayNo, SegNo, BurnPty[RayNo][SegNo]);
  }
#else
  /* The profile is being evolved, radius */
  /* was already decoded in decodechromosome() */
  BurnPtx[RayNo][SegNo] =
      cos(RayAngle[RayNo][SegNo]) * RayBurnDepth[RayNo][SegNo];
  BurnPty[RayNo][SegNo] =
      sin(RayAngle[RayNo][SegNo]) * RayBurnDepth[RayNo][SegNo];
  if (Debug) {
    fprintf(outfp,
            "\nIn DefineInitStarProfile, RayBurnDepth[%d][%d] %0.14f, RayAngle "
            "%0.14f, Burnptx[%d][%d] = %0.14f, BurnPty[%d][%d] = %0.14f",
            RayNo, SegNo, RayBurnDepth[RayNo][SegNo], RayAngle[RayNo][SegNo],
            RayNo, SegNo, BurnPtx[RayNo][SegNo], RayNo, SegNo,
            BurnPty[RayNo][SegNo]);
  }
#endif /* of NOT EVOLVINGSTARSHAPE */
}
#endif /* of STARGRAIN */

void InitializeRemainingGlobals(void) {
  int i, SegNo;
  double CylLengthForCylSegs;

  double Range, MidPt, Stepsize, FirstStep, Denom, SumSegLengths;

  int Debug = 0;

  FirstAnimationCall = 1;
  if (FirstInitGlobalsCall) {
    SetUpGeometry();
    FirstInitGlobalsCall = 0;
  }
  Alphaminus1 = Alpha - 1.;
  InvPrefAlpha = 1. / pow(RefPressure, Alpha);
  TotalPressure = RefPressure; // initial guess
  NozBurningRays = 0;
  SetupInsulationGeometry();
  /* CHANGE seg length taking into account seg 0, the central dome segment. */
  /* Other CYL segs' total is shorter by amount fabs(Qy). */
  /* Correct cyl length assignable to truly cyl segs for amt used in nozzle and
   * dome */
  /* We do NOT set a length for seg 0, as it will be calculated dynamically. */
  CylLengthForCylSegs = 0.5583 - fabs(Qy) - 0.02095;
  /* For now, setting segs at ends smaller, they determine other segs' lengths.
   */

  SegmentLength[NCYLSEGS - 1] = 0.05;
  SegmentLength[1] = 0.05;
  CylLengthForCylSegs -= 0.1;
  Range = (NCYLSEGS - 3 - 1) / 2.;
  MidPt = (double)(NCYLSEGS - 3) / 2.;
  Stepsize = 0.5;
  FirstStep = MidPt - Stepsize * Range;
  /* Sum should be half of Fibonacci number */
  SegmentLength[2] = FirstStep;
  Denom = FirstStep;
  SumSegLengths = 0.;
  for (SegNo = 3; SegNo < NCYLSEGS - 1;
       SegNo++) { /* divide the remaining length UNEQUALLY among other cyl segs
                   */
    SegmentLength[SegNo] = SegmentLength[2] + Stepsize * (SegNo - 2);
    Denom += SegmentLength[SegNo];
  }
  for (SegNo = 2; SegNo < NCYLSEGS - 1; SegNo++) {
    SegmentLength[SegNo] = (SegmentLength[SegNo] / Denom) * CylLengthForCylSegs;
    SumSegLengths += SegmentLength[SegNo];
    if (Debug)
      fprintf(outfp,
              "\nDenom %0.14f, SegmentLength[%d] %0.14f, SumSegLengths %0.14f, "
              "CylLengthForCylSegs %0.14f",
              Denom, SegNo, SegmentLength[SegNo], SumSegLengths,
              CylLengthForCylSegs);
#ifdef STARGRAIN
    for (i = 0; i < NRAYS + 2; i++) {
      LatestCatchupDepthBurned[i][SegNo] = 0.;
      CatchingUpRay[i][SegNo] = 0;
    }
#endif
  }
  for (i = 0; i < NumTimePointsCoarse; i++) {
    ThrustProfile[i] = 0.;
    mPressureProfile[i] = 0.;
  }
  MaxNumNewtonItr = 1000;
  RefPressure = 3.447e6; /* Pa */

#ifdef STARGRAIN
  FirstRayCallThisGuy[SegNo] = 1; /* Flag to control initializing print */
  for (SegNo = 2; SegNo < NCYLSEGS - 1; SegNo++) {
    CompensatingNow[SegNo] = 0;
  }
#endif /* of STARGRAIN */
  /* Initializations for each segment burned */
  for (SegNo = 0; SegNo < NTOTALSEGS; SegNo++) {
    /* SegmentRadius is the variable to update as the burn progresses,
     * incrementing each step by RDot */
    SegmentRadius[SegNo] =
        InitBurnRadius[SegNo]; /* Initialize at proper initial radius */
    CurrLayerIndex[SegNo] = 0;
    // TriangleBurn[SegNo] = 0;
#ifdef STARGRAIN
    StarBurningThisSegNow[SegNo] = 0; /* will override where it is true. */
#endif
  }

  NewtonTolerance = 1.e-8; /* Pa */
  Time = 0.;
  Thrust = 0.;
  Burning = 1; /* Initially TRUE */
  Reward = 0.; /* once per critter */

  if (WantToAnimate == neval) {
    if (!(Animatfp = fopen("./AnimationFile", "w"))) {
      printf("AnimationFile %s cannot be opened for write.\n", "AnimationFile");
    }
    WriteGraphicsOutputFile(neval);
  }
  /* End of rocket variable re-initialization */
}

void CalcBurnIntegralContributions(void) {
  double temp;
  int SegNo;

  int Debug = 0;

  IntegralContribs = 0.;
  /* Contributions to (and subtraction from) burn integral for mdotin tallied
   * here */
  for (SegNo = 0; SegNo < NTOTALSEGS; SegNo++) {
#ifdef STARGRAIN
    if (SegNo > 1 && SegNo < NCYLSEGS - 1 &&
        StarBurningThisSegNow[SegNo] == 1) {
      CalcStarBurnIntegrContribs(SegNo);
      IntegralContribs += SumOfProducts[SegNo];
      if (Debug) {
        fprintf(outfp, "\n Seg %d Ray Integral Contributions %e", SegNo,
                SumOfProducts[SegNo]);
      }
    } else if (SegNo == NTOTALSEGS - 1 && NozBurningRays) {
#else
    if (SegNo == NTOTALSEGS - 1 && NozBurningRays) {
#endif
      temp = CalcNozzleRaysTotalIntegralContribs();
      IntegralContribs += temp;
      if (Debug) {
        fprintf(outfp, "\nNozzle Integral Contributions %e, total after %e\n",
                temp, IntegralContribs);
      }
    } else {
      temp = SegmentArea[SegNo] * RefBurnRate[SegNo];
      if (RayDebug && SegNo > 1 && SegNo < 5)
        fprintf(outfp, "\nSegmentArea[%d] %0.14f, RefBurnRate[%d] %0.14f",
                SegNo, SegmentArea[SegNo], SegNo, RefBurnRate[SegNo]);
      IntegralContribs += temp;
      if (Debug) {
        fprintf(outfp, "\nSeg %d Integral Contributions %e, Total after: %e\n",
                SegNo, temp, IntegralContribs);
      }
    }
  }
}

double CalcNozzleRaysTotalIntegralContribs(void) {
  int i;
  double contribstotal = 0.;

  int Debug = 0;

  /* Nozzle segment burning rays, must add their contributions */
  for (i = 0; i < NNOZRAYS; i++) {
    contribstotal += NozRayRefBurnRate[i] * NozRaySegBurnArea[i];
    if (Debug) {
      fprintf(outfp,
              "\nNozRayRefBurnRate[%d] %0.14f, NozRaySegBurnArea[%d] %0.14f", i,
              NozRayRefBurnRate[i], i, NozRaySegBurnArea[i]);
    }
  }
  if (Debug) {
    fprintf(outfp, "\nNozzleRay integral contribs ray total %0.14f",
            contribstotal);
  }
  return contribstotal;
}

int PressureViolation(void) {
  int Debug = 0;

  int ReturnValue = 0;

  if (Debug) {
    fprintf(outfp,
            "\nIn PressureViolation(), DeltaPressure = %0.14f, TotalPressure = "
            "%0.14e",
            DeltaPressure, TotalPressure);
  }
  if (DeltaPressure != DeltaPressure) {
    /* means result was NaN not a number or out of range */
    StopCode = 4; /* means result was "funny" */
    TotalPressure = -555.55555;
    ReturnValue = 1;
  }
  /* TRY relaxing minpressure when near 10 seconds, as TARGET would violate
     minpressure. */
  /* Since they allow burning past 10 seconds, relax the MINPRESSURE constraint
     after that, as their */
  /* targets ALL require that for the first, regressive thrust profile. */
  else if ((Time > 8. && (Time <= 9.5 && (DeltaPressure < 0. &&
                                          (DeltaPressure + TotalPressure <
                                           MINPRESSURE * 0.9)))) ||
           (Time > 9.5 && (DeltaPressure < 0. &&
                           (DeltaPressure + TotalPressure <
                            MINPRESSURE * 0.9 * pow(.95, (irecord - 19))))) ||
           (Time <= 8. && (DeltaPressure < 0. &&
                           (DeltaPressure + TotalPressure < MINPRESSURE)))) {

    /* UnderPressure */
    StopCode = 2;
    TotalPressure = -555.55555;
    ReturnValue = 1;
  } else if (TotalPressure > MaxPressure ||
             (DeltaPressure > 0. &&
              TotalPressure + DeltaPressure > MaxPressure)) {
    if (Debug) {
      fprintf(
          outfp,
          "\nTotalPressure calculated IN NEWTON was out of range, = %0.14f\n",
          TotalPressure);
      /* THIS causes an early stop to this guy's eval... set fitness, returning
       * 1 will end burn */
      fprintf(outfp, "\nStopping burn at time %0.14f, bad pressure %0.14f.\n",
              Time, TotalPressure);
    }
    StopCode = 1;
    TotalPressure = -555.55555;
    ReturnValue = 1;
  }
  if (Debug && ReturnValue == 0) {
    fprintf(outfp,
            "\nneval %ld Time %0.14f Pressure is acceptable, proceed to next "
            "timestep.\n",
            neval, Time);
  }
  return (ReturnValue);
}

int CalcSegRefBurnRates() {
  int SegNo, i;
  int Debug = 0;

  for (SegNo = 0; SegNo < NTOTALSEGS; SegNo++) {
    RefBurnRate[SegNo] = FindPropellantType(SegNo);
    /* Check whether burn has hit shell; if so, finish this rocket; fitness
     * already calc'd. */
    if (RefBurnRate[SegNo] == -222.22222) {
      return (-1);
    }
    if (Debug) {
      if (SegNo == NTOTALSEGS - 1 && NozBurningRays) {
        /* Nozzle segment, rays called something different */
        for (i = 0; i < NNOZRAYS; i++) {
          fprintf(outfp, "\nNozRayRefBurnRate[%d] = %0.14f", i,
                  NozRayRefBurnRate[i]);
        }
      }
    }
    if (RefBurnRate[SegNo] == -666.66666) {
      /* Burning in final two layers of a corner seg. Value of RDotRef has
       * already been set; do not overwrite */
      // BUG fixed 4/12/19 EDG
      RefBurnRate[SegNo] = RDotRef[CurrLayerIndex[SegNo]][SegNo];

      if (Debug) {
        fprintf(outfp, "\nSpecial corner RDotRef[%d] = %0.14f", SegNo,
                RDotRef[CurrLayerIndex[SegNo]][SegNo]);
      }
    }
#ifdef STARGRAIN
    /* If in a currently STAR burning seg, a -777.77777 code is returned,
     * meaning that the RefBurnRate value */
    /* is junk for this segment, and the per-ray values of interest are stored
     * in RayRDotRef[rayno][segno]. */
    if (RefBurnRate[SegNo] == -777.77777) {
      if (Debug) {
        for (i = 0; i < NRAYS + 2; i++) {
          fprintf(outfp, "\nRayRDotRef[%d][%d] = %0.14f", i, SegNo,
                  RayRDotRef[i][SegNo]);
        }
      }
    }
#endif
    if (Debug)
      fprintf(outfp, "\nPropellantType for seg %d is %0.14f\n", SegNo,
              RefBurnRate[SegNo]);
    if (RefBurnRate[SegNo] > 1.525) {
      fprintf(outfp,
              "\nneval %ld, Bad RefBurnRate[%d] =%0.14f, time =%0.14f, "
              "SegmentRadius[%d] %0.14f, CurrLayerIndex[%d] %d\n",
              neval, SegNo, RefBurnRate[SegNo], Time, SegNo,
              SegmentRadius[SegNo], SegNo, CurrLayerIndex[SegNo]);
      printf("\nneval %ld, Bad RefBurnRate[%d] =%0.14f, time =%0.14f, "
             "SegmentRadius[%d] %0.14f, CurrLayerIndex[%d] %d\n",
             neval, SegNo, RefBurnRate[SegNo], Time, SegNo,
             SegmentRadius[SegNo], SegNo, CurrLayerIndex[SegNo]);
      exit(992);
    }
  }
  return 1;
}

void CalcSegBurnAreas() {
  double ttemp, avgperimeter;
  double r, R, thetainit, thetafinal;
  int i, SegNo;

  int Debug = 0;

  /* NOTE:  Segs now burning star rays have their burn areas calculated in
   * another function */

  /* CYLINDRICAL (non-dome) SEGMENT BURN AREA CALCULATIONS */
  for (SegNo = 1; SegNo < NCYLSEGS; SegNo++) {
#ifdef STARGRAIN
    if (SegNo == 1 || SegNo == NCYLSEGS - 1 || !StarBurningThisSegNow[SegNo]) {
#endif
      /* Now, do this for each segment on cylinder part, NOT dome segment 0 and
       */
      /* not a star seg still burning "specially". */
      /* We are getting SegmentArea[SegNo] */
      /* CALCULATE current area to multiply local burn value by */
      if (SegNo > 0 && SegNo < NCYLSEGS) {
        SegmentArea[SegNo] =
            SegmentLength[SegNo] * SegmentRadius[SegNo] * 2. * PI;
      }
      if (Debug && SegNo > 0) {
        fprintf(outfp,
                "\nCyl seg SegmentArea[%d] = %0.14e, SegmentLength[%d] = %e, "
                "SegmentRadius[%d] = %0.14e\n",
                SegNo, SegmentArea[SegNo], SegNo, SegmentLength[SegNo], SegNo,
                SegmentRadius[SegNo]);
      }
#ifdef STARGRAIN
    }
#endif
  }
  /* CENTRAL ELLIPTICAL DOME BURN AREA CALCULATIONS */
  /* Now do the dome central elliptical burn area calculation. Does a domed
   * circle depending on burn depth. */
  /* Initial central elliptical segment chordal width is 2*Qx ==
   * 2*InitBurnRadius. */
  /* Final central elliptical segment chordal width is 2*Px. Now solve for the
   * chordal length as a */
  /* function of burn depth. Then convert to arc length to calculate burn area.
   */
  /* From the spreadsheet, in general, calculated chord length as fn of burn
   * depth. */
  /* Now can calculate chord length at PRESENT depth of burn. initialchord and
   * finalchord were */
  /* calculated in cm above. */
  depth = SegmentRadius[0] - InitBurnRadius[0];
  chordlength =
      (((finalchord - initialchord) / 100.) / InitBurnRadius[0]) * depth +
      2 * InitBurnRadius[0];
  /* Now, TEMPORARILY!!!, convert chord length to length of arc along ellipse,
   * using WRONG conversion */
  SegmentLength[0] =
      chordlength * 1.111; /* WILL FIX SOON with a bisection call */
  ttemp = (SegmentLength[0] / 2.) *
          (SegmentLength[0] / 2.); /* i.e., "modified" r^2 */
  ttemp = ttemp * PI;
  SegmentArea[0] = ttemp;
  if (Debug) {
    fprintf(outfp,
            "\nCentral dome initialchord %0.14f cm, finalchord %0.14f cm, "
            "depth %0.14f, chordlength %0.14f, SegmentLength[0] %0.14f, "
            "SegmentArea[0] %0.14f\n",
            initialchord, finalchord, depth, chordlength, SegmentLength[0],
            SegmentArea[0]);
  }

  /* CORNER SEGMENT BURN AREA CALCULATIONS */
  /* Now can calculate the segment widths for Corner segments as a function of
   */
  /* burn depth using the half-angles between the boundaries. */
  for (SegNo = 0; SegNo < NCORNERSEGS; SegNo++) {
    rayangle = (BdaryAngle[SegNo] + BdaryAngle[SegNo + 1]) / 2.;
    /* But area calculated instead as a portion of area of the torus it belongs
     * to. */
    if (Debug) {
      fprintf(outfp, "\n*** cornerseg %d, rayangle %0.14f", SegNo, rayangle);
    }
    r = SegmentRadius[SegNo + NCYLSEGS] - InitBurnRadius[SegNo + NCYLSEGS];
    R = InitBurnRadius[SegNo + NCYLSEGS];
    thetainit = BdaryAngle[SegNo];
    thetafinal = BdaryAngle[SegNo + 1];
    SegmentArea[SegNo + NCYLSEGS] = 2. * PI * r *
                                    ((R * thetafinal + r * sin(thetafinal)) -
                                     (R * thetainit + r * sin(thetainit)));

    if (Debug) {
      fprintf(outfp,
              "\nCorner seg r %0.14f, R %0.14f, thetainit %0.14f, thetafinal "
              "%0.14f, SegmentArea[%d] = %e\n",
              r, R, thetainit, thetafinal, SegNo + NCYLSEGS,
              SegmentArea[SegNo + NCYLSEGS]);
    }
  }
  /* Nozzle Segment Burn Area Calculations Here */
  if (NozBurningRays) {
    NozRayBurnPty[0] = NozRayRadius[0]; /* On vertical ray */
    for (i = 1; i < NNOZRAYS; i++) {
      /* Calculate y coordinate of nozzle rays for use in calculating perimeters
       */
      NozRayBurnPty[i] = sin(NozRayAngle[i]) * NozRayRadius[i];
      if (Debug) {
        fprintf(outfp, "\nNozRayBurnPty[%d] = %0.14f", i, NozRayBurnPty[i]);
      }
    }
    for (i = 0; i < NNOZRAYS - 1; i++) {
      /* Approximating burn area when burning rays by calculating the average of
       * the perimeters of each pair */
      /* of adjacent rays, then multiplying that by the arc length between the
       * adjacent rays to get surf area. */
      avgperimeter =
          (NozRayBurnPty[i] * 2. * PI + NozRayBurnPty[i + 1] * 2. * PI) / 2.;
      /* Now area is just average perimeter times arc length between rays. */
      NozRaySegBurnArea[i] = ((NozRayAngle[i] - NozRayAngle[i + 1]) *
                              (NozRayRadius[i] + NozRayRadius[i + 1]) / 2.) *
                             avgperimeter;
      if (Debug) {
        fprintf(outfp, "\nAvg perimeter %d %0.14f, NozRaySegBurnArea %d %0.14f",
                i, avgperimeter, i, NozRaySegBurnArea[i]);
        fprintf(outfp, "\nNozRayRadius[%d] %0.14f", i, NozRayRadius[i]);
      }
    }
  } else {
    /* Burning nozzle seg, but NOT yet RAYS */
    /* Approximating burn area as for corner segs: a portion of area of the
     * torus it is on. */
    r = SegmentRadius[NTOTALSEGS - 1] - InitBurnRadius[NTOTALSEGS - 1];
    R = InitBurnRadius[NTOTALSEGS - 1];
    thetainit = NozRayAngle[8];
    thetafinal = NozRayAngle[0];
    SegmentArea[NTOTALSEGS - 1] = 2. * PI * r *
                                  ((R * thetafinal + r * sin(thetafinal)) -
                                   (R * thetainit + r * sin(thetainit)));
    if (Debug) {
      fprintf(outfp,
              "\nNOZZLE seg r %0.14f, R %0.14f, thetainit %0.14f, thetafinal "
              "%0.14f, SegmentArea[%d] = %e\n",
              r, R, thetainit, thetafinal, NTOTALSEGS - 1,
              SegmentArea[NTOTALSEGS - 1]);
    }
  }

  /* INSULATION BURN AREA CALCULATED HERE */
  /* Insul bottom burn length is current tip point to current back point. */
  /* Now using InsulBottomBurnDepth as 0. at beginning, not InitBurnRadius */
  InsulBottomBurnArea = (InsulShellPtx - InsulPtx) * 2. *
                        (InsulBottomBurnDepth + InitBurnRadius[1]) * PI;
  /* Calculating conic section volume as perimeter at midpoint times length of
   * side */
  InsulPropBdaryPerimeter = (InsulPropBdaryPty + InsulPty) / 2. * 2. * PI;
  InsulPropBdaryBurnArea = InsulPropBdaryPerimeter * InsulTopBurningLength;
  InsulBurnArea = InsulBottomBurnArea + InsulPropBdaryBurnArea;
  if (Debug) {
    fprintf(outfp,
            "\nInCalcSegBurnAreas, InsulBottomBurnDepth %0.14e, "
            "InsulBottomBurnArea %0.14e, InsulPropBdaryBurnArea %0.14e",
            InsulBottomBurnDepth, InsulBottomBurnArea, InsulPropBdaryBurnArea);
    fprintf(outfp,
            "\nIn CalcSegBurnAreas, InsulPtx %0.14e, InsulPty %0.14e, "
            "InsulShellPtx %0.14e, InsulShellPty %0.14e, "
            "InsulPropBdaryPerimeter %0.14e",
            InsulPtx, InsulPty, InsulShellPtx, InsulShellPty,
            InsulPropBdaryPerimeter);
    fprintf(outfp,
            "\nInsulBottomBurnArea %0.14e, InsulPropBdaryBurnArea %0.14e, "
            "InsulBurnArea %0.14e",
            InsulBottomBurnArea, InsulPropBdaryBurnArea, InsulBurnArea);
  }
}

int ReadGenome(int *AlreadyJammedFlag) {
  /* First will deal with already jammed data, so can just read into the working
   * array */
  /* Already jammed means that segment 0's zone 0 and zone 0 of all corner
   * segments */
  /* are identical to that of seg [NCYLSEGS-1] and the last two layers of the
   * corner segs */
  /* will be disregarded (both thicknesses and propellant types), although they
   * will be read in. */

  int SegNo, Layer, PhenoPtr;
  int Debug = 0;

  if (!*AlreadyJammedFlag) {
    printf("\n*** Sorry, that kind of input is not yet implemented... you need "
           "to put LEGAL, properly constrained data into the input.\n");
    exit(998);
  }
  if (!(PropellantToEvalfp = fopen("./PropellantToEval", "r"))) {
    printf("Input file %s cannot be read.\n", "PropellantToEval");
  }
  if (!(ThicknessesToEvalfp = fopen("./ThicknessesToEval", "r"))) {
    printf("Input file %s cannot be read.\n", "ThicknessesToEval");
  }
  /* Expect to read NTOTALSEGS * NLAYERS * 2 worth of data:  first, for each
   * segment, the NLAYERS propellant types of each layer, */
  /* and the NLAYERS layer thicknesses for that layer; these are read segment by
   * segment. They will all be read into a single */
  /* 1-D array. */

  /* With Abhiroop, modifying so reads propellant type INTS, but actual real
   * thicknesses of layers. */
  PhenoPtr = 0;
  for (SegNo = 0; SegNo < NTOTALSEGS; SegNo++) {
    for (Layer = 0; Layer < NLAYERS; Layer++) {
      fscanf(PropellantToEvalfp, "%d", &(PhenoValues[PhenoPtr++]));
    }
  }
  fclose(PropellantToEvalfp);
  for (SegNo = 0; SegNo < NTOTALSEGS; SegNo++) {
    for (Layer = 0; Layer < NLAYERS; Layer++) {
      fscanf(ThicknessesToEvalfp, "%lf", &(LayerThicknesses[Layer][SegNo]));
    }
  }
  fclose(ThicknessesToEvalfp);
  if (Debug) {
    for (SegNo = 0; SegNo < NTOTALSEGS; SegNo++) {
      for (Layer = 0; Layer < NLAYERS; Layer++) {
        printf("\nLayerThicknesses[%d][%d] %f", Layer, SegNo,
               LayerThicknesses[Layer][SegNo]);
      }
    }
  }
  if (remove("./PropellantToEval") == 0) {
    // printf("Deleted successfully");
  } else {
    printf("Unable to delete PropellantToEval file");
  }
  if (remove("./ThicknessesToEval") == 0) {
    // printf("Deleted successfully");
  } else {
    printf("Unable to delete ThicknessesToEval file");
  }
  return (1);
}

int DistributeHEEDSInputs(void) {
  int SegNo, i;
  int RayPtr;
  int PhenoPtr;
  int j;
  int Depth, DepthStep, PrevDepth;

  int Debug = 0;

  PhenoPtr = 0;
  /* REVERSING ORDER OF LOOPS */
  for (i = 0; i < NLAYERS; i++) {
    for (SegNo = 0; SegNo < NTOTALSEGS; SegNo++) {
      RDotRef[i][SegNo] = PropellChoices[PropellantToEval[i + NLAYERS * SegNo]];
      PropellType[i][SegNo] = PropellantToEval[i + NLAYERS * SegNo];
    }
  }
  /* Now decode the Star-related inputs, if any */
  PhenoPtr = NLAYERS * 2 * NTOTALSEGS;
#ifdef STARGRAIN
  for (SegNo = 2; SegNo < NCYLSEGS - 1; SegNo++) {
    for (i = 0; i < NRAYS + 1; i++) {
      PhenoValues[PhenoPtr++] =
          RayBurnDepthInput[(SegNo - 2) * (NRAYS + 1) + i];
    }
    PhenoValues[PhenoPtr++] = 0;
    /* Read Propellant type for  star portion, map to correct range */
    /* Now putting this into two locations, outer and inner */
    PreCylStarRDotRef[0][SegNo] =
        PropellChoices[PreCylStarRDotRefInput[SegNo - 2] / 6 + 5];
    PreCylStarRDotRef[1][SegNo] =
        PropellChoices[PreCylStarRDotRefInput[SegNo - 2] % 6 + 5];
    RDotRefCatchup[SegNo] =
        PropellChoices[RDotRefCatchupInput[SegNo - 2]]; /* SHOULD be a
                                                           fast-burning one */
    PhenoPtr += 2;
    if (Debug) {
      fprintf(outfp, "\nPreCylStarRDotRef[0][%d] %0.14f", SegNo,
              PreCylStarRDotRef[0][SegNo]);
      fprintf(outfp, "\nPreCylStarRDotRef[1][%d] %0.14f", SegNo,
              PreCylStarRDotRef[1][SegNo]);
      fprintf(outfp, "\nRDotRefCatchup[%d] %0.14f", SegNo,
              RDotRefCatchup[SegNo]);
    }
  }

#ifdef EVOLVINGSTARSHAPE

  RayPtr = 0;
  for (SegNo = 2; SegNo < NCYLSEGS - 1; SegNo++) {
    /* Start from *RayBurnDepthInput passed in, NRAYS+1 for this segment */
    /* Rays start at NLAYERS * 2 * NTOTALSEGS, goes up and next seg index starts
     * at that + NRAYS+4 */
    PhenoValues[NLAYERS * 2 * NTOTALSEGS + NRAYS + 1 +
                (SegNo - 2) * (NRAYS + 4)] =
        0; /* at Pt B, farthest from center of rocket */
    DepthStep = 1;
    Depth = 0;
    PrevDepth = Depth;
    if (Debug || RayDebug) {
      fprintf(
          outfp,
          "\nRay depths phenotype graph, from Pt B (ray[NRAYS+1] to [0])\n");
      fprintf(outfp, "                |\n");
    }
    for (i = NRAYS; i >= 0; i--) {
      if (RayBurnDepthInput[RayPtr] == 0 && DepthStep > 1) {
        DepthStep -= 1;
        if (Depth > DepthStep + 1) {
          /* If Depth would be taken negative, or DepthStep is 0, this acts as a
           * NOP */
          Depth -= DepthStep;
        }
      } else if (RayBurnDepthInput[RayPtr] == 1) {
      } else if (RayBurnDepthInput[RayPtr] == 2) {
        DepthStep += 1;
        Depth += DepthStep;
      } else if (RayBurnDepthInput[RayPtr] == 3) {
        if (Depth > PrevDepth) {
          DepthStep += 2;
          Depth += DepthStep;
        } else {
          /* Don't allow abrupt and dramatic reversal, narrow burn channel */
          DepthStep += 2;
          Depth += 1;
        }
      }
      if (Depth > NRAYDEPTHS - 1)
        Depth = NRAYDEPTHS - 1;

      // Added by Abhiroop Ghosh. Nov 2, 2019.
      // Ray depths can be directly supplied by the user by setting rayDepthFlag
      // = 1
      if (rayDepthFlag == 0) {
        PhenoValues[NLAYERS * 2 * NTOTALSEGS + i + (SegNo - 2) * (NRAYS + 4)] =
            Depth;
        RayDepthHEEDS[RayPtr] = Depth;
        PrevDepth = Depth;
      } else {
        PhenoValues[NLAYERS * 2 * NTOTALSEGS + i + (SegNo - 2) * (NRAYS + 4)] =
            RayBurnDepthInput[RayPtr];
        RayDepthHEEDS[RayPtr] = RayBurnDepthInput[RayPtr];
        PrevDepth = RayBurnDepthInput[RayPtr];
      }
      RayPtr++;
      if (Debug || RayDebug) {
        for (j = 0; j < NRAYDEPTHS - Depth; j++) {
          fprintf(outfp, " ");
        }
        for (j = 0; j < Depth; j++) {
          fprintf(outfp, "*");
        }
        fprintf(outfp, "|\n");
      }
    }
    if (Debug) {
      for (i = 0; i < NRAYS + 1; i++) {
        fprintf(outfp, "\n  PhenoValues[%d] %d",
                NLAYERS * 2 * NTOTALSEGS + (NRAYS + 4) * (SegNo - 2) + i,
                PhenoValues[NLAYERS * 2 * NTOTALSEGS +
                            (NRAYS + 4) * (SegNo - 2) + i]);
      }
    }
  }
#endif /* of EVOLVINGSTARSHAPE */
#endif /* of STARGRAIN */
  return (1);
}

/***********************************************************
 * DECODING OF CHROMOSOME
 * EDG Revising 10/16/18 so can do any number of segments
 * and 11/1 so can read corner segments
 ***********************************************************/
void decodechromosome(void) {
  int i, j, SegNo, intsum, ChrPtr, PhePtr, PheBase, PheCopyFrom;
  int PhePtrStartofRayDepths;
  int ChromPtrStartofRayDepths;
  int Depth, DepthStep, PrevDepth;
  double FirstZoneThickness;
  int DistFromCenterCode;
  double Factor;
  int Expon;
  int OuterLayer, InnerLayer;

  /*DEBUG HERE */
  int Debug = 0;

  /* The InitBurnRadius is known, so can do the initializations needed before
   * the layer thickness */
  /* decoding can be done.  Call for those initializations HERE. */
  InitializationsNeededForDecoding();
  /* Now translate the CODES in the Phenotype array into their actual numerical
   * values */
#ifdef EVOLVINGSTARSHAPE
  /* Star shape not available for first, second, or last Cylinder Segs */
  /* (dome, seg next to nozzle, seg next to dome) */
  /* Ray[0] is along the x axis, angle 0. radians. That depth is evolved*/
  /* Depth of Ray[NRAYS+1] is not evolved, but set to determine deepest */
  /* initial depth (closest to shell). When InitBurnRadius is later */
  /* evolved, this depth will agree with it. */
  for (SegNo = 2; SegNo < NCYLSEGS - 1; SegNo++) {
    PtBRadius[SegNo] = InitBurnRadius[SegNo];
    if (Debug)
      fprintf(outfp, "\n*** *** *** decoding init ray profile***");
    for (i = 0; i < NRAYDEPTHS; i++) {
      RayInitDepthChoices[i][SegNo] =
          PtBRadius[SegNo] -
          i * (PtBRadius[SegNo] - PtAx[SegNo]) / (NRAYDEPTHS - 1);
      if (Debug) {
        fprintf(outfp, "\n NRAYDEPTHS %d", NRAYDEPTHS);
        fprintf(outfp,
                "\nIn decodechromosome, PtBRadius[%d] = %0.14f, PtAx[%d] "
                "%0.14f, RayInitDepthChoices[%d][%d] = %0.14f",
                SegNo, PtBRadius[SegNo], SegNo, PtAx[SegNo], i, SegNo,
                RayInitDepthChoices[i][SegNo]);
      }
    }
    PhePtr = 2 * NLAYERS * NTOTALSEGS + (SegNo - 2) * (NRAYS + 4);
    for (i = 0; i < NRAYS + 1; i++) {
      RayBurnDepth[i][SegNo] =
          RayInitDepthChoices[PhenoValues[PhePtr++]][SegNo];
      if (Debug) {
        fprintf(
            outfp,
            "\nIn Decodechromosome, Segno %d, Ray depth PhenoValues[%d] = %d",
            SegNo, i, PhenoValues[PhePtr - 1]);
        fprintf(outfp,
                "\nIn decodechromosome, Initial RayBurnDepth[%d][%d] = %0.14e",
                i, SegNo, RayBurnDepth[i][SegNo]);
      }
    }
    RayBurnDepth[NRAYS + 1][SegNo] = PtBRadius[SegNo];
    if (Debug) {
      fprintf(outfp,
              "\nIn Decodechromosome, Initial RayBurnDepth[%d][%d] = %0.14f",
              (NRAYS + 1), SegNo, RayBurnDepth[NRAYS + 1][SegNo]);
    }
  }
#endif /* of EVOLVINGSTARSHAPE */
  PhePtr = 2 * NLAYERS * NTOTALSEGS;
#ifdef STARGRAIN
  PhePtr += 2 * (NCYLSEGS - 3);
#ifdef EVOLVINGSTARSHAPE
  PhePtr += (NCYLSEGS - 3) * (NRAYS + 2);
#endif /* of EVOLVINGSTARSHAPE */
#endif /* of STARGRAIN */
  /* At this point, all segs on the chromosome are decoded into two arrays for
   * each segment, */
  /* specifying the propellant design:  RDotRef[i][SegNo], which specifies the
   * reference burn */
  /* rates in each layer, and LayerStartRadius[i][SegNo], which tells at what
   * radius this layer */
  /* starts. Final layer automatically ends at the inside shell of the rocket.
   * Also decoded are */
  /* the initial RayBurnDepths if EVOLVINGSTARSHAPE is set. */
  if (Debug) {
    for (SegNo = 0; SegNo < NTOTALSEGS; SegNo++) {
      fprintf(outfp, "\n Segment  RDotRef     LayerStartRadius\n");
      for (i = 0; i < NLAYERS + 1; i++) {
        fprintf(outfp, "%d,          %e            %e\n", SegNo,
                RDotRef[i][SegNo], LayerStartRadius[i][SegNo]);
      }
    }
  }
}

/**********************************************************************************
 * Computes fitness measure for this individual
 * Returns value for for critter->init_fitness
 ***********************************************************************************/
double calc_fitness(void) {
  int i, k, lastitoprint;
  int WantInnovizationOutput = 0;
  int WantThisGuyOutputFile = 1;

  int PhenoPtrStartofRayDepths, PhePtr;
  double RemainingDistanceToBurn;
  double SumSqRemDist;
  double StdDevRemDist;
  double Variance;
#ifdef STARGRAIN
  double RayRemainder[NRAYS + 2][NTOTALSEGS];
  ;
#endif
  double Difference;
  double Term1 = 0.;
  double Term2 = 0.;
  int SegNo;

  int Debug = 0;

  WantWriteRestartData = 0;

  /* thrusts and times are in array ThrustProfile and Times */
  /* Target thrusts are in TargetThrusts[30]. Want to compare every 0.5-second
   * element */
  /* in ThrustProfile array with next element of TargetThrusts, so long as */
  /* target thrusts are given every 0.5 seconds. Will be 0 after burnout at
   * shell */

  /* Record initial time, thrust, then every ThrustInterval steps */
  i = 0;
  lastitoprint = 0;

#ifdef PHASE2
  /* Now, reward DeltaV! */
  /* It seems to range from about 1900 to shy of 2500, so kluge up a reward */
  DeltaVReward = (DeltaV - 2350.) * 100.;
  if (DeltaVReward < 0.)
    DeltaVReward = 0.;
  if (Debug) {
    fprintf(outfp, "\nThroat area %0.14f, DeltaV %0.14f, DeltaVReward %0.14f",
            ThroatArea, DeltaV, DeltaVReward);
  }
#endif

  /* 02/28/19 putting reward for long burn as part of simultaneity reward, not
   * thrust reward */
  SimultaneityReward = Reward;
  Reward = 0.;

  /* THRUSTPROFILE REWARDS HERE */
  while (i < 29) { /* Looping through half-second time increment data */
    /* NOTE: Now continuing to reward thrusts after burnout has occurred. */
    /* Accumulate errors in thrusts from target thrusts */
    if (ThrustProfile[i] !=
        0.) { /* Don't want to reward thrusts when was not burning */
      Difference = fabs(ThrustProfile[i] - TargetThrusts[i]);
      if (Difference / TargetThrusts[i] < .10) {
        /* Thrust in desirable range, reward. */
        Reward +=
            1000. + 1000. / (Difference +
                             50.); /* some for in limits, more as get closer */
        if (Difference / TargetThrusts[i] < 0.04)
          Reward += 8000.;
        if (i <= 5 && Difference / TargetThrusts[i] < 0.05)
          Reward += 5000.;
        if (i == 0 && Difference / TargetThrusts[0] < .05)
          Reward += 55000.;
        if (i == 1 && Difference / TargetThrusts[1] < .05)
          Reward += 44000.;
        if (i == 2 && Difference / TargetThrusts[2] < .05)
          Reward += 33000.;
        if (i == 3 && Difference / TargetThrusts[3] < .05)
          Reward += 22000.;
        if (i == 4 && Difference / TargetThrusts[4] < .05)
          Reward += 31000.;
        if (i == 5 && Difference / TargetThrusts[5] < .05)
          Reward += 20000.;
        if (i > 5 && i <= 10 && Difference / TargetThrusts[i] < .05)
          Reward += 10000.;
        if (i > 10 && Difference / TargetThrusts[i] < 0.05)
          Reward += 5000.;
        // Just added below 02/02/19 to work on 8.5 second and beyond thrusts
        if (i > 16 && Difference / TargetThrusts[i] < 0.02)
          Reward += 10000.;
      }
      /* Now, as an alternative to killing the rocket, we just penalize for
       * pressure below */
      /* min allowable... this helps let the rocket approach MINPRESSURE, which
       * is must to */
      /* satisfy the first (declining 25g) thrust profile. */
      // CHECK THIS OUT NOW
      if (Time <= 10. && mPressureProfile[i] < MINPRESSURE)
        Reward *= 0.9;
    }
    lastitoprint = i;
    i++;
  }
  ThrustReward = Reward;
  Reward = 0.;
  if (Debug && ThrustProfile[0] > 0.) {
#ifdef EVOLVEPTAXOFSTAR
    fprintf(outfp, "\nFor PtAx's\n");
    for (SegNo = 0; SegNo < NCYLSEGS - 3; SegNo++) {
      fprintf(outfp, "%0.14f, ", PtAx[SegNo + 2]);
    }
#endif
    fprintf(
        outfp,
        "\nneval %ld's ThrustProfile      TargetThrusts at ThroatArea %0.14f\n",
        neval, ThroatArea);
    for (i = 0; i < lastitoprint; i++) {
      if (ThrustProfile[i] > 0.)
        fprintf(outfp, "%0.14f                  %0.14f\n", ThrustProfile[i],
                TargetThrusts[i]);
    }
    fprintf(outfp, "Thrust & Burntime Reward was %0.14f\n", ThrustReward);
  }
  /* Reward a solution in which all segments burn out (approximately) at the
   * same time. */
  /* Don't bother unless/until one segment reaches the shell. */
  RemainingDistanceToBurn = 0.; /* Total of all segs, local var. */
  SumSqRemDist = 0.;
  /* LET'S try biasing all segments toward being as close to the shell as
   * possible when rocket stops, */
  /* as long as the time is beyond 5, regardless of */
  /* the reason for stopping... deeper burns better, BUT reward them MORE the
   * closer time gets to 10. */
  /* Later, can rerun with longer burn time targets, but this gives a fixed,
   * common target to all for now. */
  for (SegNo = 0; SegNo < NTOTALSEGS; SegNo++) {
    RemainingDistToBurnInSeg[SegNo] = 0.;
    /* Now total up remaining unburned distances in all segs. */
    if (SegNo < NCYLSEGS) {
#ifdef STARGRAIN
      /* Find out if cyl seg burning rays */
      if (StarBurningThisSegNow[SegNo]) {
        /* Below is done only for output, since rocket is finished */
        SegmentRadius[SegNo] = 0.;
        for (k = 0; k < NRAYS + 2; k++) {
          RemainingDistToBurnInSeg[SegNo] +=
              (FinalBurnRadius[SegNo] - RayBurnDepth[k][SegNo]) / (NRAYS + 2);
          if (Debug || RayDebug)
            fprintf(outfp, "\nRayBurnDepth[%d][%d] %0.14f", k, SegNo,
                    RayBurnDepth[k][SegNo]);
        }
        /* At this point, RemDistToBurnInSeg is set; don't really need
         * SegmentRadius or use it. */
        SegmentRadius[SegNo] += RemainingDistToBurnInSeg[SegNo];
      } else {
#endif
        RemainingDistToBurnInSeg[SegNo] =
            fabs(FinalBurnRadius[SegNo] - SegmentRadius[SegNo]);

#ifdef STARGRAIN
      }
#endif
      SumSqRemDist +=
          RemainingDistToBurnInSeg[SegNo] * RemainingDistToBurnInSeg[SegNo];
      if (Debug)
        fprintf(outfp, "\nRemainingDistToBurnInSeg[%d] = %0.14f", SegNo,
                RemainingDistToBurnInSeg[SegNo]);
    }
    /* Now must handle corner segs and nozzle */
    /* Nozzle seg could still be burning cylindrically, or as rays */
    else if (SegNo == NTOTALSEGS - 1) {
      if (NozBurningRays) {
        for (k = 0; k < NNOZRAYS; k++) {
          RemainingDistToBurnInSeg[SegNo] +=
              (NozRayFinalBurnRadius[k] - NozRayRadius[k]) / 9.;
          if (Debug)
            fprintf(
                outfp,
                "\nNozRayFinalBurnRadius[%d] %0.14f, NozRayRadius[%d] %0.14f",
                k, NozRayFinalBurnRadius[k], k, NozRayRadius[k]);
        }
      } else {
        /* Looks strange, but since each nozzle ray has a different final
         * radius, must */
        /* calculate using them even though not now burning nozzle rays, and
         * average. */
        for (k = 0; k < NNOZRAYS; k++) {
          RemainingDistToBurnInSeg[SegNo] +=
              (NozRayFinalBurnRadius[k] - SegmentRadius[SegNo]) / 9.;
          if (Debug)
            fprintf(
                outfp,
                "\nNozRayFinalBurnRadius[%d] %0.14f, SegmentRadius[%d] %0.14f",
                k, NozRayFinalBurnRadius[k], SegNo, SegmentRadius[SegNo]);
        }
      }
      SumSqRemDist +=
          RemainingDistToBurnInSeg[SegNo] * RemainingDistToBurnInSeg[SegNo];
    }
    /* All that's left is corner segs */
    else {
      RemainingDistToBurnInSeg[SegNo] =
          fabs(FinalBurnRadius[SegNo] - SegmentRadius[SegNo]);
      SumSqRemDist +=
          RemainingDistToBurnInSeg[SegNo] * RemainingDistToBurnInSeg[SegNo];
      if (Debug)
        fprintf(outfp, "\nCorner seg [%d] RemDistToBurn %0.14f", SegNo,
                RemainingDistToBurnInSeg[SegNo]);
    }
    RemainingDistanceToBurn += RemainingDistToBurnInSeg[SegNo];
  }
  /* Calc sqrt(mean of squares minus square of mean)*/
  Variance =
      SumSqRemDist / NTOTALSEGS - ((RemainingDistanceToBurn / NTOTALSEGS) *
                                   (RemainingDistanceToBurn / NTOTALSEGS));
  StdDev = sqrt(Variance);
  /* Will start with segs even... so decrease for awhile.  Counteract that by
   * increasing it as */
  /* the MEAN value decreases. */
  /* I think I need to increase it much more sharply with burn distance, because
   * std dev factor */
  /* increases fast, decreasing reward fast. */
  /* OR, could make the std dev factor increase LESS rapidly, so less to
   * compensate for initially. */
  StdDevRemDist = StdDev - ((RemainingDistanceToBurn / NTOTALSEGS) *
                            (RemainingDistanceToBurn / NTOTALSEGS));

  if (Debug)
    fprintf(outfp, "\nStdDev %0.14f, StdDevRemDist = %0.14f", StdDev,
            StdDevRemDist);
  StdDevReward = ((100. / (StdDevRemDist + 0.01)) - 1000.);
  StdDevReward /= (RemainingDistanceToBurn / NTOTALSEGS + .03);
  if (StdDevReward > 20000.) {
    StdDevReward = (StdDevReward - 10000) * 10.;
  } else {
    StdDevReward = StdDevReward * 5.;
  }
  for (SegNo = 0; SegNo < NTOTALSEGS; SegNo++) {
    if (RemainingDistToBurnInSeg[SegNo] < 0.) {
      if (SimultaneityReward > 2000. * fabs(RemainingDistToBurnInSeg[SegNo])) {
        SimultaneityReward -= 1000. * fabs(RemainingDistToBurnInSeg[SegNo]);
      } else {
        SimultaneityReward -= SimultaneityReward / 2.;
      }
      if (Debug) {
        fprintf(outfp, "\nOverburn penalty %0.14f",
                1000. * fabs(RemainingDistToBurnInSeg[SegNo]));
      }
    }
  }
  SimultaneityReward += 5000. / (RemainingDistanceToBurn / 3. + 0.1);
  /* Add another term that starts small but gets big as approach optimum */
  SimultaneityReward += 500. / (RemainingDistanceToBurn / 3. + 0.001);

  if (Debug) {
    fprintf(outfp, "\nTotal RemainingDistanceToBurn = %0.14f",
            RemainingDistanceToBurn);
  }

  if (Debug && SimultaneityReward > 0. && gen % 10 == 9) {
    fprintf(outfp,
            "\nneval %ld, RemDistToBurn %0.14f, Time %0.14f, Term1 %0.14f, "
            "Term2 %0.14f, SimultaneityReward %0.14f",
            neval, RemainingDistanceToBurn, Time, Term1, Term2,
            SimultaneityReward);
  }
  if (WantThisGuyOutputFile) {
    for (SegNo = 0; SegNo < NTOTALSEGS; SegNo++) {
      BurnDistRemaining[SegNo] =
          (FinalBurnRadius[SegNo] - SegmentRadius[SegNo]);
    }
    //     WriteThisGuyOutputFile();
  }
  if (Debug) {
    fprintf(outfp,
            "\n*** Done with this rocket, time = %0.14f, fitness = %0.14f\n",
            Time, Reward);
  }
  return (Reward);
}

void WriteThisGuyOutputFile(void) {
  int i, SegNo;

  /* Output for this individual after evaluation, when called from HEEDS or
   * other optimizer */
  if (!(Evaluationfp = fopen("./Evaluation", "w"))) {
    printf("Output file %s cannot be written.\n", "Evaluation");
  }
  /* Write the time and thrust and desired thrust arrays into the graphing file,
   * 4 entries/line */
  /* StopCodes are : 0 normal termination (burn reached shell)
   *                 1 high pressure termination
   *                 2 low pressure termination
   *                 3 time limit reached (15 seconds)
   *                 4 abnormal situation (could not complete the simulation) */
  fprintf(Evaluationfp,
          "\n   Times,     Thrusts,    Desired Thrusts, Pressures, Reason for "
          "Termination, %d, evaluation number %ld, Fitness %0.14f\n",
          StopCode, neval, Reward);
  for (i = 0; i < 29; i++) {
    fprintf(Evaluationfp, "%e, %e, %e, %e\n", Times[i], ThrustProfile[i],
            TargetThrusts[i], mPressureProfile[i]);
  }
  fprintf(Evaluationfp, "\nNTOTALSEGS %d\nSegment Number, Unburned Depth",
          NTOTALSEGS);
  for (SegNo = 0; SegNo < NTOTALSEGS; SegNo++) {
    fprintf(Evaluationfp, "\n      %d               %0.14f\n", SegNo,
            (FinalBurnRadius[SegNo] - SegmentRadius[SegNo]));
  }
  fclose(Evaluationfp);
  Evaluationfp = NULL;
  WantThisGuyOutputFile = 0;
}

double CalcCornerSpecialPropellantMix(int SegNo) {
  /* CORNER SEGMENT PROCESSING ONLY, NOT STAR BURN */
  /* For explanation, see notes of 10/31/18 */
  /* For last 2 layers of a corner segment, equalizes the burnout of the two
   * boundaries by */
  /* dividing those 2 layers with a diagonal line from the beginning of the
   * shorter boundary */
  /* in next-to-last layer to point B along longer boundary, with type 4
   * propellant on the longer */
  /* side of that diagonal line until point B, and type 0 propellant on all of
   * the shorter */
  /* side. Then the two layers' boundaries will burn out at exactly the same
   * time. The propellant */
  /* types do NOT appear on the chromosome, since cannot be changed, but are
   * coded as type "0" */
  /* after chrom is decoded, but that values has no effect... the calculation is
   * just sent here. */
  /* This layer burn obviates the need to burn in rays--the two boundary rays
   * are sufficient. */
  double TempRefBurnRate, F, Diff, FOFB;

  int Debug = 0;

  /* Ignore the two propellant types listed, mix types 0 and 5, at least for
   * now. */
  /* Find out which boundary is longer. BdaryLen indexes from 0 as first CORNER
   * seg bdary. */
  F = PropellChoices[5] / PropellChoices[0];
  /* Calculate B, the distance along longer boundary to be mixing the fuel types
   * before */
  /* going to all type 0. */
  if (BdaryLen[SegNo - NCYLSEGS] > BdaryLen[SegNo - NCYLSEGS + 1]) {
    /* Lower index is longer. Figure out where "B" is on that boundary. */
    Diff = BdaryLen[SegNo - NCYLSEGS] - BdaryLen[SegNo - NCYLSEGS + 1];
    B[SegNo] = F * Diff / (F - 1.);
    if (Debug)
      fprintf(outfp, "\nF %0.14f, Diff %0.14f, B[%d] = %0.14f", F, Diff, SegNo,
              B[SegNo]);
  } else {
    /* Upper index is longer, work on that one. */
    Diff = BdaryLen[SegNo - NCYLSEGS + 1] - BdaryLen[SegNo - NCYLSEGS];
    B[SegNo] = F * Diff / (F - 1.);
    if (Debug)
      fprintf(outfp, "\nF %0.14f, Diff %0.14f, B[%d] = %0.14f", F, Diff, SegNo,
              B[SegNo]);
  }
  /* Must make sure changing to LAST layer from next-to-last does not change
   * this distance */
  FOFB =
      (SegmentRadius[SegNo] - LayerStartRadius[NLAYERS - 2][SegNo]) / B[SegNo];

  TempRefBurnRate = PropellChoices[4] * FOFB + PropellChoices[0] * (1. - FOFB);
  /* This returns a "mixture" burn rate, although the propellants are not
   * actually mixed, just the 2 layers */
  /* divided between the two propellant types by the diagonal line to point B.
   * The actual mixture will change */
  /* at each time step until the burn reaches B, after which it will be all type
   * 0 propellant. */
  if (Debug) {
    fprintf(outfp,
            "\nIn CalcCornerSpecialPropellantMix, seg %d, diff =%0.14f, "
            "F=%0.14f, B[%d]=%0.14f, FOFB=%0.14f, TempRefBurnRate=%0.14f\n",
            SegNo, Diff, F, SegNo, B[SegNo], FOFB, TempRefBurnRate);
    fflush(outfp);
  }
  return TempRefBurnRate;
}

/**********************************************************************************
 * Computes what propellant type is about to be burned, just using global
 *variables Returns value for RefBurnRate Modified to detect "special" burn of
 *last 2 layers of corner segments, 11/2/18.
 ***********************************************************************************/
double FindPropellantType(int SegNo) {
  int i, j;

  int Debug = 0;

  if (Debug && SegNo > 1 && SegNo < NCYLSEGS - 1) {
    fprintf(outfp, "\nneval = %ld", neval);
#ifdef STARGRAIN
    fprintf(outfp, "\nSTARGRAIN, and SegNo %d, StarBurningThisSegNow[SegNo] %d",
            SegNo, StarBurningThisSegNow[SegNo]);
#endif
    if (NozBurningRays)
      fprintf(outfp, "\nNozBurningRays %d", NozBurningRays);
  }
#ifdef STARGRAIN
  if (SegNo < 2 || SegNo > NCYLSEGS - 2 || !StarBurningThisSegNow[SegNo]) {
#else
  if (1 == 1) {
#endif
    /* Comes here if normal cyl burning or corner or nozzle burning */
    if (SegNo == NTOTALSEGS - 1 && !NozBurningRays &&
        SegmentRadius[NTOTALSEGS - 1] >= NozRayBurnStartRadius) {
      /* Here ONLY if in nozzle and just STARTING to burn rays */
      if (Debug) {
        fprintf(outfp, "\n\n\n!!! Nozzle Ray Burns has been triggered");
        fprintf(outfp,
                "\nSegmentRadius[Nozzle] %0.14f, NozRayBurnStartRadius %0.14f",
                SegmentRadius[NTOTALSEGS - 1], NozRayBurnStartRadius);
      }
      for (i = 0; i < NNOZRAYS; i++) {
        NozRayRadius[i] = NozRayBurnStartRadius;
      }
      if (Debug)
        fprintf(outfp, "\nNozRayBurnStartRadius %0.14f, NozRayRadius[0] %0.14f",
                NozRayBurnStartRadius, NozRayRadius[0]);
      NozBurningRays = 1;
      FindNozzleRayPropellants();
      /* Now NozRayRefBurnRate[i]'s have been set. */
      return -777.77777;
    } else if (SegNo == NTOTALSEGS - 1 && NozBurningRays) {
      FindNozzleRayPropellants();
      /* Here for nozzle already burning rays */
      return -777.77777;
    } else {
      /* NORMAL CYLINDER OR CORNER SEGMENT OR STAR SEGMENT AFTER CIRCULARIZED
       * PROCESSING */
      /* OR NOZZLE SEGMENT BEFORE RAYS */
      /* Flag StarBurningThisSegNow should have been zeroed when no longer
       * burning in rays. */
      /* HERE, handle special case of non-cyl (STAR, for example) burns, initial
       * part. */
      /* NORMAL processing, in corner or cylinder not currently (or ever) star
       * burning: determine layer from depth. */
      /* BIG DECISION:  Using SegmentRadius as current depth, find out whether
       * have crossed into */
      /* another propellant layer since the last timestep. Big timestep could
       * skip layers, so "while" */
      /* For now, if it crosses shell, letting it be caught elsewhere ... but
       * check it out */
      while ((CurrLayerIndex[SegNo] < NLAYERS - 1) &&
             (SegmentRadius[SegNo] >
              LayerStartRadius[CurrLayerIndex[SegNo] + 1][SegNo])) {
        CurrLayerIndex[SegNo]++;
      }
      /* Okay, we're now to burn in layer CurrLayerIndex[SegNo] */
      /* Need to handle CORNER segs separately if in last 2 layers with special
       * compensation */
      /* NEED also to handle NOZZLE when NOT burning rays and not starting to
       * burn rays, handled above */
      /* 3/5/19 EDG fixing that by excluding nozzle from special corner seg
       * handling below. */
      if (SegNo >= NCYLSEGS && SegNo != NTOTALSEGS - 1 &&
          CurrLayerIndex[SegNo] > NLAYERS - 3) {
        /* i.e., only in a Corner seg, last 2 layers */
        RDotRef[CurrLayerIndex[SegNo]][SegNo] =
            CalcCornerSpecialPropellantMix(SegNo);
        /* Putting in special return so won't overwrite RDotRef when burning
         * last 2 layers of corners. */
        /* Modified 3/05/19 EDG. */
        /* This value is NOT one of std propell types, but gotten by mixing them
         * on a per-ray basis, */
        /* as described elsewhere.  It must be set here and not altered
         * elsewhere until next timestep. */
        /* Don't need to save "evolved" original values because this individual
         * will never use them. */
        return -666.66666;
      }
      return RDotRef[CurrLayerIndex[SegNo]][SegNo];
    }
  }
#ifdef STARGRAIN
  else {
    /* StarBurningThisSegNow[SegNo] is SET, burning rays, so we're still in the
     * early layers of a star seg. */
    /* NEW circularizing logic implemented here. */
    /* Find the propellant type in which each RAY is now burning. */
    /* Must check for last of all rays RayBurnDepth[RayNo] from 0 to NRAYS+1,
     * crossing */
    /* InitBurnRadius[SegNo]. Only checking cylinders eligible for star burn. */
    for (j = 0; j < NRAYS + 2; j++) {
      if (NRAYSinPreCylLayer[SegNo] > 0) {
        /* Means are NOT YET compensating: that starts when last ray enters
         * cylinder. */
        if (RayInLayer[j][SegNo] == -1) {
          if (RayBurnDepth[j][SegNo] >= InitBurnRadius[SegNo]) {
            (NRAYSinPreCylLayer[SegNo])--;
            if (Debug)
              fprintf(outfp, "\nNRAYSinPreCylLayer[%d] = %d", SegNo,
                      NRAYSinPreCylLayer[SegNo]);
            RayInLayer[j][SegNo] = 0;
            RayRDotRef[j][SegNo] =
                PropellChoices[PropellType[RayInLayer[j][SegNo]][SegNo]];
            if (NRAYSinPreCylLayer[SegNo] == 0) {
              /* Ready NOW to begin burning CATCHUP propellant in all rays but
               * shortest, */
              /* as soon as they pass the CatchupLayer boundary. */
              /* Last ray just entered cylinder layer 0, but last to enter may
               * not be the */
              /* initially furthest ray--it may have jumped past others in the
               * layer so */
              /* determine the farthest ray and closest (to shell) shortest ray
               * now. */
              /* ONLY comes here once for each segment. */
              FarRayNo[SegNo] = FindFarRay(SegNo);
              ShortestRayNo[SegNo] = FindShortestRayNo(SegNo);
              if (Debug) {
                fprintf(outfp,
                        "\n\n\n*** Started Compensation process. In SegNo %d, "
                        "shortest ray is %d, farthest ray is %d",
                        SegNo, ShortestRayNo[SegNo], FarRayNo[SegNo]);
              }
              /* Segment-wide marker */
              CompensatingNow[SegNo] = 1;
              NRaysToCatchUp[SegNo] =
                  NRAYS + 2; /* Now counting the shortest ray itself, too. */
              /* Decide on the layer number that will trigger catchup burns for
               * all rays. */
              /* Need to find the layer the currently shortest ray is in, and
               * add one to it to */
              /* get the boundary for triggering catchup burn of all other rays.
               */
              StartCatchupLayer[SegNo] =
                  RayInLayer[ShortestRayNo[SegNo]][SegNo] + 1;
              /* KLUGE HERE to speed up compensation... slow down any layers
               * that have caught up */
              /* including the shortest ray, which will be caught up as soon as
               * it passes into */
              /* the next layer. */

              if (Debug)
                fprintf(outfp,
                        "\nIn seg %d, NRaysToCatchUp is %d, StartCatchupLayer "
                        "is %d, depth %0.14f",
                        SegNo, NRaysToCatchUp[SegNo], StartCatchupLayer[SegNo],
                        LayerStartRadius[StartCatchupLayer[SegNo]][SegNo]);
            }
          } else {
            if (RayBurnDepth[j][SegNo] < StarRayOuterInnerBdary[SegNo]) {
              RayRDotRef[j][SegNo] = PreCylStarRDotRef[0][SegNo];
            } else {
              RayRDotRef[j][SegNo] = PreCylStarRDotRef[1][SegNo];
            }
            if (Debug)
              fprintf(outfp, "\nPreCyl RayRDotRef[%d][%d] = %0.14f", j, SegNo,
                      RayRDotRef[j][SegNo]);
          }
        } else {
          /* Some rays still in pre-cyl, but not this one, so assign normal
           * propellant, but */
          /* first check to see if it passed a layer boundary, advance its layer
           * if it did. */
          if (RayBurnDepth[j][SegNo] >=
              LayerStartRadius[RayInLayer[j][SegNo] + 1][SegNo]) {
            (RayInLayer[j][SegNo])++;
          }
          RayRDotRef[j][SegNo] =
              PropellChoices[PropellType[RayInLayer[j][SegNo]][SegNo]];
        }
      } else {
        /* There are no rays burning in pre-cyl, so we're either compensating
         * now or done with it. */
        /* CatchingUpRay will have 3 values:  0: not yet burning catchup
         * propellant, which means this */
        /* ray has not yet crossed StartCatchupLayer start boundary; 1: we are
         * now burning it, */
        /* and 2: completed burning it, so it has caught up. Was initialized to
         * 0 at initialization. */
        /* First track the layer the ray is burhing in, to be sure is correct
         * when stops catching up. */
        if (RayBurnDepth[j][SegNo] >=
            LayerStartRadius[RayInLayer[j][SegNo] + 1][SegNo]) {
          (RayInLayer[j][SegNo])++;
        }
        if (CompensatingNow[SegNo]) {
          if (CatchingUpRay[j][SegNo] < 2) {
            if (RayBurnDepth[j][SegNo] >=
                LayerStartRadius[StartCatchupLayer[SegNo]][SegNo]) {
              /* In or beyond StartCatchupLayer, so burn Catchup unless caught
               * up. */
              if (RayBurnDepth[j][SegNo] <
                  RayBurnDepth[ShortestRayNo[SegNo]][SegNo]) {
                /* Is in process of catching up, burning catchup propellant. */
                if (CatchingUpRay[j][SegNo] == 0)
                  CatchingUpRay[j][SegNo] = 1;
                /* Mark that is in process of catching up. Has NOT YET caught
                 * up. */
                LatestCatchupDepthBurned[j][SegNo] = RayBurnDepth[j][SegNo];
                //		        RayRDotRef[j][SegNo] =
                // PropellChoices[StarRayCatchupPropellType[SegNo]];
                RayRDotRef[j][SegNo] =
                    PropellChoices[RDotRefCatchupInput[SegNo - 2]];
                if (Debug)
                  fprintf(outfp,
                          "\nRay[%d][%d] is burning RayRDotRef[%d][%d] %0.14f",
                          j, SegNo, j, SegNo, RayRDotRef[j][SegNo]);
              } else {
                /* Since was not already marked as 2, this ray HAS just caught
                 * up. */
                /* It may be the shortest ray itself, just crossing this
                 * boundary, so should */
                /* start burning SLOWEST propellant, like any other "caught up"
                 * ray. */
                /* For IntegralContributions to work, need to track RayInLayer
                 * even when it is */
                /* burning catchup propellant. */
                /* Force caught up rays to burn type 0 while seg is still
                 * compensating. */
                RayRDotRef[j][SegNo] = PropellChoices[0];
                /* Recording this as last depth reached by burning Catchup
                 * propellant */
                LatestCatchupDepthBurned[j][SegNo] = RayBurnDepth[j][SegNo];
                if (CatchingUpRay[j][SegNo] < 2) {
                  /* This ray has JUST caught up with (or is) the shortest ray,
                   * so it can stop */
                  /* burning catchup */
                  CatchingUpRay[j][SegNo] = 2;
                  NRaysToCatchUp[SegNo]--;
                  if (Debug) {
                    fprintf(
                        outfp,
                        "\nRay %d just caught up, now SLOW RayRDotRef[%d][%d] "
                        "%0.14f, NRaysToCatchUp[%d] = %d",
                        j, j, SegNo, RayRDotRef[j][SegNo], SegNo,
                        NRaysToCatchUp[SegNo]);
                  }
                }
                if (NRaysToCatchUp[SegNo] == 0) {
                  /* Last ray of seg just caught up */
                  EndCatchUp(j, SegNo);
                  /* CurrLayerIndex[SegNo] has now been set; that's all we need
                   */
                  goto EndOfSeg;
                }
              }
            } else {
              /* Hasn't yet reached StartCatchupLayer, so burn normal propellant
               */
              RayRDotRef[j][SegNo] =
                  PropellChoices[PropellType[RayInLayer[j][SegNo]][SegNo]];
            }
          } else {
            /* CatchingUpRay is 2, so was already caught up, but still
             * compensating, so burn type 0 */
            /* Instead of burning normal fuel, have all "caught up" rays burn
             * slowest fuel type. */
            RayRDotRef[j][SegNo] = PropellChoices[0];
            if (Debug) {
              fprintf(outfp,
                      "\n *** Ray %d caught up earlier, burning "
                      "RayRDotRef[%d][%d] = %0.14f",
                      j, j, SegNo, RayRDotRef[j][SegNo]);
            }
          }
        } else {
          /* NOT compensating, so burning normal propellant */
          RayRDotRef[j][SegNo] =
              PropellChoices[PropellType[RayInLayer[j][SegNo]][SegNo]];
        }
        /* Could be compensating now or not */
        if (RayBurnDepth[j][SegNo] >=
            LayerStartRadius[RayInLayer[j][SegNo] + 1][SegNo]) {
          /* This ray already in cylinder, just crossed some other layer start
           * radius */
          (RayInLayer[j][SegNo])++;
          if (RayInLayer[j][SegNo] > NLAYERS - 1) {
            /* Burn has hit the shell in this segment, must stop the burn */
            StopCode = 0;
            Time += DeltaTime;
            return -222.22222; /* NEED to calculate fitness, this rocket is done
                                  burning, */
            /* reached shell.  Return a signal that signifies that. */
          } else {
            if (RayBurnDepth[j][SegNo] >=
                RayBurnDepth[ShortestRayNo[SegNo]][SegNo]) {
              /* May still be in compensating mode, but THIS ray has passed
               * shortest ray. */
              if (!CompensatingNow[SegNo]) {
                RayRDotRef[j][SegNo] =
                    PropellChoices[PropellType[RayInLayer[j][SegNo]][SegNo]];
              }
              /* (Otherwise, will stay at min type 0 as already set) */
            }
          }
        } else if (CatchingUpRay[j][SegNo] == 0) {
          /* Not in pre-cyl layer, didn't just cross a boundary. Could be in
           * layer 0 or */
          /* or higher before burning catchup, but not during catchup or after
           * catchup. */
          RayRDotRef[j][SegNo] =
              PropellChoices[PropellType[RayInLayer[j][SegNo]][SegNo]];
        }
      } /* end of loop when there are no rays burning in pre-cyl propellant */
    }   /* Loop over rays ends here. */
    /* Returning a special code to tell that instead of returning a reference
     * burn rate for the whole */
    /* segment, we returned ray-specific RayRDotRef values in global variables,
     * so the return value is to be ignored. */
    if (Debug) {
      for (j = 0; j < NRAYS + 2; j++) {
        fprintf(outfp, "\nRayDepth [%d][%d] = %0.14f, RayRDotRef = %0.14f", j,
                SegNo, RayBurnDepth[j][SegNo], RayRDotRef[j][SegNo]);
      }
    }
    return -777.77777;
  }
EndOfSeg:
#endif
  return RDotRef[CurrLayerIndex[SegNo]][SegNo];
}

#ifdef STARGRAIN
void EndCatchUp(int j, int SegNo) {
  int i, k;

  int Debug = 0;

  /* Burn just went beyond compensation layer, should now be roughly
   * circularized. Turn off flag saying it */
  /* needs to be burned by rays any more. Set the burn depth of the central ray
   * (the segment's burndepth) */
  /* to the average of those of all the rays when this depth is reached. */

  CompensatingNow[SegNo] = 0;
  if (Debug) {
    fprintf(outfp,
            "\n***Star burn was just unified in segment %d, layer #%d, turn "
            "off special star burn. *** \n",
            SegNo, j);
  }
  SegmentRadius[SegNo] = 0.;
  for (k = 0; k < NRAYS + 2; k++) {
    SegmentRadius[SegNo] += RayBurnDepth[k][SegNo];
    if (Debug) {
      fprintf(outfp,
              "\nNOTE:  CHECK RAY DIFFERENCES: Compensation Finished, Final "
              "Ray burn depth %d = %0.14f",
              k, RayBurnDepth[k][SegNo]);
    }
  }

  SegmentRadius[SegNo] /= (NRAYS + 2);
  NextSegRadius[SegNo] =
      SegmentRadius[SegNo]; /* THIS is where next radius is kept... */
  /* Prepare to cut 'er loose as a normal burn from here out */
  StarBurningThisSegNow[SegNo] = 0;
  CurrLayerIndex[SegNo] = 0;
  while ((CurrLayerIndex[SegNo] < NLAYERS - 1) &&
         (SegmentRadius[SegNo] >
          LayerStartRadius[CurrLayerIndex[SegNo] + 1][SegNo])) {
    (CurrLayerIndex[SegNo])++;
  }
  if (Debug) {
    fprintf(outfp, "\nInEndCatchUp(), CurrLayerIndex[%d] = %d", SegNo,
            CurrLayerIndex[SegNo]);
    fprintf(
        outfp,
        "\nSegment %d now burning evolved propellant RDotRef[%d][%d] %0.14f",
        SegNo, CurrLayerIndex[SegNo], SegNo,
        RDotRef[CurrLayerIndex[SegNo]][SegNo]);
  }
  /* Now know what layer this new arc is burning in */
  return;
}

int FindShortestRayNo(int SegNo) {
  int i;
  int ShortestNo = 0;
  double ShortestDist;

  int Debug = 0;

  ShortestNo = 0;
  ShortestDist = RayBurnDepth[0][SegNo];
  for (i = 1; i < NRAYS + 2; i++) {
    if (RayBurnDepth[i][SegNo] > ShortestDist) {
      ShortestNo = i;
      ShortestDist = RayBurnDepth[i][SegNo];
    }
  }
  if (Debug) {
    fprintf(outfp, "\nShortest Ray number %d", ShortestNo);
  }
  return ShortestNo;
}

int FindFarRay(int SegNo) {
  int i;
  int FurthestRayNo = 0;
  double FurthestDist;

  int Debug = 0;

  FurthestRayNo = 0;
  FurthestDist = RayBurnDepth[0][SegNo];
  for (i = 1; i < NRAYS + 2; i++) {
    if (RayBurnDepth[i][SegNo] < FurthestDist) {
      FurthestRayNo = i;
      FurthestDist = RayBurnDepth[i][SegNo];
    }
  }
  if (Debug) {
    fprintf(outfp, "\nFurthestRayNo %d", FurthestRayNo);
  }
  return FurthestRayNo;
}
#endif /* of STARGRAIN */

/*********************************************************************************
 * Newton solver.
 * Returns new TotalPressure value
 * @param aSegmentArea current chamber area
 * @param TotalPressure total pressure at current time step
 * @param ThroatArea current throat area
 * Revised Oct 16, 2018, EDG, to sum up several segments' contributions to the
 *integral. It happens in residual(), which the Newton function calls.
 **********************************************************************************/
double newton(double lTotalPressure, double lThroatArea) {
  _Bool tDone = 0;
  size_t tIteration = 1;

  /* Define a new local variable NewTotalPressure with value set by passed
   * parameter, to manipulate locally. */
  double tMyResidualEvaluation;
  double tMyJacobianEvaluation;
  double residual();
  double jacobian();
  double LocalTotalPressure;

  LocalTotalPressure = lTotalPressure;
  while (tDone == 0) {
    tMyResidualEvaluation = residual(LocalTotalPressure, lThroatArea);
    tMyJacobianEvaluation = jacobian(LocalTotalPressure, lThroatArea);
    if (tMyJacobianEvaluation == -9999.) {
      fprintf(outfp,
              "\n****** tMyJacobianEvaluation was NaN, reset to -9999.\n");
      StopCode = 4;
      return (tMyJacobianEvaluation);
    }
    if (fabs(tMyJacobianEvaluation) < 1.e-60) {
      if (tMyJacobianEvaluation < 0.) {
        tMyJacobianEvaluation = -1.e-60;
      } else {
        tMyJacobianEvaluation = 1.e-60;
      }
      fprintf(outfp,
              "\n *** *** Jacobian was nearly zero, %e , reset away so don't "
              "divide by 0.\n",
              tMyJacobianEvaluation);
    }
    DeltaPressure = -1. * tMyResidualEvaluation / tMyJacobianEvaluation;
    if (LocalTotalPressure + DeltaPressure < 0.) {
      return (-444.44444); /* code for Jacobian failed */
    }
    LocalTotalPressure +=
        DeltaPressure; /* UPDATING the local var, not the pressure passed in. */
    tIteration += 1;
    if (fabs(DeltaPressure) < NewtonTolerance || tIteration > MaxNumNewtonItr) {
      tDone = 1;
    }
  } /* end of while tDone loop */
  return (LocalTotalPressure);
}

/*********************************************************************************
 * @brief Jacobian evaluation.
 * @param TotalPressure total pressure at current time step
 * @param ThroatArea current throat area
 **********************************************************************************/
double jacobian(double lTotalPressure, double lThroatArea) {
  double Value;
  int i;

  /* Must calculate first term per segment and sum them, then subtract second
     term to get proper gradient
     Are calculating gradient of (mdotin - mdotout); mdotin is integral (so sum)
     over segments
     BUT the segment area x burn rates for all segments, adjusted for lateral
     burns, etc., are already
     summed in IntegralContribs, so should just use those, not sum here by segments.
  */

  Value = 0.;
  Value += PropellantDensity * IntegralContribs * Alpha * InvPrefAlpha *
           pow(lTotalPressure, Alphaminus1);
  Value += (-lThroatArea / CharacteristicVelocity);
  if (Value != Value) {
    /* means result was NaN not a number */
    StopCode = 4;
    return (-9999.);
  }
  return Value;
}

/*********************************************************************************
 * Residual evaluation
 * EDG: This function returns mdotin - mdotout, for use in newton algorithm
 * to determine actual pressure Pc to enable calculating rdot(x,y,z).
 * @param aSegmentArea current chamber area
 * @param TotalPressure total pressure at current time step
 * @param ThroatArea current throat area
 * Revised 10/18/18 by EDG so can handle burn of several differently configured
 *segments.
 **********************************************************************************/
double residual(double lTotalPressure, double lThroatArea) {
  double Value = 0.;

  Value += PropellantDensity * IntegralContribs * InvPrefAlpha *
           pow(lTotalPressure, Alpha);
  Value += InsulDensity * InsulBurnArea * InsulRefBurnRate *
           pow(lTotalPressure / RefPressure, 0.8);
  Value -= lThroatArea * lTotalPressure / CharacteristicVelocity;
  return Value;
}
