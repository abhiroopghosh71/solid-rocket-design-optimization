/* test_rocket.c v1.07
Author: Abhiroop Ghosh
Last modified: Sep 18, 2019

This file provides a file I/O interface to the solid-rocket simulator. This is useful if
a file I/O based optimization software like HEEDS is being used. It also allows the user
to treat the simulator as a black box model without being concerned about the correct interfacing.
*/
#include "approcket.c"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

FILE* outputfp;

int main(int argc, char *argv[])
{
    char* input_file;
    char* output_file;
    char* radius_file;

    char* INPUT_FILE_DEFAULT = "input_vec";
    char* OUTPUT_FILE_DEFAULT = "output_vec";
    char* LAYER_RADIUS_DEFAULT_FILE = "layer_radius";
    
    // Read inputs and write outputs to the files specified by the user through command-line arguments
    // First argument is the input file name
    // Second argument is the output file name
    // Third argument is the start radii of each layer for all the segments
    // If number of args is less than 3, the missing ones are set to default values
    if(argc > 3)
    {
        input_file = argv[1];
        output_file = argv[2];
        radius_file = argv[3];
    }
    else if(argc == 3)
    {
        input_file = argv[1];
        output_file = argv[2];
        radius_file = LAYER_RADIUS_DEFAULT_FILE;
        printf("if 1\n");
    }
    else if(argc == 2)
    {
        input_file = argv[1];
        output_file = OUTPUT_FILE_DEFAULT;
        radius_file = LAYER_RADIUS_DEFAULT_FILE;
        printf("if 2\n");
    }
    else
    {
        input_file = INPUT_FILE_DEFAULT;
        output_file = OUTPUT_FILE_DEFAULT;
        radius_file = LAYER_RADIUS_DEFAULT_FILE;
        printf("if 3\n");
    }
    printf("Input file name %s\n", input_file);
    printf("Output file name %s\n", output_file);
    printf("Radius file name %s\n", radius_file);

   int    tReasonStopped;
   double tThrustReward, tSimultaneityReward;

    // Array to store data at 0.5 second intervals till burnout time
   double tTimes[NumTimePointsCoarse];
   double tThrustProfile[NumTimePointsCoarse];
   double tPressureProfile[NumTimePointsCoarse];

   // Array to store data at 0.025 second intervals till burnout time
   double tTimesFine[NumTimePointsFine];
   double tThrustProfileFine[NumTimePointsFine];
   double tPressureProfileFine[NumTimePointsFine];

   double tBurnDistRemaining[NTOTALSEGS];
   double tThroatAreaInput = .001297; /* JUST A LEGAL EXAMPLE */
   double tDeltaV = 0.;
   double tMaxThrust = 7532.6;

   double tTargetThrusts[] = {7532.6, 7190.6, 6864.1, 6552.7, 6255.1, 5971.3, 5700.0, 5441.1, 5194.2, 4958.4, 4733.4,
                                4518.5, 4313.4, 4117.3, 3930.4, 3752.1, 3581.7, 3419.3, 3264.1, 3115.5, 2974.1, 2839.6,
                                2711.2, 2588.7, 2471.7, 2360.0, 2253.3, 2151.4, 2054.2, 1961.3};
//   double tTargetThrusts[] = {5000., 5000., 5000., 5000., 5000., 5000., 5000., 5000., 5000.,
//                                5000., 5000., 5000., 5000., 5000., 5000., 5000., 5000., 5000., 5000., 5000., 5000.,
//                                5000., 5000., 5000., 5000., 5000., 5000., 5000., 5000., 5000.};
//   double tTargetThrusts[] = {4906.31189966839, 4906.31189966839, 4906.31189966839, 4906.31189966839, 4906.31189966839,
//                              4906.31189966839, 4906.31189966839, 4906.31189966839, 4906.31189966839, 4906.31189966839,
//                              4906.31189966839, 4906.31189966839, 4906.31189966839, 4906.31189966839, 4906.31189966839,
//                              4906.31189966839, 4906.31189966839, 4906.31189966839, 4906.31189966839, 4906.31189966839};

    // Store actual depths of the rays in the star portion
   int tRayBurnDepths[(NCYLSEGS-3)*(NRAYS+1)];
   for(int i = 0; i < (NCYLSEGS-3)*(NRAYS+1); i++)
   {
        tRayBurnDepths[i] = 0;
   }
   
   FILE* inputfp;
   FILE* radiusfp;
   if (!(inputfp = fopen(input_file, "r")))
   {
       printf("Input file %s cannot be read.\n", INPUT_FILE_DEFAULT);
       exit(0); 
   }
   if (!(radiusfp = fopen(radius_file, "r")))
   {
       printf("Input file %s cannot be read.\n", LAYER_RADIUS_DEFAULT_FILE);
       exit(0); 
   }
   
   int tPropellantToEval[NTOTALSEGS*NLAYERS];
   double tLayerStartRadius[NTOTALSEGS*NLAYERS];
   int tRayBurnDepthInput[(NCYLSEGS-3)*(NRAYS+1)];
   int tPreCylStarRDotRefInput[NCYLSEGS-3];
   int tRDotRefCatchupInput[NCYLSEGS-3];
   char line[100];
   
   for(int i = 0; i < (NTOTALSEGS*NLAYERS); i++)
   {
        //printf("Here %d\n", i);
       fgets(line, sizeof(line), inputfp);
       tPropellantToEval[i] = (int) atof(line);
   }
   
   for(int i = 0; i < (NTOTALSEGS*NLAYERS); i++)
   {
        //printf("Here %d\n", i);
       fgets(line, sizeof(line), radiusfp);
       tLayerStartRadius[i] = atof(line);
   }
   fclose(radiusfp);
   
   for(int i = 0; i < (NCYLSEGS - 3); i++)
   {
        //printf("Here %d\n", i);
        for(int j = 0; j < NRAYS+1; j++)
        {
           fgets(line, sizeof(line), inputfp);
           tRayBurnDepthInput[i*(NRAYS+1) + j] = (int) atof(line);
       }
       fgets(line, sizeof(line), inputfp);
       tPreCylStarRDotRefInput[i] = (int) atof(line);
       fgets(line, sizeof(line), inputfp);
       tRDotRefCatchupInput[i] = (int) atof(line);
   }  
   fclose(inputfp);

   int tSegStarBurnStatus[NTOTALSEGS];
   double tTimestarBurnFinish[NTOTALSEGS];
   double tSegResidualPerTimestepCoarse[NTOTALSEGS][NumTimePointsCoarse];
   double tSegResidualPerTimestepFine[NTOTALSEGS][NumTimePointsFine];
   double tSegBurnAreaPerTimestepCoarse[NTOTALSEGS][NumTimePointsCoarse];
   double tSegBurnAreaPerTimestepFine[NTOTALSEGS][NumTimePointsFine];
   double tBurnLayerPerTimestepCoarse[NTOTALSEGS][NumTimePointsCoarse];
   double tBurnLayerPerTimestepFine[NTOTALSEGS][NumTimePointsFine];

   // Write to output file
   FILE* outputfp;
   if (!(outputfp = fopen(output_file, "w")))
   {
       printf("Output file %s cannot be written.\n", "output_vec"); 
       exit(0);
   }

   // int tActualRayBurnDepths[] = {1, 1, 2, 7, 8, 15, 2, 5, 9, 14, 15, 15, 1, 5, 10, 15, 15, 15};
    int rayDepthFlag = 0;  // If 1, means the exact ray depths are being supplied instead of the codes
   int tActualRayBurnDepths[] = {1, 1, 2, 7, 8, 15, 2, 5, 9, 14, 15, 15, 1, 5, 10, 15, 15, 15};

   /* Modify the for loop to run the simulator for the same inputs multiple tTimes.
   This tests the simulator for a bug where same run gives different results.*/
   for(int loop_count = 1; loop_count<=3; loop_count++)
   {
        // Call the simulator

        if (rayDepthFlag == 0)
        {
        objectivefunction(&tReasonStopped, &tThrustReward, &tSimultaneityReward, tTimes, tThrustProfile, tPressureProfile,
                          tTimesFine, tThrustProfileFine, tPressureProfileFine,
                          tBurnDistRemaining, tPropellantToEval, tLayerStartRadius, tRayBurnDepthInput, tPreCylStarRDotRefInput,
                          tRDotRefCatchupInput, tThroatAreaInput, &tDeltaV, tMaxThrust, tRayBurnDepths, tSegStarBurnStatus,
                          tTimestarBurnFinish, tSegResidualPerTimestepCoarse, tSegResidualPerTimestepFine, tTargetThrusts,
                          rayDepthFlag, tSegBurnAreaPerTimestepCoarse, tSegBurnAreaPerTimestepFine,
                          tBurnLayerPerTimestepCoarse, tBurnLayerPerTimestepFine);
      }
      else if (rayDepthFlag == 1)
      {
        objectivefunction(&tReasonStopped, &tThrustReward, &tSimultaneityReward, tTimes, tThrustProfile, tPressureProfile,
                          tTimesFine, tThrustProfileFine, tPressureProfileFine,
                          tBurnDistRemaining, tPropellantToEval, tLayerStartRadius, tActualRayBurnDepths, tPreCylStarRDotRefInput,
                          tRDotRefCatchupInput, tThroatAreaInput, &tDeltaV, tMaxThrust, tRayBurnDepths, tSegStarBurnStatus,
                          tTimestarBurnFinish, tSegBurnAreaPerTimestepCoarse, tSegResidualPerTimestepFine, tTargetThrusts,
                          rayDepthFlag, tSegBurnAreaPerTimestepCoarse, tSegBurnAreaPerTimestepFine,
                          tBurnLayerPerTimestepCoarse, tBurnLayerPerTimestepFine);
      }
      
      printf("\nThrust profile fine\n");
      for (int i = 0; i < 500; i++)
      {
            printf("%f,     %f\n", tTimesFine[i], tThrustProfileFine[i]);
        }

        // Print output to stdout
        printf("\ntTimes(s)  Thrust (N)  TargetThrust (N)   Pressure (Pa)\n");
        for(int i = 0; i < 21; i++)
        {
           printf("%f  %f  %f  %f\n", tTimes[i], tThrustProfile[i], TargetThrusts[i], tPressureProfile[i]);
        }

        printf("SegNo   Burn Dist. Remaining (m)\n");
        for(int i = 0; i < NTOTALSEGS; i++)
        {
           printf("%d       %f\n", i, tBurnDistRemaining[i]);
        }

        printf("Precyl Inputs\n");
        for(int i = 0; i < (NCYLSEGS-3); i++)
        {
            printf("%d ", tPreCylStarRDotRefInput[i]);
        }
        printf("\n");
        
        printf("RDotRef Catchup Inputs\n");
        for(int i = 0; i < (NCYLSEGS-3); i++)
        {
            printf("%d ", tRDotRefCatchupInput[i]);
        }
        printf("\n");
        
        printf("Ray Burn Depth Inputs\n");
        for(int i = 0; i < (NCYLSEGS-3)*(NRAYS+1); i++)
        {
            printf("%d ", tRayBurnDepthInput[i]);
        }
        printf("\n");

        printf("Actual Ray Burn Depths\n");
        for(int i = 0; i < (NCYLSEGS-3)*(NRAYS+1); i++)
        {
            printf("%d ", tRayBurnDepths[i]);
        }
        printf("\n");
        
        printf("Segment wise star-burn status\n");
        for(int i = 0; i < NTOTALSEGS; i++)
        {
            printf("%d ", tSegStarBurnStatus[i]);
        }
        printf("\n");

        printf("Segment-wise star circularization time\n");
        for(int i = 0; i < NTOTALSEGS; i++)
        {
            printf("%f ", tTimestarBurnFinish[i]);
        }
        printf("\n");

       printf("\nThrustReward %f", tThrustReward);
       printf("\nSimultaneityReward %f", tSimultaneityReward);

       printf("\nThroat Area %f", tThroatAreaInput);
       printf("\ntDeltaV %f", tDeltaV);
       printf("\nReason Stopped %d\n", tReasonStopped);

       /* StopCodes are : 0 normal termination (burn reached shell)
           *                 1 high pressure termination
           *                 2 low pressure termination
           *                 3 time limit reached (15 seconds)
           *                 4 abnormal situation (could not complete the simulation) */
           
           
        // Write to output file
        fprintf(outputfp, "%d\n", tReasonStopped);
        fprintf(outputfp, "%f\n", tThrustReward);
        fprintf(outputfp, "%f\n", tSimultaneityReward);
        fprintf(outputfp, "%f\n", tThroatAreaInput);
        fprintf(outputfp, "%f\n", tDeltaV);

        for(int i = 0; i < NTOTALSEGS; i++)
        {
            for(int j = 0; j < NumTimePointsFine; j++)
            {
                fprintf(outputfp, "%f ", tSegBurnAreaPerTimestepCoarse[i][j]);
            }
            fprintf(outputfp, "\n");
        }

        for(int i = 0; i < NumTimePointsCoarse; i++)
        {
            fprintf(outputfp, "%f ", tTimes[i]);
        }
        fprintf(outputfp, "\n");

        for(int i = 0; i < NumTimePointsCoarse; i++)
        {
            fprintf(outputfp, "%f ", tThrustProfile[i]);
        }
        fprintf(outputfp, "\n");

        for(int i = 0; i < NumTimePointsCoarse; i++)
        {
            fprintf(outputfp, "%f ", tPressureProfile[i]);
        }
        fprintf(outputfp, "\n");

        for(int i = 0; i < NTOTALSEGS; i++)
        {
           fprintf(outputfp, "%f ", tBurnDistRemaining[i]);
        }
        fprintf(outputfp, "\n");
        
        for(int i = 0; i < (NCYLSEGS-3)*(NRAYS+1); i++)
        {
            fprintf(outputfp, "%d ", tRayBurnDepthInput[i]);
        }
        fprintf(outputfp, "\n");

        for(int i = 0; i < (NCYLSEGS-3)*(NRAYS+1); i++)
        {
            fprintf(outputfp, "%d ", tRayBurnDepths[i]);
        }
        fprintf(outputfp, "\n");

        for(int i = 0; i < NTOTALSEGS; i++)
        {
            fprintf(outputfp, "%d ", tSegStarBurnStatus[i]);
        }
        fprintf(outputfp, "\n");

        // Print data in much finer time interval
        for(int i = 0; i < NumTimePointsFine; i++)
        {
            fprintf(outputfp, "%f ", tTimesFine[i]);
        }
        fprintf(outputfp, "\n");

        for(int i = 0; i < NumTimePointsFine; i++)
        {
            fprintf(outputfp, "%f ", tThrustProfileFine[i]);
        }
        fprintf(outputfp, "\n");

        for(int i = 0; i < NumTimePointsFine; i++)
        {
            fprintf(outputfp, "%f ", tPressureProfileFine[i]);
        }
        fprintf(outputfp, "\n");
        /*
        for(int i = 0; i < NTOTALSEGS; i++)
        {
            for(int j = 0; j < NumTimePointsCoarse; j++)
            {
                fprintf(outputfp, "%f ", tSegResidualPerTimestepCoarse[i][j]);
            }
            fprintf(outputfp, "\n");
        }
        fprintf(outputfp, "\n");*/

        //fprintf(outputfp, "\n");
        for(int i = 0; i < NTOTALSEGS; i++)
        {
            for(int j = 0; j < NumTimePointsFine; j++)
            {
                fprintf(outputfp, "%f ", tSegResidualPerTimestepFine[i][j]);
            }
            fprintf(outputfp, "\n");
        }
        fprintf(outputfp, "\n");

        /*
        for(int i = 0; i < NTOTALSEGS; i++)
        {
            for(int j = 0; j < NumTimePointsCoarse; j++)
            {
                fprintf(outputfp, "%f ", tSegBurnAreaPerTimestepCoarse[i][j]);
            }
            fprintf(outputfp, "\n");
        }
        fprintf(outputfp, "\n");*/

        //fprintf(outputfp, "\n");

        for(int i = 0; i < NTOTALSEGS; i++)
        {
            for(int j = 0; j < NumTimePointsFine; j++)
            {
                fprintf(outputfp, "%f ", tSegBurnAreaPerTimestepFine[i][j]);
            }
            fprintf(outputfp, "\n");
        }
        fprintf(outputfp, "\n");

        for(int i = 0; i < NTOTALSEGS; i++)
            {
                for(int j = 0; j < NumTimePointsFine; j++)
                {
                    fprintf(outputfp, "%f ", tBurnLayerPerTimestepFine[i][j]);
                }
                fprintf(outputfp, "\n");
            }
            fprintf(outputfp, "\n");


            for(int i = 0; i < NTOTALSEGS; i++)
            {
                for(int j = 0; j < NumTimePointsFine; j++)
                {
                    fprintf(outputfp, "%f ", tBurnLayerPerTimestepFine[i][j]);
                }
                fprintf(outputfp, "\n");
            }
            fprintf(outputfp, "\n");

    }



    return 0;
}
