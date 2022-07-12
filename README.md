# Solid rocket motor design using multi-objective optimization
Multi-objective optimization of a solid rocket motor [1] to fit a pre-specified thrust profile and minimize insulation. 
This package uses Pymoo [2] for running NSGA-II [3].

<h2>Note to macOS and Windows users</h4>
The code currently may not work properly on MAC-OS. One known issue is when you run ```make``` it might say:

```xcrun: error: invalid active developer path (/Library/Developer/CommandLineTools), missing xcrun at: /Library/Developer/CommandLineTools/usr/bin/xcrun```

Try to follow the instructions given in https://apple.stackexchange.com/questions/254380/why-am-i-getting-an-invalid-active-developer-path-when-attempting-to-use-git-a but its not guaranteed to fix the issue.

The code should run properly on Windows, but it has not been extensively tested. Please contact Abhiroop Ghosh (ghoshab1@msu.edu) for any issues faced.

<h2>Running unit tests to ensure everything is in order:</h4>
1. Clone the repo:

    ```https://github.com/abhiroopghosh71/solid-rocket-design-optimization.git```
    
2. Change into the working directory ```cd solid-rocket-design-optimization```

3. To install the dependencies use ```pip install requirements.txt```. If you are using Anaconda, it is recommended to create a new virtual environment using ```conda create --name <envname> --file requirements.txt```

4. Run ```make clean```

5. Run ```make```

6. Run ```python -m unittest tests/test_rocket_evaluator.py tests/test_pymoo_optimization.py```

<h2>A sample case for the standalone rocket burn simulator</h4>
From the ```solid-rocket-design-optimization``` directory, run ```python rocket_eval_demo.py```. The file has a predefined input in the ```__main__``` block.

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
The code uses NSGA-II [3] optimization algorithm to find the solid fuel configuration to match a given thrust profile.

1. Change into the solid-rocket-design-optimization directory.

    ```cd solid-rocket-design-optimization```

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

2. ```--target-thrust```: Thrust profile to be matched by the optimizer. Default is 'baseline'. Other thrust profiles 
may be used by giving the corresponding file name under ```problems/rocket_propellant_design/thrust_profiles/```.

4. ```--ngen```: Number of generations of NSGA-II.

5. ```--popsize```: Population size of NSGA-II.

For example,
```python optimize.py --ngen 100 --popsize 200 --target-thrust baseline --seed 1234```

This runs the optimization for 100 generations and a population size of 200 with a seed of 1234.

More details about the solid rocket burn simulator written in C is provided in the file 
```librocket/README_ROCKET.md```. Please report issues to me, Abhiroop Ghosh, at ghoshab1@msu.edu.

<h2>References:</h4>
1. A. Ghosh, E. Goodman, K. Deb, R. Averill, A. Diaz, "A Large-scale Bi-objective Optimization of Solid Rocket Motors Using Innovization," _2020 IEEE Congress on Evolutionary Computation (CEC)_, 1–8. https://doi.org/10.1109/CEC48606.2020.9185861

2. J. Blank and K. Deb, pymoo: Multi-Objective Optimization in Python, in IEEE Access, vol. 8, pp. 89497-89509, 2020, https://doi.org/10.1109/ACCESS.2020.2990567

3. K. Deb, A. Pratap, S. Agarwal and T. Meyarivan, "A fast and elitist multiobjective genetic algorithm: NSGA-II," in _IEEE Transactions on Evolutionary Computation_, vol. 6, no. 2, pp. 182-197, April 2002. https://doi.org/10.1109/4235.996017
