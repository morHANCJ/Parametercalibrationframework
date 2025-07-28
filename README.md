# Parametercalibrationframework
A parameter calibration framework for urban flood modeling

This calibration framework consists of the CA-LLSO algorithm developed by Feng-Feng Wei (available at: https://github.com/CarrieWei/CA-LLSO_Code), and an urban flood model capable of using four different methods: BB, BH, BR, and BP. 

The CA-LLSO algorithm is written in Python, while the urban flood model is written in Fortran.

To run the urban flood model, users need to provide a condition file and the required input files contained in Inputfile.zip, which should be prepared according to the specific case.

In the testing.zip, you can find a simple test case to verify that the program runs correctly.

Currently, the MainProgram.f90 file is written based on the BP method. To use other methods, you need to manually adjust the parameters. In particular, for the BB and BH methods, a Buildingmark parameter is provided in the file to set the boundary conditions.

Corresponding paper: Feng-Feng Wei, Wei-Neng Chen, Qiang Yang, Jeremiah Deng, Xiao-Nan Luo, Hu Jin and Jun Zhang, "A Classifier-Assisted Level-Based Learning Swarm Optimizer for Expensive Optimization", IEEE Transactions on Evolutionary Computation, vol.25, no.2, pp.219-233, 2021.
