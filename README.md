# Parametercalibrationframework
A parameter calibration framework for urban flood modeling

This calibration framework consists of the CA-LLSO algorithm developed by Feng-Feng Wei (available at: https://github.com/CarrieWei/CA-LLSO_Code), and an urban flood model capable of using four different methods: BB, BH, BR, and BP. 

The CA-LLSO algorithm is written in Python, while the urban flood model is written in Fortran.

To run the urban flood model, users need to provide a condition file and the required input files contained in Input.zip, which should be prepared according to the specific case.
In the Input.zip, you can find a simple test case to verify that the program runs correctly.

Corresponding paper: Feng-Feng Wei, Wei-Neng Chen, Qiang Yang, Jeremiah Deng, Xiao-Nan Luo, Hu Jin and Jun Zhang, "A Classifier-Assisted Level-Based Learning Swarm Optimizer for Expensive Optimization", IEEE Transactions on Evolutionary Computation, vol.25, no.2, pp.219-233, 2021.
