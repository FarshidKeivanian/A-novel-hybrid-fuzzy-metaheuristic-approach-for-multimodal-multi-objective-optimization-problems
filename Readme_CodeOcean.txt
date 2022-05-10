As advised by the CodeOcean team, the code was a bit adjusted to reduce the computational time. To return to the experimental setup section of the paper that was conducted by High Performance Computing (HPC) at the University of Newcastle, please change the stopping criteria, MaxRun and MaxFEs, as the experimental setup in the paper.

In the paper, we used HPC, Intel® Xeon® Gold 6150 Processor, 24.75M Cache, 2.70 GHz. For 22 single-objective benchmarks, the stopping criterion was set as the following: MFEs = 20,000 for benchmarks #f1, f3, f4 and f6; MFEs = 50,000 for benchmarks #f2, f5, f7, f8, f9 and f10; MFEs = 50,000 for {#f11 to #f16}; MFEs = 200,000 for {#f17 to #f21}; and MFEs = 400,000 for #f22. For 25 multi-objective benchmarks, F1-F25, the stopping criterion was set as MFEs = 40,000.

Using a common hardware configuration, Dual Core i5 processor with clock frequency 2.6 GHz×2 and 8 GB RAM, the compute time for the main code 'RunAll_in_Less_ComputeTime.m' will be around 73 minutes. To get better results such as those were reported in the paper, the stopping criteria should be adjusted according to the section experimental setup in the paper.

The main file is "RunAll_Less_ComputeTime.m" which runs the codes, individually.

The functions f1, f2, ..., f22 are single-objective functions, and the functions F1, F2 ...., F25 are multi-objective functions.

In the paper, we used HPC, Intel® Xeon® Gold 6150 Processor, 24.75M Cache, 2.70 GHz. We conducted the experimental study with the following stopping criteria for the benchmarks: MFEs = 20,000 for benchmarks #f1, f3, f4 and f6; MFEs = 50,000 for benchmarks #f2, f5, f7, f8, f9 and f10; MFEs = 50,000 for {#f11 to #f16}; MFEs = 200,000 for {#f17 to #f21}; and MFEs = 400,000 for #f22. For 25 multi-objective benchmarks, the stopping criterion was set as MFEs = 40,000.

For a new work, you can model a real-worl problem in your field into a single or multi-objective optimisation problem, and then solve it by using the proposed fuzzy adaptive metaheuristic methods, FAEICA and MOFAEICA, respectively. Indeed, these methods can solve any optimisation and/or prediction problems such as optimum engineering design, precise health-related predictive framework (machine learning models), stock portfolio optimization models, and data analysis problems, and so on, from single-objective and multi-objective optimisation perspective.

The list of references for the benchmarks:

Zitzler, E., K. Deb and L. Thiele (2000). "Comparison of multiobjective evolutionary algorithms: Empirical results." Evolutionary computation 8(2): 173-195

Suganthan, P. N., N. Hansen, J. J. Liang, K. Deb, Y.-P. Chen, A. Auger and S. Tiwari (2005). "Problem definitions and evaluation criteria for the CEC 2005 special session on real-parameter optimization." KanGAL report 2005005: 2005.

Zhang, Q., A. Zhou, S. Zhao, P. N. Suganthan, W. Liu and S. Tiwari (2008). "Multiobjective optimization test instances for the CEC 2009 special session and competition." University of Essex, Colchester, UK and Nanyang technological University, Singapore, special session on performance assessment of multi-objective optimization algorithms, technical report 264.

Li, X., A. Engelbrecht and M. G. Epitropakis (2013). "Benchmark functions for CEC’2013 special session and competition on niching methods for multimodal function optimization." RMIT University, Evolutionary Computation and Machine Learning Group, Australia, Tech. Rep.