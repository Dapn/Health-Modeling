**************************************************************************************************
*************************Health Simulation - modeling for health decision*************************
***MicroSimulation of a hypothetical disease based on Markov model in R Programming (R Studio)**** 
********************Simulations Sick_Sicker_model_1 and Sick_Siker_model_2************************
Based on : 
Jalal H, et al. An Overview of R in Health Decision Sciences. Med. Dec. Mak. 2017; 37(3): 735-746.
Krijkamp EM, et al. Microsimulation modeling for health decision sciences using R: A tutorial. Med Decis Making. 2018;38(3):400-22.


Files:
Sick_Sicker_model_1 -  Micro Simulation Model with Monte Carlo Sampler without memory nor adaptive individual characteristics
Sick_Sicker_model_2 - model 1 updated with memory and ind features
TR_Dist_Graph - Population Progression in each Cycle Graph, writes into file TR_Population_Distribution_no_trt and TR_Population_Distribution_trt


Model 1 and 2 Assumtions and Characteristics:

- Comparing 2 Conditions : Treating and non-Treating
- Sick Sicker type modelation
- 4 Health States (H, S1, S2, D)
- In Treating, individuals in the state S1 and S2 receive treatment
- Only individuals in S1 will beneficiate from the treatment 
- All individual are healthy in the bigginig of simulation  (Cycle O)
- Results : Costs, QALYs,Incremental Costs, QALYs Gained, ICER. All with MCSE: Monte Carlo Standard Error    


For Model 2 - Memory and Advanced/Ind adjustment
- Each individual was a unique feature of baseline treating effect Between 0,95 - 1,05  (receptivity of a treatment)
- Increasing mortality rate every additional cycle/year in state S1 and S2
- Decreasing the utility(QALYs) of treated sick individuals with every additional year being in S1/S2 (note S2 does not beneficiate from treatment)


----------------------------------------------------------
Simulation Results based on random seed=1 and n=1000000 :

                   Costs  *  QALYs     * Incremental Costs  * QALYs Gained *   ICER
No Treatment       62728 38  15.29 0.005                                           
Treatment         117557 73  15.797 0.006             54830 37        0.507 0 108111
* are MCSE values 


 








