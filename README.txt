To run main_IPOPT, first 

(1) Install IPOPT and the MATLAB toolbox ebertolazzi/mexIPOP​T
(2) Set the polynomial degree for voltage with state-of-charge (open_circuit_deg)
(3) Set the polynomial degree for internal resistance with state-of-charge (resistance_deg)
(4a) Set savecsv = true if you want to save the optimisation data to a csv file, otherwise set to false.
(4b) If savecsv = true, set tag = 'what you want your file to be named'.

To run main_MPC, first

(1) Install IPOPT and the MATLAB toolbox ebertolazzi/mexIPOP​T
(2) Ensure that the parameters U_max, R_max, SOC_max and d are the same in test_IPOPT_MPC
(3) Set the desired values for L and K
(4a) Set savecsv = true if you want to save the optimisation data to a csv file, otherwise set to false.
(4b) If savecsv = true, set tag = 'what you want your file to be named'.

