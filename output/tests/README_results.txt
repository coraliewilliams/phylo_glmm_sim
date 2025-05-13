

sim_results_test_1
	brms iter=2000 (default)
	MCMCglmm nitt=13000 (default)


sim_results_test_2
	brms iter=4000 (default x2)
	MCMCglmm nitt=53000 (default x5)


sim_results_test_3
	brms iter=4000 (default x2)
	MCMCglmm nitt=303000 (default x30)
	
	
~~~~~~~~~~~~
sim_results_test_3a
	MCMCglmm settings
	nitt=103000 x10
	thin=100 x10

sim_results_test_3b
	MCMCglmm settings
	nitt=503000 x50
	thin=200 x20
	
	
sim_results_test_3c
	MCMCglmm settings
	nitt=203000 x20
	thin=150 x15
	
sim_results_test_3d
	MCMCglmm settings
	nitt=403000 x40
	thin=10 (default)
~~~~~~~~~~~



sim_results_test_4	
	brms iter=10000 (default x5)
	MCMCglmm nitt=13000 (default x30)