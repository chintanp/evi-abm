/**
* Name: seedtestint
* Based on the internal empty template. 
* Author: Chintan Pathak
* Tags: 
*/


model seedtestext

/* Insert your model definition here */
global {
	
	file seedfile <- csv_file("../seed.csv", false);
	matrix seedm <- matrix(seedfile);
	float seed <- float(seedm[0, 0]);
	
	// float seed <- 123.0; 
	
	init {
		write(seed);
		float prob <- 0.512;
		
		write(bool(binomial(1, prob)));
	}
}


experiment nothing {
	
}