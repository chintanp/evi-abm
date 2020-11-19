/**
* Name: seedtestint
* Based on the internal empty template. 
* Author: Chintan Pathak
* Tags: 
*/


model seedtestint

/* Insert your model definition here */
global {
	float seed <- 123.0; 
	
	init {
		write(seed);
		float prob <- 0.512;
		
		write(bool(binomial(1, prob)));
	}
}


experiment nothing {
	
}