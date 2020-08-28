/**
* Name: readtest
* Based on the internal empty template. 
* Author: chintan
* Tags: 
*/


model readtest

/* Insert your model definition here */

global skills: [SQLSKILL] {
	
	file envfile <- csv_file("../.env", true);
	matrix dbcredsm <- matrix(envfile);
	
	int analysis_id <- 75;
	string db_host <- dbcredsm[0, 0];
	string db_type <- dbcredsm[1, 0];
	string db_name <- dbcredsm[2, 0];
	string db_port <- dbcredsm[3, 0];
	string db_user <- dbcredsm[4, 0];
	string db_pwd <- dbcredsm[5, 0];
	
	map<string, string> DBPARAMS <- [ 'host'::db_host, 'dbtype'::db_type, 'database'::db_name, 'port'::db_port, 'user'::db_user, 'passwd'::db_pwd];
	
	list<list<list>> sim_start_time <- select(DBPARAMS, 'select sim_start_time from analysis_record where analysis_id = ' + analysis_id);
	
	date starting_date <- date(sim_start_time[2][0][0]);
	init {
		write (db_host + db_type + db_name + db_port + db_user + db_pwd);
		
		write(string(current_date));
		
		
	}
	
}

experiment test {
	
} 