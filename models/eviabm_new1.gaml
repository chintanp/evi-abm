/***
* Name: network
* Author: chintan
* Description: This file runs the simulation using the finite state machine control scheme
* Tags: fsm, final 
***/
//
model network_fsm

/* Insert your model definition here */
global skills: [SQLSKILL] {

/* Global Variables */

// Batch parameters
	string simulation_date <- "2019-07-01";

	// Get analysis ID
	file aidfile <- csv_file("../analysis_id", false);
	matrix aidm <- matrix(aidfile);
	int analysis_id <- int(aidm[0, 0]);
	

	// Database credentials
	file envfile <- csv_file("../.env", true);
	matrix dbcredsm <- matrix(envfile);
	string db_host <- dbcredsm[0, 0];
	string db_type <- dbcredsm[1, 0];
	string db_name <- dbcredsm[2, 0];
	string db_port <- dbcredsm[3, 0];
	string db_user <- dbcredsm[4, 0];
	string db_pwd <- dbcredsm[5, 0];

	// Simulation level

	// DB Details 
	map<string, string> DBPARAMS <- ['host'::db_host, 'dbtype'::db_type, 'database'::db_name, 'port'::db_port, 'user'::db_user, 'passwd'::db_pwd];
	date cd <- date("now");
	int simulation_time <- 2 * 24 * 60 * 60; // time for which this simulation will run
	float step <- 1 #mn; // Time step of 1 minute for the simulation
	// start the simulation based on the sim_start_time from the trip generation 
	list<list<list>> sim_start_time <- select(DBPARAMS, 'select sim_start_time from analysis_record where analysis_id = ' + analysis_id);
	date starting_date <- date(sim_start_time[2][0][0]);
	float lookup_distance <- 10 / 0.000621371; // convert miles to m - this is the distance from the road that one looks for charging stations
	float reconsider_charging_time <- 10.0; // Time in minutes to reconsider charging decision
	// files
	file roads_shapefile <- file("../includes/SpatialJoin/WA_roads_4326.shp"); // WA state shapefile, with road speed limit, major roads only - transformed to EPSG:4326 in QGIS
	string start_time_file <- "../results/logs/simulation_time_" + simulation_date + "_" + string(cd, 'yyyy-MM-dd-HH-mm-ss') + ".csv";

	// DB queries 
	string
	bevse_query <- 'SELECT dcfc_count, ev_network, longitude, latitude, ev_connector_types, bevse_id, connector_code, dcfc_fixed_charging_price, dcfc_var_charging_price_unit, dcfc_var_charging_price, dcfc_fixed_parking_price, dcfc_var_parking_price_unit, dcfc_var_parking_price, zip FROM built_evse where dcfc_count >= 1;';
	string
	nevse_query <- 'SELECT nevse_id, longitude, latitude, connector_code, dcfc_plug_count, dcfc_fixed_charging_price, dcfc_var_charging_price_unit, dcfc_var_charging_price, dcfc_fixed_parking_price, dcfc_var_parking_price_unit, dcfc_var_parking_price from new_evses where dcfc_plug_count > 0 and analysis_id = ' + analysis_id;
	string evtrip_query <- 'select e.veh_id, e.origin_zip, z1.latitude as olat, z1.longitude as olng, e.destination_zip, 
		z2.latitude as dlat, z2.longitude as dlng, e.soc, e.trip_start_time, b.range_fe, 
		b.capacity, b.fuel_consumption, b.connector_code
	from evtrip_scenarios e	
	inner join zipcode_record z1 on cast(e.origin_zip as text) = z1.zip
	inner join zipcode_record z2 on cast(e.destination_zip as text) = z2.zip
	inner join wa_bevs b on e.veh_id = b.veh_id
	where e.analysis_id =' + analysis_id + " and origin_zip = '98177' and destination_zip = '98671' order by e.veh_id::int";
	string param_query <- "select ap.param_id, param_name, ap.param_value from analysis_params ap
							join sim_params sp on sp.param_id = ap.param_id
							where ap.analysis_id = " + analysis_id + " and sp.param_type IN ('global', 'eviabm') 
							order by ap.param_id asc;";

// matrices 

	// lists
	list<road> roadsList;
	list<charging_station> all_chargers;
	list<charging_station> all_chargers_chademo;
	list<charging_station> all_chargers_combo;
	list<list<list>> params <- nil;
path plot_shortest_path;

	// spatial
	geometry shape <- envelope(roads_shapefile);
	graph road_network;
	graph road_network_weighted;
	graph road_network_driving;

	// Other variables
	int SOC_MAX <- 80; // SOC at which to stop charging
	int SOC_MIN <- 0; // SOC at which to stop driving if no charging is available
	int MIN_SOC_CHARGING <- 60; // Minimum SOC to consider charging - above this SOC the user must never charge
	int MAX_SOC_CHARGING <- 20; // maximum SOC to consider charging - below this SOC, the user must charge
	int connector_code;
	float BLOCK_SIZE <- 200.0; // Size of block in meters, where relocation preferred over waiting
	
		list<charging_station> plot_chargers_nearby;
	list<charging_station> plot_compat_chargers;
	list<point> plot_cpts;
	geometry my_circle;
	geometry plot_my_circle;
	road my_path;
	
	
	
	init {
		write (seed);
		save string(date("now")) type: csv header: false to: start_time_file rewrite: false;
		write ("Initiating new simulation");
		// Read in the params from the DB
		params <- list<list<list>>(select(DBPARAMS, param_query));
		// Assign the parameters as read from the DB
		// seed <- float(params[2][0][2]);
		lookup_distance <- float(params[2][2][2]) / 0.000621371;
		step <- float(params[2][3][2]) #mn;
		reconsider_charging_time <- float(params[2][4][2]);
		SOC_MAX <- int(params[2][5][2]);
		SOC_MIN <- int(params[2][6][2]);
		MIN_SOC_CHARGING <- int(params[2][7][2]);
		MAX_SOC_CHARGING <- int(params[2][8][2]);
		BLOCK_SIZE <- float(params[2][9][2]);
		write ("Analysis_id: " + analysis_id);
		write ("Seed: " + seed);
		write ("lookup_distance: " + lookup_distance);
		write ("step: " + step);
		write ("reconsider_charging_time: " + reconsider_charging_time);
		write ("SOC_MAX: " + SOC_MAX);
		write ("SOC_MIN: " + SOC_MIN);
		write ("MIN_SOC_CHARGING: " + MIN_SOC_CHARGING);
		write ("MAX_SOC_CHARGING: " + MAX_SOC_CHARGING);
		write ("BLOCK_SIZE: " + BLOCK_SIZE);

		// Create roads from the shapefile 
		// Reading attributes 'Spd' and 'ID', to add more attributes, go to R, or QGIS		
		create road from: roads_shapefile with: [maxspeed::float(read('spd')) #miles / #hour, road_ID::int(read('id'))];
		write ("Roads created");
		// create a list of roads - this is needed for finding the current road etc.
		roadsList <- (road as list);

		// Road network as a graph from road agents
		road_network <- as_edge_graph(road);

		// Add time to travel as weights to the map
		map<road, float> map_weights <- road as_map (each::each.shape.perimeter / each.maxspeed);
		road_network_weighted <- copy(road_network) with_weights map_weights;
		road_network_driving <- as_edge_graph(road);

		// Create charging station from built_evse
		list<list<list>> bevses <- list<list<list>>(select(DBPARAMS, bevse_query));
		loop ii from: 0 to: length(bevses[2]) - 1 {
			create charging_station with:
			[shape::to_GAMA_CRS({float(bevses[2][ii][2]), float(bevses[2][ii][3])}, 'EPSG:4326'), station_id::string(bevses[2][ii][5]), zip::int(bevses[2][ii][13]), dcfc_plug_count::int(bevses[2][ii][0]), ev_network::bevses[2][ii][1], ev_connector_types::bevses[2][ii][4], ev_connector_code::int(bevses[2][ii][6]), dcfc_fixed_charging_price::float(bevses[2][ii][7]), dcfc_var_charging_price_unit::string(bevses[2][ii][8]), dcfc_var_charging_price::float(bevses[2][ii][9]), dcfc_fixed_parking_price::float(bevses[2][ii][10]), dcfc_var_parking_price_unit::string(bevses[2][ii][11]), dcfc_var_parking_price::float(bevses[2][ii][12])];
		}

		write ("Built charging stations created");
		// Add a prefix of b to built_evse
		int built_evse_count <- length(charging_station);
		loop i from: 0 to: built_evse_count - 1 {
			charging_station[i].station_id <- "b" + charging_station[i].station_id;
		}

		list<list<list>> nevses <- list<list<list>>(select(DBPARAMS, nevse_query));
		// Create charging station from new_evses
		if (length(nevses[2]) > 0) {
			loop ii from: 0 to: length(nevses[2]) - 1 {
				create charging_station with:
				[shape::to_GAMA_CRS({float(nevses[2][ii][1]), float(nevses[2][ii][2])}, 'EPSG:4326'), station_id::string(nevses[2][ii][0]), dcfc_plug_count::int(nevses[2][ii][4]), ev_connector_code::int(nevses[2][ii][3]), dcfc_fixed_charging_price::float(nevses[2][ii][5]), dcfc_var_charging_price_unit::string(nevses[2][ii][6]), dcfc_var_charging_price::float(nevses[2][ii][7]), dcfc_fixed_parking_price::float(nevses[2][ii][8]), dcfc_var_parking_price_unit::string(nevses[2][ii][9]), dcfc_var_parking_price::float(nevses[2][ii][10])];
			}

			// Add a prefix of n to new_evse 
			loop i from: built_evse_count to: length(charging_station) - 1 {
				charging_station[i].station_id <- "n" + charging_station[i].station_id;
			}

		}

		write ("New charging stations created");
		// Create aggregations for charging stations		
		all_chargers <- (charging_station as list);
		all_chargers_chademo <- all_chargers where (each.ev_connector_code = 1 or each.ev_connector_code = 3);
		all_chargers_combo <- all_chargers where (each.ev_connector_code = 2 or each.ev_connector_code = 3);

		// Create EV agents   
		list<list<list>> evs <- list<list<list>>(select(DBPARAMS, evtrip_query));
//		loop ii from: 0 to: length(evs[2]) - 1 {
//			create EVs with:
//			[shape::to_GAMA_CRS({float(evs[2][ii][3]), float(evs[2][ii][2])}, "EPSG:4326"), veh_ID::string(evs[2][ii][0]), range::float(evs[2][ii][9]), capacity::float(evs[2][ii][10]), fuel_consumption::float(evs[2][ii][11]), connector_code::int(evs[2][ii][12]), the_target::point(to_GAMA_CRS({float(evs[2][ii][6]), float(evs[2][ii][5])}, "EPSG:4326")), origin_zip::int(evs[2][ii][1]), destination_zip::int(evs[2][ii][4]), SOC::float(evs[2][ii][7]), trip_start_time::date(string(evs[2][ii][8]), '%Y-%M-%D %h:%m:%s')];
//		}

		loop ii from: 0 to: 0 {
			create EVs with: [shape::road_network_driving.vertices closest_to
			to_GAMA_CRS({float(evs[2][ii][3]), float(evs[2][ii][2])}, "EPSG:4326"), veh_ID::string(evs[2][ii][0]), range::float(evs[2][ii][9]), capacity::float(evs[2][ii][10]), fuel_consumption::float(evs[2][ii][11]), connector_code::int(evs[2][ii][12]), the_target::road_network_driving.vertices
			closest_to
			point(to_GAMA_CRS({float(evs[2][ii][6]), float(evs[2][ii][5])}, "EPSG:4326")), origin_zip::int(evs[2][ii][1]), destination_zip::int(evs[2][ii][4]), SOC::float(evs[2][ii][7]), trip_start_time::date(string(evs[2][ii][8]), '%Y-%M-%D %h:%m:%s')];
		}

		write ("EVs created");
		string ev_ids <- EVs[0].veh_ID;
		if (length(EVs) > 1) {
		}

	}

//	reflex log_power_use {
//		string pd_query <- 'INSERT INTO evse_power_draw (simulation_ts, evse_id, analysis_id, power_val) VALUES';
//		string valq <- "('" + string(current_date) + "', '" + all_chargers[0].station_id + "', " + analysis_id + ", " + all_chargers[0].current_power_draw + ")";
//		loop ii from: 1 to: length(all_chargers) - 1 {
//			valq <- valq + " , " + "('" + string(current_date) + "', '" + all_chargers[ii].station_id + "', " + analysis_id + ", " + all_chargers[ii].current_power_draw + ")";
//		}
//
//		pd_query <- pd_query + valq;
//		do executeUpdate params: DBPARAMS updateComm: pd_query;
//	}
//
//	reflex log_EV_agents {
//		if (length(EVs) > 0) {
//			string
//			info_query <- 'INSERT INTO ev_info (simulation_ts, analysis_id, veh_id, lat_val, lng_val, soc_val, prob_val, speed_val, state_val, tocharge_val, chargers_nearby, nearest_evse_id, nearest_evses, charging_decision_time) VALUES';
//			point l_4326 <- point(CRS_transform(EVs[0].location, "EPSG:4326"));
//			float latitude <- l_4326.y;
//			float longitude <- l_4326.x;
//			string
//			valq_info <- "('" + string(current_date) + "'," + analysis_id + ", '" + EVs[0].veh_ID + "', " + latitude + ", " + longitude + ", " + with_precision(EVs[0].SOC, 3) + ", " + with_precision(EVs[0].prob_charging, 3) + ", " + with_precision(EVs[0].veh_speed, 1) + ", '" + EVs[0].state + "', " + EVs[0].to_charge + ", '" + EVs[0].chargers_onpath + "', '" + EVs[0].nearest_evse + "', '" + EVs[0].nearest_evses + "', '" + string(EVs[0].charging_decision_time) + "')";
//			if (length(EVs) > 1) {
//				loop jj from: 1 to: length(EVs) - 1 {
//					l_4326 <- point(CRS_transform(EVs[jj].location, "EPSG:4326"));
//					latitude <- l_4326.y;
//					longitude <- l_4326.x;
//					valq_info <-
//					valq_info + ", " + "('" + string(current_date) + "'," + analysis_id + ", '" + EVs[jj].veh_ID + "', " + latitude + ", " + longitude + ", " + with_precision(EVs[jj].SOC, 3) + ", " + with_precision(EVs[jj].prob_charging, 3) + ", " + with_precision(EVs[jj].veh_speed, 1) + ", '" + EVs[jj].state + "', " + EVs[jj].to_charge + ", '" + EVs[jj].chargers_onpath + "', '" + EVs[jj].nearest_evse + "', '" + EVs[jj].nearest_evses + "', '" + string(EVs[jj].charging_decision_time) + "')";
//				}
//
//			}
//
//			info_query <- info_query + valq_info;
//			do executeUpdate params: DBPARAMS updateComm: info_query;
//		}
//
//	}
	///////////////////////////////////////////
	// Let the simulation go on untill all the agents die
	///////////////////////////////////////////

	// stop the simulation since time is up
	reflex halting_timeup when: cycle * step = simulation_time {
		string evse_util_query <- 'INSERT INTO evse_util (analysis_id, evse_id, util_val) VALUES';
		string valq_evse_util <- "(" + analysis_id + ", '" + charging_station[0].station_id + "', " + charging_station[0].energy_consumed + ")";
		loop ii from: 1 to: length(charging_station) - 1 {
			valq_evse_util <- valq_evse_util + ", " + "(" + analysis_id + ", '" + charging_station[ii].station_id + "', " + charging_station[ii].energy_consumed + ")";
		}

		evse_util_query <- evse_util_query + valq_evse_util;
		// do executeUpdate params: DBPARAMS updateComm: evse_util_query;
		save string(date("now")) type: csv header: false to: start_time_file rewrite: false;
		string status_upd_query <- "update analysis_record set status = 'solved' where analysis_id = " + analysis_id;
		// do executeUpdate params: DBPARAMS updateComm: status_upd_query;
		write ("Time is up");
		do die;
	}

	// stop the simulation since all EVs are done
	reflex halting_alldone when: empty(EVs) {
		string evse_util_query <- 'INSERT INTO evse_util (analysis_id, evse_id, util_val) VALUES';
		string valq_evse_util <- "(" + analysis_id + ", '" + charging_station[0].station_id + "', " + charging_station[0].energy_consumed + ")";
		loop ii from: 1 to: length(charging_station) - 1 {
			valq_evse_util <- valq_evse_util + ", " + "(" + analysis_id + ", '" + charging_station[ii].station_id + "', " + charging_station[ii].energy_consumed + ")";
		}

		evse_util_query <- evse_util_query + valq_evse_util;
		// do executeUpdate params: DBPARAMS updateComm: evse_util_query;
		save string(date("now")) type: csv header: false to: start_time_file rewrite: false;
		string status_upd_query <- "update analysis_record set status = 'solved' where analysis_id = " + analysis_id;
		// do executeUpdate params: DBPARAMS updateComm: status_upd_query;
		write ("all EVs are done");
		do die;
	}

}

species road schedules: [] {
	geometry display_shape <- line(shape.points, 2.0);
	float maxspeed;
	int road_ID;
	rgb road_color <- #red;

	aspect geom {
		draw shape color: road_color; // (maxspeed >= 60 #miles / #hour) charging_decision_time? #yellow : (((maxspeed >= 50 #miles / #hour) ? road_color : ((maxspeed >= 30 #miles / #hour) ? #green : #blue))); 
	}

	aspect geom3D {
		draw display_shape color: #black;
	}

}

species EVs skills: [moving, SQLSKILL] control: fsm {
	string veh_ID; // identifying number, don't know the utility of this yet
	float capacity; // battery capacity in kWhr
	float range; // EV range in miles
	float fuel_consumption; // Fuel i.e. energy consumption in kWhr per 100 mile  
	float SOC; // the current State of Charge of the vehicle - expressed as a % 
	int connector_code;
	point the_target;
	road currentRoad <- nil;
	float veh_speed <- 5 #miles / #hour;
	int origin_zip;
	int destination_zip;
	float speed <- 1 #miles / #hour; // speed of motion
	float distance_travelled <- 0.0;
	float dist;
	float dist_in_miles;
	date trip_start_time; // Minutes since 0000 hrs when the trip starts
	float trip_distance;
	list<charging_station> chargers_onpath <- [];
	list<charging_station> compat_chargers;
	point charger_target;
	float remaining_range <- range * SOC / 100; // vehicles do not always start with the full SOC
	float starting_SOC <- -1.0;
	string charge_start_time <- nil;
	path shortest_path;
	bool to_charge <- false;
	date charging_decision_time <- date(1, 1, 1);
	list<point> pts_on_path;
	point cpt_on_path;
	float prob_charging <- 0.0;
	point prev_loc <- location;
	map<charging_station, list> charger_dists <- [];
	map<charging_station, list> charger_dists_old <- [];
	// float dist_dest <- 0.0;
	road my_path;
	float dist_to_road;
	geometry my_circle;
	//	cs_pt_on_path csp_nearest <- nil;
	//	cs_pt_on_path csp_next_nearest <- nil;
	//	list<cs_pt_on_path> csp_others <- nil;
	//	list<cs_pt_on_path> csp_nearby <- nil;
	list<point> cpts_on_path;
	list<float> cs_dists;
	list<charging_station> chargers_nearby <- [];
	charging_station next_nearest_evse;
	
	charging_station nearest_evse <- nil;
	point snap_point;

	init {

	// The list of compatible chargers based on teh connector code of the EV
	// 1. Chademo 2. Combo 4. Tesla		
		if (connector_code = 1) {
			compat_chargers <- all_chargers_chademo;
		} else if (connector_code = 2) {
			compat_chargers <- all_chargers_combo;
		} else if (connector_code = 4) {
			compat_chargers <- all_chargers;
		}
plot_compat_chargers <- compat_chargers;
		// Find the shortest path on the graph between the EV and the target, 
		// this should happen only once when the trip starts and then the vehicle will traverse on this path	
		shortest_path <- path_between(road_network_driving, self, the_target);
		plot_shortest_path <- shortest_path;
		// Get a list of points on path
		pts_on_path <- geometry(shortest_path.segments) points_on (10);
		location <- pts_on_path closest_to (self);
		trip_distance <- shortest_path.shape.perimeter * 0.000621371;

		// Find the chargers that are within the 'lookup_distance' of the path. 
		// This isnt the 'on road' distance to the charging station, but just a buffer around the path 
		chargers_onpath <- (compat_chargers overlapping (shortest_path.shape + lookup_distance));
		plot_chargers_nearby <- chargers_onpath;
		if (empty(chargers_onpath)) {
		} else {
			loop cs over: chargers_onpath {
				snap_point <- pts_on_path closest_to (cs);
				add snap_point to: cpts_on_path;
			}

			using topology(road_network_driving) {
				loop cp over: cpts_on_path {
					add with_precision(distance_to(self, cp), 1) to: cs_dists;
				}

			}

		} }

		// Action (read function in traditional programming langauges) to consume fuel 
	// Start with avg fuel consumption
	// replace with speed, elevation dependent
	action update_states {
		float time_elapsed <- step;
		float fuel_consumed <- 0.0;
		if (location = the_target) {
		} else {

		// as long as as the vehicle is approaching the target, the SOC should be updated
			dist <- veh_speed * time_elapsed; // changed from veh_speed to check the difference
			dist_in_miles <- dist * 0.000621371;
			fuel_consumed <- fuel_consumption * dist_in_miles / 100; // Fuel consumed per mile = fuel_consumption / 100
			SOC <- SOC - fuel_consumed * 100 / capacity;
			// trip_length_remaining <- trip_length_remaining - dist; 
			remaining_range <- remaining_range - dist_in_miles;
			distance_travelled <- distance_travelled + dist_in_miles;
			// dist_dest <- (trip_distance - distance_travelled) / 0.000621371; // this is dist to dest in m
		}

	}

	action update_chargers_nearby {

	// Variables to keep track of the old state
		list<float> cs_dists_old <- copy(cs_dists);
		list<point> cpts_on_path_old <- copy(cpts_on_path);
		list<charging_station> chargers_onpath_old <- copy(chargers_onpath);
		list<float> cs_dists_new;
		int cs_near_old <- length(chargers_onpath_old);
		if (cs_near_old > 0) {
		// Get the distances from the EV current location (self) to each of the chargers
			using topology(road_network_driving) {
				loop cn from: 0 to: cs_near_old - 1 {
				// write("Within");
					float dist_cn <- with_precision(distance_to(self, cpts_on_path_old[cn]), 1);
					// write(dist_cn);
					// If the distance to a charger is increasing, then we are past it and so it should not be 
					// part of our chargers_nearby list and all associated data-structures need to updated
					// ************* To test the sensitivity of this ************************************
					if (dist_cn > cs_dists_old[cn]) {
					// write("remove");
					// write(cs_dists_old[cn]);
						remove chargers_onpath_old[cn] from: chargers_onpath;
						remove cpts_on_path_old[cn] from: cpts_on_path;
					}

				}
				// Find the new distances from the self to the charging stations ahead
				// maybe can be combined with the previous loop
				if (length(chargers_onpath) > 0) {
					loop cn from: 0 to: length(chargers_onpath) - 1 {
						add with_precision(distance_between(topology(road_network_driving), [self, cpts_on_path[cn]]), 1) to: cs_dists_new;
					}

					cs_dists <- copy(cs_dists_new);
				}

			}

		}

		do search_charger;
	}

	action search_charger {
		if (length(chargers_onpath) = 0) {
			nearest_evse <- nil;
		} else {
		// using topology(road_network_driving) {
		// Find the charger closest to the current location
			nearest_evse <- chargers_onpath closest_to location;
			// nearest_evses <- chargers_onpath closest_to (location, 2);
			ask nearest_evse {
				myself.chargers_nearby <- myself.chargers_onpath at_distance (2 * myself.dist);
			}
			
			next_nearest_evse <- (chargers_onpath - chargers_nearby - nearest_evse) closest_to location;
			// }
			// write (nearest_evses);
		}

	}

	// This updates the EV and EVSE attrbutes during charging
	action charge {
	// if (SOC < 100) {
		nearest_evse.energy_consumed <- nearest_evse.energy_consumed + (nearest_evse.max_power * step / 60.0 / 60.0); // Energy in kWhr
		SOC <- SOC + (nearest_evse.max_power * step * 100.0 / capacity / 60.0 / 60.0);
		remaining_range <- remaining_range + (nearest_evse.max_power * step * 100.0 / 60.0 / 60.0 / fuel_consumption);
		// }

	}

	////////////////////////////////////////////
	// Finite state machine control architecture
	////////////////////////////////////////////

	////////////////////////////////////////////
	// 1. Initial state 
	// Resting - to denote the time spent at home before starting the trip 
	////////////////////////////////////////////
	state resting initial: true {
		transition to: driving when: current_date >= trip_start_time;
	}

	///////////////////////////////////////////
	// 2. State - 2 
	// Driving - to denote the driving on the road.
	////////////////////////////////////////////
	state driving {
		int temp_ID;
		bool must_charge_now <- false;
		bool dont_charge_now <- false;
		// Find the current road on which the vehicle is traveling
		currentRoad <- (roadsList select (each != currentRoad)) with_min_of (each distance_to self);
		// Find the speed of traffic on this road
		ask currentRoad {
			myself.veh_speed <- self.maxspeed; // Use maxpseed until the speed of traffic is found
		}
		// Goto the destination 
		path path_travelled <- goto(target: the_target, on: road_network_driving, speed: veh_speed, return_path: true);
		do update_states; // This is where SOC gets updated
		do update_chargers_nearby; // This updates the list of chargers that are still relevant, as all those behind in the route will be removed
		if (nearest_evse != nil) {
			if (min(cs_dists) <= 2 * dist) {
				using topology(road_network_driving) {
					if (next_nearest_evse != nil) {
					//					write ("In if");code 
						// list<charging_station> other_evse <- nearest_evses - [nearest_evse];
						float dist_next_charger <- with_precision(distance_to(self, next_nearest_evse), 1);
						float dist_dest <- with_precision(distance_to(self, the_target), 1);
						if (remaining_range * 0.95 / 0.000621371 < dist_dest) {
							if (remaining_range * 0.95 / 0.000621371 < dist_next_charger) {
								must_charge_now <- true;
							}

						} else {
							dont_charge_now <- true;
						}

					} else {
						float dist_dest <- with_precision(distance_to(self, the_target), 1);
						if (0.95 * remaining_range / 0.000621371 < dist_dest) {
							must_charge_now <- true;
						} else {
							dont_charge_now <- true;
						}

					}

				}

				if (dont_charge_now = false) {
				// always charge if you cant make it to the destination or the next station
					if (must_charge_now = true) {
						to_charge <- true;
						charging_decision_time <- current_date;
					} else if (SOC <= MIN_SOC_CHARGING and charging_decision_time = date(1, 1, 1)) { //and charging_decision_time = date(1, 1, 1)
					// Further only think about charging if SOC <= MIN_SOC_CHARGING - 
					// this is the deterministic SOC-based charging choice to avoid high SOC charging
						to_charge <- charge_makes_sense();
						charging_decision_time <- current_date;
					} else if (SOC <= MIN_SOC_CHARGING and (current_date - charging_decision_time) / 60 >= reconsider_charging_time) {
						to_charge <- charge_makes_sense();
						charging_decision_time <- current_date;
					} else if (SOC <= MAX_SOC_CHARGING) {
						to_charge <- true;
						charging_decision_time <- current_date;
					} } } }

					// It is time to rest, if we have reached our destination
		transition to: finished when: location = the_target;

		// It is time to charge if charge 'makes sense'
		transition to: drive_to_charger when: (to_charge = true and nearest_evse != nil) {
			charger_target <- point(pts_on_path closest_to nearest_evse.location);
			charging_decision_time <- date(1, 1, 1); 
			must_charge_now <- false; // reset to origin, as otherwise this causes infinite loop at chargers
		}

		transition to: stranded when: SOC <= SOC_MIN;
	}

	///////////////////////////////////////////
	// 3. State 3
	// Drive_To_Charger - the vehicle is driving to a charging station
	////////////////////////////////////////////
	state drive_to_charger {
		int temp_ID;
		// Find the current road on which the vehicle is traveling
		currentRoad <- (roadsList select (each != currentRoad)) with_min_of (each distance_to self);
		// Find the speed of traffic on this road
		ask currentRoad {
			myself.veh_speed <- self.maxspeed; // Use maxpseed until the speed of traffic is found
			myself.speed <- self.maxspeed; // Use maxpseed until the speed of traffic is found

		}
		// Goto the charger 
		path path_travelled <- goto(target: charger_target, on: road_network_driving, speed: veh_speed, return_path: true);
		do update_states; // This is where SOC gets updated
		transition to: stranded when: SOC <= SOC_MIN;
		// It is time to charge, if we have reached a charger
		transition to: queue_for_charging when: location = charger_target;
	}

	///////////////////////////////////////////
	// 4. State 4
	// queue_for_charging - Create a queue for charging
	////////////////////////////////////////////
	state queue_for_charging {
		bool wait_for_charger <- false;
		bool ok_to_charge <- false;
		bool drive_on <- false;
		bool goto_wait <- false;
		bool relocate_nearby <- false;
		if (nearest_evse.plugs_in_use = nearest_evse.dcfc_plug_count) { // if the charging station is occupied
			if (distance_between(topology(road_network_driving), [self, the_target]) < remaining_range / 0.000621371) {
				drive_on <- true; // keep on driving if you can make it to the destination
			} else {
				if (next_nearest_evse != nil) {
				//					write ("In if");
					// list<charging_station> other_evse <- nearest_evses - [nearest_evse];
					// write(other_evse);
					float dist_next_charger <- with_precision(distance_between(topology(road_network_driving), [self, next_nearest_evse]), 1);
					if (dist_next_charger > BLOCK_SIZE or dist_next_charger > remaining_range / 0.000621371) {
						goto_wait <- true;
					} else {
						nearest_evse <- next_nearest_evse;
						relocate_nearby <- true;
					}

				} else {
					goto_wait <- true;
				}

			}

		} else if (nearest_evse.plugs_in_use < nearest_evse.dcfc_plug_count) {
			if (length(nearest_evse.waiting_evs) = 0) {
				ok_to_charge <- true;
			} else if (length(nearest_evse.waiting_evs) > 0) {
				goto_wait <- true;
			}

		}

		transition to: waiting when: goto_wait = true {
			add self to: nearest_evse.waiting_evs;
		}

		transition to: charging when: ok_to_charge = true;
		transition to: driving when: drive_on = true;
		transition to: drive_to_charger when: relocate_nearby = true {
			charger_target <- point(pts_on_path closest_to nearest_evse.location);
		}
	}

	///////////////////////////////////////////
	// 5. State 5
	// Charging - the vehicle is charging at a charging station 
	////////////////////////////////////////////
	state charging {
		enter {
			bool ok_to_charge;
			if (nearest_evse.plugs_in_use < nearest_evse.dcfc_plug_count) {

			// when this EV comes from waiting state
				if (nearest_evse.waiting_evs contains self) {
					ok_to_charge <- true;
					remove self from: nearest_evse.waiting_evs;
					starting_SOC <- SOC;
					charge_start_time <- string(current_date);
					// make this plug in use
					nearest_evse.plugs_in_use <- nearest_evse.plugs_in_use + 1;
				} else if (length(nearest_evse.waiting_evs) > 0) { // if from some other state, but there are other EVs waiting
					ok_to_charge <- false;
				} else { // if there are no EVs waiting
					starting_SOC <- SOC;
					charge_start_time <- string(current_date);
					// make this plug in use
					nearest_evse.plugs_in_use <- nearest_evse.plugs_in_use + 1;
					ok_to_charge <- true;
				}

			} else {
				ok_to_charge <- false;
			}

		}

		if (ok_to_charge = true) {

		// Perform charge
			do charge;
		}

		// If charged enough, go to driving
		transition to: driving when: SOC >= SOC_MAX {
			string charge_end_time <- string(current_date);
			float ending_SOC <- SOC;
			string evse_id <- nearest_evse.station_id;
			string cs_query <- 'INSERT INTO evse_charging_session (charge_start_time, charge_end_time, veh_id, starting_soc, ending_soc, evse_id, analysis_id) VALUES';
			string
			valq_cs <- "('" + charge_start_time + "', '" + charge_end_time + "', '" + veh_ID + "', " + starting_SOC + ", " + ending_SOC + ", '" + evse_id + "', " + analysis_id + ")";
			cs_query <- cs_query + valq_cs;
			// do executeUpdate params: DBPARAMS updateComm: cs_query;
			starting_SOC <- -1.0;
			charge_start_time <- nil;
			nearest_evse.plugs_in_use <- nearest_evse.plugs_in_use - 1;
			to_charge <- false;
			// remove nearest_evse from: chargers_onpath;
		}

		transition to: waiting when: ok_to_charge = false {
			if (!(nearest_evse.waiting_evs contains self)) {
				add self to: nearest_evse.waiting_evs;
			}

		}

	}

	/////////////////////////////////////////
	// 6. State 6
	// Finished - the vehicle is resting at the destination 
	////////////////////////////////////////////
	state finished {
		string fin_query <- 'INSERT INTO ev_finished (fin_ts, veh_id, origin_zip, destination_zip, analysis_id, trip_distance, distance_travelled) VALUES';
		string
		valq_fin <- "('" + string(current_date) + "', '" + veh_ID + "'," + origin_zip + "," + destination_zip + "," + analysis_id + ", " + with_precision(trip_distance, 0) + ", " + with_precision(distance_travelled, 0) + ")";
		fin_query <- fin_query + valq_fin;
		// do executeUpdate params: DBPARAMS updateComm: fin_query;
		transition to: done;
	}

	/////////////////////////////////////////
	// 7. State 7
	// stranded - the vehicle is not able to complete the trip and is out of charge enroute
	////////////////////////////////////////////
	state stranded {
		point l_4326 <- point(CRS_transform(location, "EPSG:4326"));
		float latitude <- l_4326.y;
		float longitude <- l_4326.x;
		string stranded_query <- 'INSERT INTO ev_stranded (stranded_ts, veh_id, stranded_lat, stranded_lng, origin_zip, destination_zip, analysis_id) VALUES';
		string
		valq_stranded <- "('" + string(current_date) + "', '" + veh_ID + "', " + latitude + ", " + longitude + ", " + origin_zip + ", " + destination_zip + ", " + analysis_id + ")";
		stranded_query <- stranded_query + valq_stranded;
		// do executeUpdate params: DBPARAMS updateComm: stranded_query;
		transition to: done;
		// do die;
	}

	/////////////////////////////////////////
	// 8. State 7
	// waiting - the vehicle is waiting to charge as the charging station is fully occupied
	////////////////////////////////////////////
	state waiting {
		bool ok_to_charge <- false;
		enter {
			string wait_start_time <- string(current_date);
			string ev_waiting_query <- 'INSERT INTO evse_evs_waiting (wait_start_time, veh_id, soc_val, evse_id, analysis_id) VALUES';
			string valq_ev_waiting <- "('" + wait_start_time + "', '" + veh_ID + "', " + SOC + ", '" + nearest_evse.station_id + "', " + analysis_id + ");";
			ev_waiting_query <- ev_waiting_query + valq_ev_waiting;
			// do executeUpdate params: DBPARAMS updateComm: ev_waiting_query;
		}

		if (nearest_evse.plugs_in_use < nearest_evse.dcfc_plug_count and first(nearest_evse.waiting_evs) = self) {
			ok_to_charge <- true;
			string ev_waiting_query <- 'UPDATE evse_evs_waiting set wait_end_time = ';
			string
			valq_ev_waiting <- "'" + string(current_date) + "' where veh_id = '" + veh_ID + "' and analysis_id = " + analysis_id + " and evse_id = '" + nearest_evse.station_id + "'; ";
			ev_waiting_query <- ev_waiting_query + valq_ev_waiting;
			// do executeUpdate params: DBPARAMS updateComm: ev_waiting_query;
		}

		transition to: charging when: ok_to_charge = true;
	}

	state done {
		do die; // kill the agents after they are done
	}
	// This is where we integrate Yan's CCDM 
	// charge will make sense depending on following factors: 
	// SOC, time_in_car, charging_cost,  charging_time, access_time, amenity_level (?), deviation(?) 
	bool charge_makes_sense {
	// Considering SDCM4 cofficients 
		float intercept <- 2.034;
		float c_soc <- -4.584;
		float c_dev <- 2.440;
		float c_time_in_car <- -0.069;
		float c_charging_cost <- -0.010;
		float c_charging_time <- -0.242;
		float c_access_time <- -0.025;
		float c_amenity_restroom <- 0.049;
		float c_amenity_more <- 0.213;
		float soc <- SOC / 100; // SOC in percent
		int dev <- 0; // still to find out how to calculate dev for a trip, but either 1 or 0
		float time_in_car <- (current_date - trip_start_time) / 3600; // time in hours since driving - need to check how this variable is affected after charging once 
		float charging_cost_in_dollar; // this is the charging cost in dollar - will be depedent on the price model of each EVSE
		float charging_price <- 0.0;
		float parking_price <- 0.0;
		float charging_time;
		// Still need to find values of following covariates
		float access_time;
		path path_to_cs;
		int amenity_restroom <- 1; // al charging stations have restrooms, ** this assumption needs validation **
		int amenity_more;
		
		do search_charger;
		
		
		if (nearest_evse != nil) {
			
			// Talk to the nearest EVSE to find out VSE specific parameters
			ask nearest_evse {
				charging_time <- (80 - myself.SOC) * myself.capacity / 100 / self.max_power; // energy used / power = time
				if (self.dcfc_var_parking_price_unit = "min") {
					parking_price <- self.dcfc_fixed_parking_price + charging_time * 60 * self.dcfc_var_parking_price;
				}

				if (self.dcfc_var_charging_price_unit = "min") {
					charging_price <- self.dcfc_fixed_charging_price + charging_time * 60 * self.dcfc_var_charging_price;
				} else if (self.dcfc_var_charging_price_unit = "kWh") {
					charging_price <- self.dcfc_fixed_charging_price + ((80 - myself.SOC) * myself.capacity * self.dcfc_var_charging_price / self.max_power);
				}

				charging_cost_in_dollar <- parking_price + charging_price;
				path_to_cs <- path_between(road_network_weighted, self, myself);
				access_time <- (path_to_cs.weight) / 60; // this converts the time to minutes
				if (self.restaurents > 0) {
					amenity_more <- 1;
				} else {
					amenity_more <- 0;
				}

			}

			float
			u_charging <- intercept + (c_soc * soc) + (c_dev * dev) + (c_time_in_car * time_in_car) + (c_charging_cost * charging_cost_in_dollar) + (c_charging_time * charging_time) + (c_access_time * access_time) + (c_amenity_restroom * amenity_restroom) + (c_amenity_more * amenity_more);
			float odds_charging <- exp(u_charging);
			prob_charging <- odds_charging / (1 + odds_charging);
			// set the seed every time to ensure repeatability
			// seed <- float(params[2][0][2]);
			// Make a random draw using the probability from a binomial distribution					
			return bool(binomial(1, prob_charging));
		}

	}

	// aspect definitions
	aspect circle {
		draw square(2000) color: #green;
	} }

species charging_station {
	string station_id; // station_id as in AFDC dataset
	int dcfc_plug_count; // Number of DCFC plugs at the location
	string ev_network; // Network company that the charging station belongs to, for ex: Blink, ChargePoint etc.
	string ev_connector_types; // Type of connector - CHADEMO (1), COMBO (2), BOTH (3), TESLA (4)
	float max_power <- 50.0; // max charging power per charging station in kWhr
	float charging_cost; // charging cost in $ per kWhr
	int ev_connector_code; // Codes as indicated in connector_type parenthesis
	rgb cs_color <- #red;
	int zip;
	int restaurents <- 1; // Denoting the number of restaurents in the vicinity of the CS
	float energy_consumed <- 0.0;
	bool charger_available <- true;
	int plugs_in_use <- 0;
	float current_power_draw <- plugs_in_use * max_power;
	int evs_passed;
	list<EVs> waiting_evs;
	list<EVs> queued_evs;
	float dcfc_fixed_charging_price;
	string dcfc_var_charging_price_unit;
	float dcfc_var_charging_price;
	float dcfc_fixed_parking_price;
	string dcfc_var_parking_price_unit;
	float dcfc_var_parking_price;

	// aspect definitions
	aspect circle {
		draw square(1000) color: cs_color;
	}

	reflex update_power_draw {
		current_power_draw <- plugs_in_use * max_power;
	}

}

experiment no_gui_exp {
	output {
	}

}


experiment gui_exp {

// Define parameters here if necessary
// parameter "My parameter" category: "My parameters" var: one_global_attribute;
	// float seed <- 123.0;
	float seed <- float(aidm[0, 1]) parameter: true;
	parameter var: analysis_id name: "analysis_id" category: "My parameters";
	// parameter var: seed name: "seed" category: "My parameters";
	parameter var: db_host name: "db_host" category: "My parameters";
	parameter var: db_type name: "db_type" category: "My parameters";
	parameter var: db_name name: "db_name" category: "My parameters";
	parameter var: db_port name: "db_port" category: "My parameters";
	parameter var: db_user name: "db_user" category: "My parameters";
	parameter var: db_pwd name: "db_pwd" category: "My parameters";
	// Define attributes, actions, a init section and behaviors if necessary
	
//	 reflex getseed {
//        ask simulations {write seed;}
//    } 
//	init { 
//		create simulation with:[seed::seedValue];
//		
//	}
	output {
		display WA_network type: opengl {
			species road aspect: geom refresh: false;
			species EVs aspect: circle;
			species charging_station aspect: circle refresh: false;
			graphics output_overlay {
				draw string(current_date) at: {world.shape.width / 1.2, world.shape.height} color: #blue; // shows the clock at the bottom right edge of display
				loop ii from: 0 to: length(plot_compat_chargers collect each.location) {
				// draw square(1000) color: #magenta at: plot_compat_chargers[ii].location;
				}

				loop jj from: 0 to: length(plot_chargers_nearby collect each.location) {
				// draw square(1000) color: #blue at: plot_chargers_nearby[jj].location;
				}

				//				draw square(5000) color: #yellow at: nearest_evses[0];
				//				draw square(5000) color: #brown at: nearest_evses[1];
				loop jj from: 0 to: length(plot_chargers_nearby collect each.location) {
					draw circle(1100) color: #yellow at: plot_chargers_nearby[jj].location;
					draw square(1100) color: #turquoise at: plot_cpts[jj];
					draw string(int(plot_chargers_nearby[jj])) at: plot_chargers_nearby[jj].location size: 1 #m color: #black;
				}

                
				draw geometry(plot_shortest_path.segments) color: #green width: 4;
				draw geometry(my_path) color: #blue width: 50;
				// draw circle(1100) color: #yellow at: plot_chargers_nearby[0].location;
				draw geometry(plot_my_circle) color: #darkorange empty: true width: 10;
				// draw string(int(nearest_evse)) at: nearest_evse.location size: 5#m color: #black;

				//                loop rr from: 0 to: length(cpts_on_path) - 1 {
				//                    draw circle(500) at: cpts_on_path[rr] color: #orange ;
				//                }
			}

		}

	}

}

