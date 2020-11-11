/**
* Name: pointtest
* Based on the internal empty template. 
* Author: Chintan Pathak
* Tags: 
*/
model pointtest

/* Insert your model definition here */

/* Insert your model definition here */
global {

	init {
		list<geometry> roads_list <- list(polyline({0, 0}, {100, 100}));
		create road from: roads_list;
        graph road_network;
        road_network <- as_edge_graph(road);
		//
		create EV_station number: 5;
		loop ev over: EV_station {
			create EV_station_on_network number: 1 {
				location <- ev.location;
				ev.dist_to_road <- distance_to(self, road closest_to ev);
				ev.my_circle <- circle(ev.dist_to_road + 0.1);
				write(ev.my_road);
				location <- point(one_of(ev.my_circle inter ev.my_road));
				write(location);
				my_parent <- ev.name;
			}

		}
        write(distance_between(topology(road_network), [EV_station_on_network[0], EV_station_on_network[1]]));
	}

}

species road {

	aspect default {
		draw shape + 0.05 color: #black;
	}

}

species EV_station {
	rgb color <- #green;
	geometry my_road;
	float dist_to_road;
	geometry my_circle;

	init {
		my_road <- road closest_to self;
	}

	aspect default {
		draw circle(0.5) color: color;
		draw string((self)) color: #black at: self.location + {-1.0, -1.0, -1.0};
		draw my_circle color: #lightgray empty: true;
	}

}

species EV_station_on_network {
	rgb color <- #red;
	string my_parent;

	aspect default {
		draw square(0.5) color: color;
		draw "proj_" + my_parent color: #black at: self.location + {1.0, 1.0, 1.0};
	}

}

experiment name type: gui {
	output {
		display "My display" {
			species EV_station;
			species EV_station_on_network;
			species road;
		}

	}

}