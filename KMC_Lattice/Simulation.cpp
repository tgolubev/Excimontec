// Copyright (c) 2018 Michael C. Heiber
// This source file is part of the KMC_Lattice project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The KMC_Lattice project can be found on Github at https://github.com/MikeHeiber/KMC_Lattice

#include "Simulation.h"

using namespace std;
using namespace Utils;

Simulation::Simulation() {

}

Simulation::~Simulation() {

}

void Simulation::init(const Parameters_Simulation& params, const int id) {
	Id = id;
	// General Parameters
	Enable_logging = params.Enable_logging;
//FOR NOW MANUALLY ENABLE LOGGING, until figure out how to use the params
//Enable_logging = true;
	// Event recalculation method parameters
	Enable_FRM = params.Enable_FRM;
	Enable_selective_recalc = params.Enable_selective_recalc;
	Recalc_cutoff = params.Recalc_cutoff;
	Enable_full_recalc = params.Enable_full_recalc;
	// Lattice Parameters
	Parameters_Lattice params_lattice;
	params_lattice.Enable_periodic_x = params.Enable_periodic_x;
	params_lattice.Enable_periodic_y = params.Enable_periodic_y;
	params_lattice.Enable_periodic_z = params.Enable_periodic_z;
	params_lattice.Length = params.Length;
	params_lattice.Width = params.Width;
	params_lattice.Height = params.Height;
	params_lattice.Unit_size = params.Unit_size;
	temperature = params.Temperature;
	// Initialize data structures
	lattice.init(params_lattice, &generator);
	object_ptrs.clear();
	event_ptrs.clear();
	generator.seed((int)time(0)*(id + 1));
	// Output files
	Logfile = params.Logfile;
}


list<Event*>::const_iterator Simulation::addEvent(Event* event_ptr) {
	event_ptrs.push_back(event_ptr);
	return --event_ptrs.end();
}

void Simulation::addObject(Object* object_ptr) {
	// Add an event for the object to the event list and link the event to the object
	event_ptrs.push_back(nullptr);
	object_ptr->setEventIt(--event_ptrs.end());
	// Add new object to the object vector and link the object to the event
	object_ptrs.push_back(object_ptr);
	// Set occupancy of site
	lattice.setOccupied(object_ptr->getCoords());
	(*lattice.getSiteIt(object_ptr->getCoords()))->setObjectPtr(*(--object_ptrs.end()));
	// Update counters
	N_objects_created++;
	N_events_executed++;
}

list<Event*>::const_iterator Simulation::chooseNextEvent() {
	auto event_target_it = event_ptrs.begin();
	if (event_ptrs.size() > 1) {
		for (auto it = ++event_ptrs.begin(); it != event_ptrs.end(); ++it) {
			if ((*it) != nullptr && ((*event_target_it) == nullptr || ((*it)->getExecutionTime() < (*event_target_it)->getExecutionTime()))) {
				event_target_it = it;
			}
		}
	}
	return event_target_it;
}

vector<Object*> Simulation::findRecalcObjects(const Coords& coords_start, const Coords& coords_dest) const {
	vector<Object*> object_recalc_ptrs;
	object_recalc_ptrs.reserve(10);
	if (Enable_FRM) {
		for (auto const &item : object_ptrs) {
			// Recalculate event for an Object located at coords_start or coords_dest
			if (item->getCoords() == coords_start || item->getCoords() == coords_dest) {
				object_recalc_ptrs.push_back(item);
			}
			// Recalculate event for an Object that has a valid queued event that targets coords_start or coords_dest
			else if (*(item->getEventIt())!=nullptr && ((*(item->getEventIt()))->getDestCoords() == coords_start || (*(item->getEventIt()))->getDestCoords() == coords_dest)) {
				object_recalc_ptrs.push_back(item);
			}
		}
	}
	else if (Enable_selective_recalc) {
		object_recalc_ptrs = findRecalcNeighbors(coords_start);
		if (coords_start != coords_dest) {
			vector<Object*> object_recalc_ptrs2 = findRecalcNeighbors(coords_dest);
			object_recalc_ptrs.insert(object_recalc_ptrs.end(), object_recalc_ptrs2.begin(), object_recalc_ptrs2.end());
			removeDuplicates(object_recalc_ptrs);
		}
	}
	else if (Enable_full_recalc) {
		for (auto const &item : object_ptrs) {
			object_recalc_ptrs.push_back(item);
		}
	}
	return object_recalc_ptrs;
}

vector<Object*> Simulation::findRecalcNeighbors(const Coords& coords) const {
	vector<Object*> neighbor_ptrs;
	neighbor_ptrs.reserve(10);
	Coords coords2;
	const static int recalc_cutoff_sq_lat = (int)((Recalc_cutoff / lattice.getUnitSize())*(Recalc_cutoff / lattice.getUnitSize()));
	int distance_sq_lat;
	for (auto const &item : object_ptrs) {
		coords2 = item->getCoords();
		distance_sq_lat = lattice.calculateLatticeDistanceSquared(coords, coords2);
		if (distance_sq_lat <= recalc_cutoff_sq_lat) {
			neighbor_ptrs.push_back(item);
		}
	}
	return neighbor_ptrs;
}

vector<Object*> Simulation::getAllObjectPtrs() const {
	vector<Object*> object_ptrs_vec;
	object_ptrs_vec.reserve(object_ptrs.size());
	for (auto const &item : object_ptrs) {
		object_ptrs_vec.push_back(item);
	}
	return object_ptrs_vec;
}

string Simulation::getErrorMessage() const {
	return error_msg;
}

int Simulation::getId() const {
	return Id;
}

int Simulation::getN_events() const {
	int count = 0;
	for (auto const &item : event_ptrs) {
		if (item!=nullptr && item->getExecutionTime() > 0) {
			count++;
		}
	}
	return count;
}

long int Simulation::getN_events_executed() const {
	return N_events_executed;
}

int Simulation::getTemp() const {
	return temperature;
}

double Simulation::getTime() const {
	return time_sim;
}

double Simulation::getVolume() const {
	return lattice.getVolume();
}

bool Simulation::isLoggingEnabled() const {
	return Enable_logging;
}

void Simulation::moveObject(Object* object_ptr, const Coords& coords_dest) {
	Coords coords_initial = object_ptr->getCoords();
	// Clear occupancy of initial site
	lattice.clearOccupancy(coords_initial);
	// Check for periodic boundary crossing
	object_ptr->incrementDX(-lattice.calculateDX(coords_initial, coords_dest));
	object_ptr->incrementDY(-lattice.calculateDY(coords_initial, coords_dest));
	object_ptr->incrementDZ(-lattice.calculateDZ(coords_initial, coords_dest));
	// Set object coords to new site
	object_ptr->setCoords(coords_dest);
	// Set occupancy of new site
	lattice.setOccupied(coords_dest);
	(*lattice.getSiteIt(coords_dest))->setObjectPtr(object_ptr);
	// Update counter
	N_events_executed++;
}

double Simulation::rand01() {
	return generate_canonical<double, std::numeric_limits<double>::digits>(generator);
}

void Simulation::removeEvent(Event* event_ptr) {
	bool success = false;
	for (auto it = event_ptrs.begin(); it != event_ptrs.end(); ++it) {
		if (*it == event_ptr) {
			event_ptrs.erase(it);
			success = true;
			break;
		}
	}
	if (!success) {
		cout << "Error! The Event pointer could not be found in the event list and could not be removed." << endl;
		Error_found = true;
	}
}

void Simulation::removeObject(Object* object_ptr) {
	// Clear occupancy of site
	lattice.clearOccupancy(object_ptr->getCoords());
	// Delete the event pointer
	removeEvent(*object_ptr->getEventIt());
	// Delete the object pointer
	bool success = false;
	for (auto it = object_ptrs.begin(); it != object_ptrs.end(); ++it){
		if (*it == object_ptr) {
			object_ptrs.erase(it);
			success = true;
			break;
		}
	}
    // Update counters
    N_events_executed++;
	// Error checking
	if (!success) {
		cout << "Error! The Object pointer could not be found in the object list and could not be removed." << endl;
		Error_found = true;
	}
}

void Simulation::setErrorMessage(const string& input_msg) {
	error_msg = input_msg;
}

void Simulation::setObjectEvent(const Object* object_ptr,Event* event_ptr){
	*(object_ptr->getEventIt()) = event_ptr;
}

void Simulation::setGeneratorSeed(int seed) {
	generator.seed(seed);
}

void Simulation::setTime(const double new_time){
    time_sim = new_time;
}

