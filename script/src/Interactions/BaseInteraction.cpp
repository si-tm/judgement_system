/*
 * BaseInteraction.cpp
 *
 *  Created on: 11 ott 2019
 *      Author: lorenzo
 */

#include "BaseInteraction.h"

BaseInteraction::BaseInteraction() {
	_energy_threshold = (number) 100.f;
	_is_infinite = false;
	_box = NULL;
	_generate_consider_bonded_interactions = true;
	_generate_bonded_cutoff = 2.0;
	_temperature = 1.0;

	_rcut = 2.5;
	_sqr_rcut = SQR(_rcut);
}

BaseInteraction::~BaseInteraction() {

}

number BaseInteraction::pair_interaction_term(int name, BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	// _int_map.find(name) returns an iterator which points to the
	// requested function pointer. We dereference it, apply it to the
	// that, which is a pointer to an instance of a class which inherits
	// from BaseInteraction. We then pass in the required parameters and
	// return the interaction energy.
	typename interaction_map::iterator interaction = _interaction_map.find(name);
	if(interaction == _interaction_map.end()) {
		throw oxDNAException("%s, line %d: Interaction term '%d' not found", __FILE__, __LINE__, name);
	}

	return interaction->second(p, q, false, update_forces);
}

std::map<int, number> BaseInteraction::get_system_energy_split(std::vector<BaseParticle *> &particles, BaseList *lists) {
	begin_energy_computation();

	std::map<int, number> energy_map;

	for(auto it = _interaction_map.begin(); it != _interaction_map.end(); it++) {
		int name = it->first;
		energy_map[name] = (number) 0.f;
	}

	for(auto p: particles) {
		std::vector<BaseParticle *> neighs = lists->get_all_neighbours(p);

		for(unsigned int n = 0; n < neighs.size(); n++) {
			BaseParticle *q = neighs[n];
			if(p->index > q->index) {
				for(auto it = _interaction_map.begin(); it != _interaction_map.end(); it++) {
					int name = it->first;
					energy_map[name] += pair_interaction_term(name, p, q);
				}
			}
		}
	}

	return energy_map;
}

void BaseInteraction::get_settings(input_file &inp) {
	getInputString(&inp, "topology", _topology_filename, 1);
	getInputNumber(&inp, "energy_threshold", &_energy_threshold, 0);
	getInputNumber(&inp, "T", &_temperature, 1);
	getInputBool(&inp, "generate_consider_bonded_interactions", &_generate_consider_bonded_interactions, 0);
	if(_generate_consider_bonded_interactions) {
		getInputNumber(&inp, "generate_bonded_cutoff", &_generate_bonded_cutoff, 0);
		OX_LOG(Logger::LOG_INFO, "The generator will try to take into account bonded interactions by choosing distances between bonded neighbours no larger than %lf", _generate_bonded_cutoff);
	}
}

int BaseInteraction::get_N_from_topology() {
	std::ifstream topology;
	topology.open(_topology_filename, std::ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename);
	int ret;
	topology >> ret;
	topology.close();
	return ret;
}

void BaseInteraction::read_topology(int *N_strands, std::vector<BaseParticle *> &particles) {
	*N_strands = particles.size();
	allocate_particles(particles);
	int idx = 0;
	for(auto p : particles) {
		p->index = idx;
		p->type = 0;
		p->strand_id = idx;
		idx++;
	}
}

void BaseInteraction::begin_energy_computation() {
	if(has_custom_stress_tensor()) {
		reset_stress_tensor();
	}
}

void BaseInteraction::_update_stress_tensor(LR_vector r_p, LR_vector group_force) {
	_stress_tensor[0] += r_p[0] * group_force[0];
	_stress_tensor[1] += r_p[1] * group_force[1];
	_stress_tensor[2] += r_p[2] * group_force[2];
	_stress_tensor[3] += r_p[0] * group_force[1];
	_stress_tensor[4] += r_p[0] * group_force[2];
	_stress_tensor[5] += r_p[1] * group_force[2];
}

void BaseInteraction::reset_stress_tensor() {
	std::fill(_stress_tensor.begin(), _stress_tensor.end(), 0.);
}

number BaseInteraction::get_system_energy(std::vector<BaseParticle *> &particles, BaseList *lists) {
	begin_energy_computation();

	double energy = 0.;
	std::vector<ParticlePair> pairs = lists->get_potential_interactions();
	for(auto &pair : pairs) {
		BaseParticle *p = pair.first;
		BaseParticle *q = pair.second;
		energy += (double) pair_interaction(p, q);
		if(get_is_infinite()) {
			return energy;
		}
	}

	return (number) energy;
}

number BaseInteraction::get_system_energy_term(int name, std::vector<BaseParticle *> &particles, BaseList *lists) {
	begin_energy_computation();

	number energy = (number) 0.f;
	for(auto p : particles) {
		std::vector<BaseParticle *> neighs = lists->get_all_neighbours(p);

		for(unsigned int n = 0; n < neighs.size(); n++) {
			BaseParticle *q = neighs[n];
			if(p->index > q->index) energy += pair_interaction_term(name, p, q);
			if(get_is_infinite()) return energy;
		}
	}

	return energy;
}

bool BaseInteraction::generate_random_configuration_overlap(BaseParticle *p, BaseParticle *q) {
	_computed_r = _box->min_image(p, q);

	if(_computed_r.norm() >= _sqr_rcut) {
		return false;
	}

	number energy = pair_interaction_nonbonded(p, q, false, false);

	// in case we create an overlap, we reset the interaction state
	set_is_infinite(false);

	// if energy too large, reject
	return energy > _energy_threshold;
}

void BaseInteraction::generate_random_configuration(std::vector<BaseParticle *> &particles) {
	Cells c(particles, _box);
	c.init(_rcut);

	for(auto p : particles) {
		p->pos = LR_vector(drand48() * _box->box_sides().x, drand48() * _box->box_sides().y, drand48() * _box->box_sides().z);
	}

	c.global_update();

	int N = particles.size();
	for(int i = 0; i < N; i++) {
		BaseParticle *p = particles[i];
		bool same_strand = (_generate_consider_bonded_interactions && i > 0 && p->is_bonded(particles[i - 1]));

		bool inserted = false;
		do {
			if(same_strand) {
				p->pos = particles[i - 1]->pos + LR_vector((drand48() - 0.5), (drand48() - 0.5), (drand48() - 0.5)) * _generate_bonded_cutoff;
			}
			else {
				p->pos = LR_vector(drand48() * _box->box_sides().x, drand48() * _box->box_sides().y, drand48() * _box->box_sides().z);
			}
			// random orientation
			//p->orientation = Utils::get_random_rotation_matrix (2.*M_PI);
			p->orientation = Utils::get_random_rotation_matrix_from_angle(acos(2. * (drand48() - 0.5)));
			p->orientation.orthonormalize();
			p->orientationT = p->orientation.get_transpose();

			p->set_positions();
			p->set_ext_potential(0, _box);
			c.single_update(p);

			inserted = true;

			// we take into account the bonded neighbours
			for(unsigned int n = 0; n < p->affected.size(); n++) {
				BaseParticle *p1 = p->affected[n].first;
				BaseParticle *p2 = p->affected[n].second;
				if(p1->index <= p->index && p2->index <= p->index) {
					number e = pair_interaction_bonded(p1, p2);
					if(std::isnan(e) || e > _energy_threshold) inserted = false;
				}
			}

			// here we take into account the non-bonded interactions
			std::vector<BaseParticle *> neighs = c.get_complete_neigh_list(p);
			for(unsigned int n = 0; n < neighs.size(); n++) {
				BaseParticle *q = neighs[n];
				// particles with an index larger than p->index have not been inserted yet
				if(q->index < p->index && generate_random_configuration_overlap(p, q)) inserted = false;
			}

			// we take into account the external potential
			number boltzmann_factor = exp(-p->ext_potential / _temperature);
			if(std::isnan(p->ext_potential) || drand48() > boltzmann_factor) {
				inserted = false;
			}

		} while(!inserted);

		if(i > 0 && N > 10 && i % (N / 10) == 0) OX_LOG(Logger::LOG_INFO, "Inserted %d%% of the particles (%d/%d)", i*100/N, i, N);
	}
}
