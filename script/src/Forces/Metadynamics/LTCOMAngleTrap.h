#ifndef LTCOMANGLETRAP_H_
#define LTCOMANGLETRAP_H_

#include "../BaseForce.h"
#include <vector>

class LTCOMAngleTrap: public BaseForce {
private:
	std::vector<int> _p1a;
	std::vector<int> _p2a;
	std::vector<int> _p3a;

public:
	std::vector<BaseParticle*> _p1a_ptr;
	std::vector<BaseParticle*> _p2a_ptr;
	std::vector<BaseParticle*> _p3a_ptr;
	number xmin;
	number xmax;
	int N_grid;
	number dX;
	std::vector<number> potential_grid;
	int _mode = 0;
	bool PBC = false;

	LTCOMAngleTrap();
	virtual ~LTCOMAngleTrap() {

	}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	LR_vector value(llint step, LR_vector &pos) override;
	number potential(llint step, LR_vector &pos) override;

protected:
	LR_vector _distance(LR_vector u, LR_vector v);
};

#endif // LTCOMANGLETRAP_H_
