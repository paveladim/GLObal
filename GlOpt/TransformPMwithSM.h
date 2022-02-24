#ifndef TRANSFORMPMWITHSM_H
#define TRANSFORMPMWITHSM_H

#include "TransformPM.h"
#include "SimplexMethod.h"

class TransformPMwithSM : public TransformPM {
private:
	CoordinatesValues _point;
	std::vector<double> _incs;
	std::vector<double> _non_proj_incs;
	Matrix _st;
	Vec _b;
	Vec _c;
private:
	void calculate_localLipshConst(const uint& id_hyp) override;
	void calculate_characteristic(const uint& id_hyp) override;
private:
	void decode_and_save(const uint& pos, const uint& order);
	void projection(const uint& order,
					std::vector<double>& e1,
					std::vector<double>& e2);
	void calculate_and_project(const uint& axis);
	double scalar_product(const std::vector<double>::iterator& a,
						  const std::vector<double>::iterator& b);
	void generate_simplex_table();
	void generate_right_part(const uint& function,
							 const uint& eval_a,
							 const uint& eval_v,
							 const uint& eval_u,
							 const uint& eval_b);
	size_t index(const size_t& ind);
public:
	TransformPMwithSM(const uint& dimension,
					  const uint& constraints,
					  Parameters& parameters,
					  Problem& problem);
	~TransformPMwithSM() {}
};

#endif // TRANSFORMPMWITHSM_H