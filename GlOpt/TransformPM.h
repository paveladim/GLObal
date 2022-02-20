#ifndef TRANSFORMPM_H
#define TRANSFORMPM_H

#include "DivideByThree.h"
#include "SimplexMethod.h"

class TransformPM : public DivideByThree {
protected:
    bool _areAllCharInfty;
    bool _doesGlobalChange;
    FunctionsValues _localLipshEval;
    FunctionsValues _globalLipshEval;
private:
	CoordinatesValues _points;
	std::vector<double> _incs;
	std::vector<double> _non_proj_incs;
	Matrix _st;
	Vec _b;
	Vec _c;
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
private:
	double calculate_residual(const double& t, const uint& id_hyp);
    void calculate_localLipshConst(const uint& id_hyp);
	void calculate_characteristic(const uint& id_hyp) override;
    FunctionValue give_residual(const FunctionsValues& evals);
    void calculate_globalLipshConst(const uint& id_hyp);
    double mixedLipEval(const Hyperinterval& hyp, const uint& i);
    void update_all_charact();
    uint optimal_to_trisect() override;
    uint iterate(const uint& id_hyp) override;
public:
	TransformPM(const uint& dimension,
                const uint& constraints,
                Parameters& parameters,
                Problem& problem);
    ~TransformPM() {}
};

#endif // TRANSFORMPM_H

