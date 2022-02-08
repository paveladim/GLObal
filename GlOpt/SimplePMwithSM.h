#ifndef SIMPLEPMWITHSM_H
#define SIMPLEPMWITHSM_H

#include "SimplePM.h"
#include "SimplexMethod.h"

class SimplePMwithSM : public SimplePM {
private:
	void calculate_localLipshConst(const uint& id_hyp) override;
	void calculate_characteristic(const uint& id_hyp) override;
private:
	void calculate_and_project(const CoordinatesValues& out,
							   std::vector<double>& incs,
							   const uint& axis);
	double scalar_product(const std::vector<double>& a,
						  const std::vector<double>& b);
	void generate_simplex_table(std::vector<std::vector<double>>& A,
								const std::vector<double>& incs);
	void generate_right_part(std::vector<double>& b,
							 const uint& function,
							 const uint& eval_a,
							 const uint& eval_v,
							 const uint& eval_u,
							 const uint& eval_b);
public:
	SimplePMwithSM(const uint& dimension,
				   const uint& constraints,
				   Parameters& parameters,
				   Problem& problem);
	~SimplePMwithSM() {}
};

#endif // SIMPLEPMWITHSM_H
