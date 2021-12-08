#ifndef UNIFORM_H
#define UNIFORM_H

#include "synonymous_types.h"
#include "DivideByThree.h"

class Uniform : public DivideByThree {
    void calculate_characteristic(const uint& id_hyp) override;
    uint optimal_to_trisect() override;
    uint iterate(const uint& id_hyp) override;
public:
    Uniform(const uint& dimension,
            const uint& constraints,
            Parameters& parameters,
            Problem& problem);
    void solve() override;
    ~Uniform() {}
};

#endif // UNIFORM_H