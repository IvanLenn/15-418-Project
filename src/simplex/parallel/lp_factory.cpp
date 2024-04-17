#include "lp_par.h"
#include "lp_par1.h"

LinearProgramming* createLP(const std::string type) {
    return new LinearProgramming1();
}
