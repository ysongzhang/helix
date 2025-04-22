#pragma once
#include "Math/gfpMatrix.h"

namespace hmmpc
{
void funcCheckMaliciousMatMulFalcon(const gfpMatrix &x, const gfpMatrix &y, const gfpMatrix &z);
void funcCheckMaliciousWiseMulFalcon(const gfpMatrix &x, const gfpMatrix &y, const gfpMatrix &z);

void funcCheckMaliciousMatMulLN17(const gfpMatrix &x, const gfpMatrix &y, const gfpMatrix &z);
}