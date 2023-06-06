/**
 * \file		shared/src/potential_AGAMA.cpp
 * \brief		Potential output from AGAMA.
 * \author	Rimpei Chiba
 * \date		2022-
 */
/*
#include "potential_AGAMA.hpp"

namespace RC
{
	AGAMA::AGAMA(const std::string in):
		multipole(in),
		cartGrid(in)
	{
	}

	double AGAMA::pot(const Polar &pos)
	{
		return multipole.pot(pos) + cartGrid.pot(pos);
	}

	Force AGAMA::force(const Polar &pos)
	{
		return multipole.force(pos) + cartGrid.force(pos);
	}

	double AGAMA::pot(double r)
	{
		return multipole.pot(r) + cartGrid.pot(r);
	}

	double AGAMA::dPhidr(double r)
	{
		return multipole.dPhidr(r) + cartGrid.dPhidr(r);
	}

	double AGAMA::d2Phidr2(double r)
	{
		return multipole.d2Phidr2(r) + cartGrid.d2Phidr2(r);
	}
}
*/
