#include "ma.h"

int main(int argc, char* argv[])
{
	std::vector <ST_RESULT_DATA> RESULTS;

	auto results = met_aromatic("1rcy.pdb");   /* can use auto in different namespace */

	RESULTS = std::get<1>(results);
	std::cout << "Status byte: " << std::get<0>(results) << '\n';
	std::cout << "Size of output: " << RESULTS.size() << '\n';
	std::cout << "Total execution time: " << std::get<2>(results) << " s" << '\n';

	/* do something like print to console or export to .csv */
	for (int i = 0; i < RESULTS.size(); ++i)
	{
		std::cout << RESULTS[i].MET_RES << " ";
		std::cout << RESULTS[i].ARO_RES << " ";
		std::cout << RESULTS[i].NORM_V << " ";
		std::cout << RESULTS[i].MET_THETA << " ";
		std::cout << RESULTS[i].MET_PHI << " ";
		std::cout << std::endl;
	}

	return 0;
}