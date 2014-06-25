#define BOOST_TEST_MODULE "Hogehoge"

#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>


#include "pficommon/text/json.h"
#include "pficommon/text/csv.h"
#include "GillespieWorld.hpp"
#include "GillespieSolver.hpp"
#include "serialize.hpp"

#include <vector>
#include <string>


using namespace std;

BOOST_AUTO_TEST_CASE( gillespie_solver_test )
{
	string json_file_content(read_file_all("./data/init.json"));
	World *w = init_world_from_json( string_to_json(json_file_content));
	std::cout << *w << std::endl;

	string json_model(read_file_all("./data/reaction.json"));
	Model *m = init_model_from_json( string_to_json(json_model));


	GillespieSolver gs(*w, *m);
	for(double elapse = 0.0; elapse < 10.0; elapse += gs.run(1.0)) {
		std::cout << *w << std::endl;
	}
}

