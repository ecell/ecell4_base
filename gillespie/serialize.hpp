
#include "pficommon/text/json.h"
#include "pficommon/text/csv.h"
#include "GillespieWorld.hpp"
#include "GillespieSolver.hpp"

#ifndef	INCLUDE_GUARD_SERIALIZE
# define	INCLUDE_GUARD_SERIALIZE

using namespace pfi::text::json;
using namespace std;

json string_to_json(std::string str);
string read_file_all(const char *json_filename);

World *init_world_from_json(json js_world);

World *init_world_from_csv(string csv_str);

ReactionRule *init_reaction_from_json(json js_reaction);

Model *init_model_from_json(json js_model);

Model *init_model_from_csv(string &csv_model);


#endif
