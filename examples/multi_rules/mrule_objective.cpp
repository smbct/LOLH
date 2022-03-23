#include "mrule_objective.hpp"


/*-----------------------------------------------------------------------------*/
MRuleObjective::MRuleObjective(const std::vector< ghost::Variable >& var) : ghost::Maximize(var,  "Multi rule scores" )
{


}

/*-----------------------------------------------------------------------------*/
double MRuleObjective::required_cost( const std::vector< ghost::Variable*>& variables ) const {

  // Notice the minus here.
  // GHOST's solver tries to minimize any objective function.
  // Thus, for maximization problems like this one, outputing '- returned_value' does the trick.

  return 0;
}
