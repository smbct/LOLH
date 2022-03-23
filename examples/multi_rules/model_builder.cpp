#include "model_builder.hpp"


using namespace ghost;

/*-----------------------------------------------------------------------------*/
LPOBuilder::LPOBuilder(int p_rules, double t, DataFrame<uint>& dataset) : _p_rules(_p_rules), _t(t), ghost::ModelBuilder() {


}

/*-----------------------------------------------------------------------------*/
void LPOBuilder::declare_variables() {

  // for(int i = 0; i < _n*_n; i ++) {
  //   variables.emplace_back( 1, _n*_n, " ");
  // }

}


/*-----------------------------------------------------------------------------*/
void LPOBuilder::declare_objective() {

  // for(int i = 0; i < _n*_n; i ++) {
  //   variables.emplace_back( 1, _n*_n, " ");
  // }

}
