#pragma once

#include <ghost/objective.hpp>

class MRuleObjective : public ghost::Maximize {

    double required_cost( const std::vector< ghost::Variable*>& ) const override;

  public:

    MRuleObjective(const std::vector< ghost::Variable >&);

};
