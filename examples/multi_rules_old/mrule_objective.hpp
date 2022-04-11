#pragma once

#include <ghost/objective.hpp>

#include "instance_data.hpp"

class MRuleObjective : public ghost::Maximize {

    double required_cost(const std::vector< ghost::Variable*>& ) const override;



    instance_data& _instance;

    std::vector<std::pair<uint,uint>> _atoms;

    /* for each rule, errors of all the atoms */
    std::vector<std::vector<std::pair<uint,uint>>> _atom_rule_errors;

    /* values of the previous variable assignment */
    std::vector<uint> _prev_values;


  private: /* private methods */

    /*!
     * \brief compute the score of a rule
     * \param the rule index, between 0 and _instance.p_rules-1
     * \return the score of the rule (normalized by the number of selected atoms)
     */
    double compute_rule_score(const std::vector< ghost::Variable*>& variables, uint rule_id) const;


  protected:
    
    void conditional_update_data_structures(const std::vector<ghost::Variable*> &variables, int index, int new_value) override;

  public:

    /*!
     * \brief default constructor
    */
    MRuleObjective(instance_data& instance, const std::vector< ghost::Variable >&);

};
