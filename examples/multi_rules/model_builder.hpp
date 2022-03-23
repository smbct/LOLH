// model_builder.hpp

#include <ghost/model_builder.hpp>

#include "../../c++/src/DataFrame.hpp"

class LPOBuilder : public ghost::ModelBuilder {

		int _p_rules;
		double _t;


	public:

		/*!
		 * \brief default constructor
		*/
		LPOBuilder(int p_rules, double t, DataFrame<uint>& dataset);

		void declare_variables() override;
		void declare_objective() override;

};
