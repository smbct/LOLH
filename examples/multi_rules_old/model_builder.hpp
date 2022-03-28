// model_builder.hpp

#include <ghost/model_builder.hpp>

#include "../../c++/src/DataFrame.hpp"

#include "instance_data.hpp"

class LPOBuilder : public ghost::ModelBuilder {

		instance_data _instance;

	public:

		/*!
		 * \brief default constructor
		*/
		LPOBuilder();

		void declare_variables() override;
		void declare_constraints() override;
		void declare_objective() override;

};
