#ifndef INCOMPRESSCOUETTEFLOW_H
#define INCOMPRESSCOUETTEFLOW_H

#include <iostream>
#include <set>

void ImplicitIncompressCouetteFlow(double E,
	                               double ReD,
								   size_t N,
								   size_t IterateTimes,
								   std::set<size_t> SelectedIterateTimes,
								   std::ostream &out);

void ExplicitIncompressCouetteFlow(double E,
	                               double ReD,
	                               size_t N,
	                               size_t IterateTimes,
	                               std::set<size_t> SelectedIterateTimes,
	                               std::ostream &out);

#endif // INCOMPRESSCOUETTEFLOW_H