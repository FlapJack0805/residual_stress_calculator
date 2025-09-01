#pragma once

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Kd_tree.h>

using Kernel = CGAL::Simple_cartesian<double>;
using Traits = CGAL::Search_traits_2<Kernel>;
using KD_tree =  CGAL::Kd_tree<Traits>;
using Search = CGAL::Orthogonal_k_neighbor_search<Traits>;

struct unit_vector
{
	double x, y;
};
