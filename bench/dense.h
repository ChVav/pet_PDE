#pragma once

#include <vector>
#include <string>

#include "common.h"

std::vector<double> assemble_matrix(
	const std::vector<point>& points,
	const std::vector<triangle>& triangles,
	const std::vector<line>& lines,
	const std::vector<double>& areas);

std::vector<double> one_timestep(
	double dt,
	const std::vector<double>& B,
	const std::vector<double>& u_n);