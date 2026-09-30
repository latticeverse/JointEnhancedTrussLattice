#pragma once
#include "Truss.h"
#include <string>

// Configure PTL before calling. The original objectives are retained in plsm.cpp.
void plsmtoSolid(Truss& PTL, double trd, int id, double begin_smooth_factor,
                 const std::string& skeleton_path, int max_iterations = 150);
void plsmtoSolidETH(Truss& PTL, double trd, int id, double begin_smooth_factor,
                    const std::string& skeleton_path, int max_iterations = 150);
void plsmtoSolidModified(Truss& PTL, double trd, int id, double begin_smooth_factor,
                         const std::string& skeleton_path, int max_iterations = 120);
void plsmtoSolidModifiedETH(Truss& PTL, double trd, int id, double begin_smooth_factor,
                            const std::string& skeleton_path, int max_iterations = 120);
