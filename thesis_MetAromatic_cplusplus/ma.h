#ifndef MA_H
#define MA_H

/*
Written by dsw7@sfu.ca
Program works properly as of Oct 13, 2018
I only tested against 1rcy though
*/

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <regex>
#include <cmath>
#include <map>
#include <set>
#include <ctime>
#include <tuple>

#define DIST_CUTOFF 6.0
#define ANGLE_CUTOFF 109.5
#define PI 3.14159
#define HC 180.000
#define PAT_MET "CE|SD|CG"
#define PAT_PHE_TYR "CG|CD2|CE2|CZ|CE1|CD1"
#define PAT_TRP "CD2|CE3|CZ3|CH2|CZ2|CE2"
#define SUCCESS_CODE 0
#define NO_RESIDUE_ERR -1
#define NO_RESULT_ERR -2


/* ------- */
/* structs */
/* ------- */
struct rodrigues_rotation_output
{
	// needed for rodrigues function output
	std::vector <float> vector_a;
	std::vector <float> vector_g;
};


struct Line
{
	// holds all ATOM lines directly from PDB file
	std::string ATOM;
	int ATOM_NUMBER;
	std::string ATOM_ID;
	std::string RES;
	std::string CHAIN;
	int RES_POS;
	std::vector <float> v;
};


struct MET_LP_STRUCT
{
	// holds data for methionine lone pairs
	int RES_POS;                        // the MET residue position
	std::vector <float> v;              // the MET SD coordinates
	std::vector <float> lone_pair_a;    // the lone pair a vector mapped to origin
	std::vector <float> lone_pair_g;    // the lone pair g vector mapped to origin
};


struct ST_RESULT_DATA
{
	// a struct for exporting final data
	int MET_RES;
	int ARO_RES;
	float NORM_V;
	float MET_THETA;
	float MET_PHI;
};


/* ------------------- */
/* function prototypes */
/* ------------------- */
template <class ANG1>
ANG1 deg2rad(ANG1 deg);
template <class ANG2>
ANG2 rad2deg(ANG2 rad);
std::vector<std::string> split(std::string str);
template <class T1>
T1 norm_vector(std::vector<T1> u);
template <class T2>
std::vector<T2> unit_vector(std::vector<T2> u);
template <class T3>
T3 dot_vector(std::vector <T3> u, std::vector <T3> v);
template <class T4>
T4 angle_between_vectors(std::vector <T4> u, std::vector <T4> v);
template <class T5>
std::vector<std::vector<T5> > mmult(std::vector<std::vector<T5> > A, std::vector<std::vector<T5> > B);
template <class T6>
std::vector<std::vector<T6> > add_matrices(std::vector<std::vector<T6> > M1, std::vector<std::vector<T6> > M2);
rodrigues_rotation_output rodrigues_rotation(std::vector<float> v1, std::vector<float> v2, std::vector<float> v3);
std::vector < std::vector < float > > hexagon_midpoints(std::vector < std::vector < float > > b);


std::tuple<int, std::vector<ST_RESULT_DATA>, double> met_aromatic(std::string pdb_code){
	/*
	Parameters:
		pdb_code                    -> "1rcy.pdb"
	Returns: an std::tuple of data
		int                         -> some error/success code {0 : success, -1 : no MET/ARO, -2 : no interactions}
		std::vector<ST_RESULT_DATA> -> the output data (non-empty if successful acquisition)
		double                      -> execution time in seconds
	*/
	
	/* ----------- */
	/* begin timer */
	/* ----------- */
	std::clock_t t_start = std::clock();
    double t_execution;


	/* ------------------ */
	/* get all input data */
	/* ------------------ */
	std::fstream pdb_input;
	std::string line;
	std::vector <std::string> raw_input;

	pdb_input.open(pdb_code);  // pass path through main? raw path doesn't lead to result
	if(pdb_input.is_open()) 
	{
		while(getline(pdb_input, line)) 
		{
			raw_input.push_back(line);
		}
	}
	pdb_input.close();


	/* ------------------- */
	/* split by whitespace */
	/* ------------------- */
	std::vector <std::vector <std::string> > S_I;

	for (unsigned int i = 0; i < raw_input.size(); ++i) 
	{
		S_I.push_back(split(raw_input[i]));
	}


 	/* ------------------------------------------------------------ */
	/* filter in ATOM records / A chains and route into Line struct */
	/* ------------------------------------------------------------ */
	std::vector <Line> ATOM_RECORDS;

	for (unsigned int i = 0; i < S_I.size(); ++i) 
	{
		if (S_I[i][0] == "ATOM" && S_I[i][4] == "A") 
		{
			Line t_stru_1 = {
				S_I[i][0],            // ATOM
			    stoi(S_I[i][1]),      // 1
			    S_I[i][2],            // N
			    S_I[i][3],            // THR
			    S_I[i][4],            // A
			    stoi(S_I[i][5]),      // 5
			    {
			    	stof(S_I[i][6]),  // x
			    	stof(S_I[i][7]),  // y
			    	stof(S_I[i][8])   // z
			    }
			};
			ATOM_RECORDS.push_back(t_stru_1);
		}

		// only parse first model in multi model entries
		else if(S_I[i][0] == "ENDMDL")
		{
			break;
		}
	}


	/* ----------------------------------------------------------- */
	/* push back specific residue data into appropriate containers */
	/* ----------------------------------------------------------- */
	std::regex regex_met     (PAT_MET);       // MET
    std::regex regex_phe_tyr (PAT_PHE_TYR);   // PHE, TYR
    std::regex regex_trp     (PAT_TRP);       // TRP
    std::vector <Line> D_M;  // methonine Line struct
    std::vector <Line> D_A;  // methonine Line struct

	for (unsigned int i = 0; i < ATOM_RECORDS.size(); ++i) 
	{
		// MET -> { CG / SD / CE }
		if (ATOM_RECORDS[i].RES == "MET" && 
			std::regex_match(ATOM_RECORDS[i].ATOM_ID, regex_met))
		{
			D_M.push_back(ATOM_RECORDS[i]);
		}

		// PHE -> { CG / CD2 / CE2 / CZ / CE1 / CD1 }
		else if (ATOM_RECORDS[i].RES == "PHE" && 
			     std::regex_match(ATOM_RECORDS[i].ATOM_ID, regex_phe_tyr))
		{
			D_A.push_back(ATOM_RECORDS[i]);
		}

		// TYR -> { CG / CD2 / CE2 / CZ / CE1 / CD1 }
		else if (ATOM_RECORDS[i].RES == "TYR" && 
			     std::regex_match(ATOM_RECORDS[i].ATOM_ID, regex_phe_tyr))
		{
			D_A.push_back(ATOM_RECORDS[i]);
		}

		// TRP -> { CD2 / CE3 / CZ3 / CH2 / CZ2 / CE2 }
		else if (ATOM_RECORDS[i].RES == "TRP" && 
			     std::regex_match(ATOM_RECORDS[i].ATOM_ID, regex_trp))
		{
			D_A.push_back(ATOM_RECORDS[i]);
		}
	}


	/* ----------------------------------- */
	/* exit function call if no METs or AROs
	/* ----------------------------------- */
	if (!D_M.size() || !D_A.size()) 
	{
		std::vector <ST_RESULT_DATA> RESULTS = {};
		t_execution = (std::clock() - t_start) / (double) CLOCKS_PER_SEC;
		return std::make_tuple(NO_RESIDUE_ERR, RESULTS, t_execution);
	}
	

	/* ---------------------------------- */
	/* swap ATOM IDs with ascii_uppercase */
	/* ---------------------------------- */
	std::map<std::string, std::string> MAP_ATOMS_PHE;
	MAP_ATOMS_PHE["CG"]  = std::string("A");
	MAP_ATOMS_PHE["CD2"] = std::string("B");
	MAP_ATOMS_PHE["CE2"] = std::string("C");
	MAP_ATOMS_PHE["CZ"]  = std::string("D");
	MAP_ATOMS_PHE["CE1"] = std::string("E");
	MAP_ATOMS_PHE["CD1"] = std::string("F");

	std::map<std::string, std::string> MAP_ATOMS_TYR;
	MAP_ATOMS_TYR["CG"]  = std::string("A");
	MAP_ATOMS_TYR["CD2"] = std::string("B");
	MAP_ATOMS_TYR["CE2"] = std::string("C");
	MAP_ATOMS_TYR["CZ"]  = std::string("D");
	MAP_ATOMS_TYR["CE1"] = std::string("E");
	MAP_ATOMS_TYR["CD1"] = std::string("F");

	std::map<std::string, std::string> MAP_ATOMS_TRP;
	MAP_ATOMS_TRP["CD2"] = std::string("A");
	MAP_ATOMS_TRP["CE3"] = std::string("B");
	MAP_ATOMS_TRP["CZ3"] = std::string("C");
	MAP_ATOMS_TRP["CH2"] = std::string("D");
	MAP_ATOMS_TRP["CZ2"] = std::string("E");
	MAP_ATOMS_TRP["CE2"] = std::string("F");

	std::vector <int> to_dd; 
	for (unsigned int i = 0; i < D_A.size(); ++i)
	{
		to_dd.push_back(D_A[i].RES_POS);
		if (D_A[i].RES == "PHE")
		{
			D_A[i].ATOM_ID = MAP_ATOMS_PHE[D_A[i].ATOM_ID];
		}
		else if (D_A[i].RES == "TYR")
		{
			D_A[i].ATOM_ID = MAP_ATOMS_TYR[D_A[i].ATOM_ID];
		}
		else
		{
			D_A[i].ATOM_ID = MAP_ATOMS_TRP[D_A[i].ATOM_ID];
		}
	}


	/* ------------------------- */
	/* get all hexagon midpoints */
	/* ------------------------- */
	std::set <int> drop_duplicates(to_dd.begin(), to_dd.end()); // get sorted array of residue positions
	std::vector <int> residues(drop_duplicates.begin(), drop_duplicates.end()); // cast back to vector
	std::vector <Line> MIDPOINTS;

	for (std::vector <int>::iterator r = residues.begin(); r != residues.end(); ++r)
	{
		std::vector <Line> gb;  /* groupby object needs to be destroyed every iteration */
		std::vector < std::vector <float> > to_mp;  /* data to midpoint calculation */
		std::vector < std::vector <float> > mp;  /* the midpoint coordinates */

		/* groupby residue position number in residues */
		for (unsigned int i = 0; i < D_A.size(); ++i)
		{
			if (D_A[i].RES_POS == *r) { gb.push_back(D_A[i]); }
		}

		/* then sort in groupby object */
		std::sort (std::begin(gb), std::end(gb), [&](const Line& first, const Line& second) -> bool
		{
			return (first.ATOM_ID[0] < second.ATOM_ID[0]); 
		});

		/* prepare vector of coordinates for midpoint calculation */
		for (unsigned int m = 0; m < gb.size(); ++m) { to_mp.push_back(gb[m].v); }
		
		/* get the midpoints */
		mp = hexagon_midpoints(to_mp);

		/* load all midpoints into MIDPOINTS struct */
		/* the midpoints have been compared against data from Python MetAromatic */
		for (unsigned int i = 0; i < gb.size(); ++i)
		{
			Line t_stru_2 = 
			{
				gb[i].ATOM,
				gb[i].ATOM_NUMBER,
				gb[i].ATOM_ID,
				gb[i].RES,
				gb[i].CHAIN,
				gb[i].RES_POS,
				{mp[i][0], mp[i][1], mp[i][2]}
			};
			MIDPOINTS.push_back(t_stru_2);
		}
	}
	

	// ----------------------------------------
	// get all MET lone pair vectors in advance
	// ----------------------------------------
	MET_LP_STRUCT MET_LP;
	std::vector <MET_LP_STRUCT> met_lp;
	rodrigues_rotation_output R_O;          // adapter between rodrigues function and met_lp

	for (unsigned int i = 0; i < D_M.size(); ++i)
	{
		/*
			Atoms should always be ordered CG, SD, then CE
			according to PDB formatting standards. We can
			take advantage of this and just index at the
			i - 1, i and i + 1 D_M[].ATOM_D indices to get 
			the correctly ordered coordinate data.
		*/
		if (D_M[i].ATOM_ID == "SD")
		{
			R_O = rodrigues_rotation(D_M[i].v, D_M[i + 1].v, D_M[i - 1].v);
			MET_LP.RES_POS = D_M[i].RES_POS;
			MET_LP.v = D_M[i].v;
			MET_LP.lone_pair_a = R_O.vector_a;
			MET_LP.lone_pair_g = R_O.vector_g;
			met_lp.push_back(MET_LP);
		}
	}

	
	// ----------------------------------------------
	// find closely spaced methionine / aromatic pairs
	// then compute lone pair / vector v angles
	// ----------------------------------------------
	float norm_vector_v;
	float met_theta;
	float met_phi;

	std::vector <ST_RESULT_DATA> RESULTS;

	for (unsigned int i = 0; i < met_lp.size(); ++i)
	{
		for (unsigned int j = 0; j < MIDPOINTS.size(); ++j)
		{
			std::vector <float> vector_v;  // has to be constructed and destructed in the loop
			for (unsigned int k = 0; k < 3; ++k) { vector_v.push_back(MIDPOINTS[j].v[k] - met_lp[i].v[k]); };
			norm_vector_v = norm_vector(vector_v);
			if (norm_vector_v < DIST_CUTOFF)  // the distance condition -> see DSW thesis
			{
				met_theta = rad2deg(angle_between_vectors(vector_v, met_lp[i].lone_pair_a));
				met_phi = rad2deg(angle_between_vectors(vector_v, met_lp[i].lone_pair_g));

				if (met_theta < ANGLE_CUTOFF || met_phi < ANGLE_CUTOFF)  // the angular condition -> see DSW thesis
				{
					ST_RESULT_DATA t_stru_3 = 
					{
						met_lp[i].RES_POS,
						MIDPOINTS[j].RES_POS,
						norm_vector_v,
						met_theta,
						met_phi
					};
					RESULTS.push_back(t_stru_3);
				}
			}
		}
	}


	/* ----------- */
	/* export data */
	/* ----------- */
	t_execution = (std::clock() - t_start) / (double) CLOCKS_PER_SEC;
	int retval;
	RESULTS.size() == 0 ? retval = NO_RESULT_ERR : retval = SUCCESS_CODE;
	return std::make_tuple(retval, RESULTS, t_execution);
}


/* --------------- */
/* other functions */
/* --------------- */

template <class ANG1>
ANG1 deg2rad(ANG1 deg) 
{ 
	// convert to radians
	return (deg / HC) * PI; 
}


template <class ANG2>
ANG2 rad2deg(ANG2 rad) 
{
	// convert to degrees
	return (rad / PI) * HC; 
}


std::vector<std::string> split(std::string str) 
{
	// splits a string by whitespace
	std::vector<std::string> result;
	std::istringstream iss(str);
	for (std::string s; iss >> s; ) {
		result.push_back(s);
	}
	return result;
}


template <class T1>
T1 norm_vector(std::vector<T1> u)
{
	// gets norm of vector
	T1 norm = 0.0;
	for(int i = 0; i < u.size(); ++i)
	{
		norm += pow(u[i], 2);
	}
	return sqrt(norm);
}


template <class T2>
std::vector<T2> unit_vector(std::vector<T2> u)
{
	// gets unit vector
	T2 norm = norm_vector(u);
	std::vector<T2> unit;
	for(int i = 0; i < u.size(); ++i) 
	{
		unit.push_back(u[i] / norm);
	}
	return unit;
}


template <class T3>
T3 dot_vector(std::vector <T3> u, std::vector <T3> v)
{
	// gets dot product of two vectors
	unsigned int it = u.size(); 
	T3 dot = 0.0;
	for(int i = 0; i < it; ++i)
	{
	    dot += u[i] * v[i];
	}
	  return dot;
}


template <class T4>
T4 angle_between_vectors(std::vector <T4> u, std::vector <T4> v)
{
	// finds angle between a set of vectors
  	T4 dot_vectors = dot_vector(u, v);
  	T4 norm_u = norm_vector(u);
  	T4 norm_v = norm_vector(v);
  	T4 r = dot_vectors / (norm_u * norm_v);
  	return acos(r);
}


template <class T5>
std::vector<std::vector<T5> > mmult(std::vector<std::vector<T5> > A, std::vector<std::vector<T5> > B) 
{
	// product of two matrices or a vector and a matrix
	// https://www.programiz.com/cpp-programming/examples/matrix-multiplication

	/*
	Examples:
		(1) Matrix by matrix:
			std::vector<std::vector<float> > M1 = {{10.0, 20.0}, 
	                                               {30.0, 40.0}};
			std::vector<std::vector<float> > M2 = {{4.00, 2.00}, 
	                                               {3.00, 4.00}};
	        std::vector<std::vector<float> > M3 = mmult(M1, M2);
	    (2) Matrix by vector:
	    	std::vector<std::vector<float> > M1 = {{10.0, 20.0}, 
	                                               {30.0, 40.0}};
			std::vector<std::vector<float> > M2 = {{4.00}, 
	                                               {3.00}};
	        std::vector<std::vector<float> > M3 = mmult(M1, M2);
	*/
	
	// it might be advisable to assign these values in main() 
	unsigned int rows1 = A.size();      // number of rows in first matrix
	unsigned int cols1 = A[0].size();   // number of columns in first matrix
	unsigned int cols2 = B[0].size();   // number of columns in second matrix

    // reserve space for output matrix
    // (m x n)(n x k) = m x k matrix
    std::vector<std::vector<T5> > matrix_output(rows1);
	for (unsigned int i = 0 ; i < rows1 ; i++) {
		matrix_output[i].resize(cols2);
	}

	// compute elements of output matrix
    for (unsigned int i = 0; i < rows1; i++) 
    {
        for (unsigned int j = 0; j < cols2; j++) 
        {
            T5 init = 0;
            for (unsigned int k = 0; k < cols1; k++)
            {
                init += A[i][k] * B[k][j];
            }
            matrix_output[i][j] = init;
        }
    }
	
	return matrix_output;
}


template <class T6>
std::vector<std::vector<T6> > add_matrices(std::vector<std::vector<T6> > M1, std::vector<std::vector<T6> > M2)
{
	// a function for adding two m x n matrices
	// it is assumed that both matrices have dimensions m x n
	std::vector<std::vector<T6> > r;
	for (int row = 0; row < M1.size(); ++row)
	{
		std::vector<T6> s;
		for (int col = 0; col < M1[0].size(); ++col)
		{
			s.push_back(M1[row][col] + M2[row][col]);
		}
		r.push_back(s);
	}
	return r;
}


rodrigues_rotation_output rodrigues_rotation(std::vector<float> v1, std::vector<float> v2, std::vector<float> v3)
{
	/* ------------------------------------------- */
	/*
	This function uses Rodrigues' rotation method to rotate a set of vectors in 3-dimensional
	space about their mean vector.
	v1 -> The vector which describes the origin of both vectors v2 and v3
	   -> This is the SD coordinate in this case
	v2 -> The vector CE or CG
	v3 -> The vector CE or CG
	*/
	/* ------------------------------------------- */

	std::vector<float> a;
	std::vector<float> g;
	std::vector<float> r;

	// map to origin then flip
	for (unsigned int i = 0; i < v1.size(); ++i)
	{
		a.push_back((-1) * (v2[i] - v1[i]));
		g.push_back((-1) * (v3[i] - v1[i]));
	}

	// find axis of rotation
	for (unsigned int i = 0; i < 3; ++i) { r.push_back(0.5 * (a[i] + g[i])); }

	// overwrite into unit axis of rotation
	r = unit_vector(r);

	// yield the "cross product matrix" K from r
	std::vector<std::vector<float> > K = 
	{
		{0.0,  -r[2],  r[1]},
		{r[2],   0.0, -r[0]},
		{-r[1], r[0],   0.0}
	};

	// the identity matrix I
	std::vector<std::vector<float> > I = 
	{
		{1.0, 0.0, 0.0},
		{0.0, 1.0, 0.0},
		{0.0, 0.0, 1.0}
	};

	// compute the rotation matrix R from I and K
	std::vector<std::vector<float> > R_A = add_matrices(I, K);
	std::vector<std::vector<float> > R_B = mmult(K, K);
	std::vector<std::vector<float> > R = add_matrices(R_A, R_B);

	// rotate the C - S bond vectors to align with lone pairs
	std::vector<std::vector<float> > a_p = mmult(R, { {a[0]}, {a[1]}, {a[2]} });
	std::vector<std::vector<float> > g_p = mmult(R, { {g[0]}, {g[1]}, {g[2]} });

	// set up output
	rodrigues_rotation_output vectors;
	vectors.vector_a = {a_p[0][0], a_p[1][0], a_p[2][0]};
	vectors.vector_g = {g_p[0][0], g_p[1][0], g_p[2][0]};
	return vectors;
}


std::vector < std::vector < float > > hexagon_midpoints(std::vector < std::vector < float > > coords)
{
	/*
	A function for getting hexagon midpoints using row offsetting.
	
	Input
	-----
	std::vector < std::vector < float > > [] = 
	{
		{1.0, 1.0, 1.0},
		{2.0, 2.0, 2.0},
		{3.0, 3.0, 3.0},
		{4.0, 4.0, 4.0},
		{5.0, 5.0, 5.0},
		{6.0, 6.0, 6.0}
	};

	Output
	------
	{
		{1.5, 1.5, 1.5},
		{2.5, 2.5, 2.5},
		{3.5, 3.5, 3.5},
		{4.5, 4.5, 4.5},
		{5.5, 5.5, 5.5},
		{3.5, 3.5, 3.5}
	}
	*/
	std::vector < std::vector < float > > outer; 
	for (int i = 0; i < 6; ++i)
	{
		std::vector <float> inner;
		switch (i)
		{
			default: for (int j = 0; j < 3; ++j) 
			{ 
				inner.push_back(0.5 * (coords[i][j] + coords[i + 1][j])); 
			}
			case 5:  for (int j = 0; j < 3; ++j) 
			{ 
				inner.push_back(0.5 * (coords[5][j] + coords[0][j])); 
			}
		}
		outer.push_back(inner);
	}
	return outer;
}

#endif