#ifndef __DRUSILLA_SELECT_H
#define __DRUSILLA_SELECT_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>

#include "def.h"
#include "util.h"
#include "pri_queue.h"

class MaxK_List;

// -----------------------------------------------------------------------------
//  Drusilla_Select: data structure of Drusilla_Select for c-AFN search
// -----------------------------------------------------------------------------
class Drusilla_Select {
public:
	Drusilla_Select(				// default constructor
		int   n,						// number of data objects
		int   d,						// number of dimensions
		int   l,						// number of projections
		int   m,						// number of candidates on each proj
		const float **data);			// data objects
	
	// -------------------------------------------------------------------------
	~Drusilla_Select();				// destrcutor

	// -------------------------------------------------------------------------
	void display();					// display parameters	

	// -------------------------------------------------------------------------
	int kfn(						// c-k-AFN search via Drusilla Select
		const float *query,				// query point
		MaxK_List *list);				// top-k results (return)

protected:
	int   n_pts_;					// number of data objects
	int   dim_;						// dimensionality
	int   l_;						// number of projections
	int   m_;						// number of candidates for each proj
	int   *cand_;					// furthest neighbor candidates	
	const float **data_;			// data objects

	// -------------------------------------------------------------------------
	void bulkload();				// build hash tables

	// -------------------------------------------------------------------------
	void calc_shift_data(			// calculate shift data objects
		int   &max_id,					// data id with max l2-norm (return)
		float &max_norm,				// max l2-norm (return)
		float *norm,					// l2-norm of shift data (return)
		float **shift_data); 			// shift data (return)
};

#endif // __DRUSILLA_SELECT_H