#include "headers.h"

// -----------------------------------------------------------------------------
RQALSH::RQALSH()					// default constructor
{
	n_pts_       = -1;
	dim_         = -1;
	beta_        = -1.0f;
	delta_       = -1.0f;
	appr_ratio_  = -1.0f;
	data_        = NULL;

	w_           = -1.0f;
	p1_          = -1.0f;
	p2_          = -1.0f;
	alpha_       = -1.0f;;
	m_           = -1;
	l_           = -1;
	a_array_     = NULL;
	tables_      = NULL;

	freq_        = NULL;
	lpos_        = NULL;
	rpos_        = NULL;
	checked_     = NULL;
	flag_        = NULL;
	q_val_       = NULL;
}

// -----------------------------------------------------------------------------
RQALSH::~RQALSH()					// destructor
{
	if (a_array_ != NULL) {
		delete[] a_array_; a_array_ = NULL;
	}
	if (tables_ != NULL) {
		for (int i = 0; i < m_; ++i) {
			delete[] tables_[i]; tables_[i] = NULL;
		}
		delete[] tables_; tables_ = NULL;
	}

	if (freq_ != NULL) {
		delete[] freq_; freq_ = NULL;
	}
	if (lpos_ != NULL) {
		delete[] lpos_; lpos_ = NULL;
	}
	if (rpos_ != NULL) {
		delete[] rpos_; rpos_ = NULL;
	}
	if (checked_ != NULL) {
		delete[] checked_; checked_ = NULL;
	}
	if (flag_ != NULL) {
		delete[] flag_; flag_ = NULL;
	}
	if (q_val_ != NULL) {
		delete[] q_val_; q_val_ = NULL;
	}
}

// -----------------------------------------------------------------------------
int RQALSH::build(					// build index
	int   n,							// cardinality
	int   d,							// dimensionality
	int   beta,							// false positive percentage
	float delta,						// error probability
	float ratio,						// approximation ratio
	const float **data)					// data objects
{
	// -------------------------------------------------------------------------
	//  init the input parameters
	// -------------------------------------------------------------------------
	n_pts_      = n;
	dim_        = d;
	beta_       = (float) beta / (float) n;
	delta_      = delta;
	appr_ratio_ = ratio;
	data_       = data;
	
	// -------------------------------------------------------------------------
	//  calc parameters and generate hash functions
	// -------------------------------------------------------------------------
	calc_params();
	gen_hash_func();
	
	// -------------------------------------------------------------------------
	//  bulkloading
	// -------------------------------------------------------------------------
	bulkload();
	display();
	
	return 0;
}

// -------------------------------------------------------------------------
void RQALSH::calc_params()			// calcobjects selection by largest norm
{
	// -------------------------------------------------------------------------
	//  init <w_> <p1_> and <p2_> (auto tuning-w)
	// -------------------------------------------------------------------------
	w_  = sqrt((8.0f*log(appr_ratio_)) / (appr_ratio_*appr_ratio_-1.0f));
	
	p1_ = calc_l2_prob(w_ / 2.0f);
	p2_ = calc_l2_prob(w_ * appr_ratio_ / 2.0f);

	// -------------------------------------------------------------------------
	//  init <alpha_> <m_> and <l_>
	// -------------------------------------------------------------------------
	float para1 = sqrt(log(2.0f / beta_));
	float para2 = sqrt(log(1.0f / delta_));
	float para3 = 2.0f * (p1_ - p2_) * (p1_ - p2_);

	float eta = para1 / para2;
	alpha_ = (eta * p1_ + p2_) / (1.0f + eta);

	m_ = (int) ceil((para1 + para2) * (para1 + para2) / para3);
	l_ = (int) ceil(alpha_ * m_);
	
	// -------------------------------------------------------------------------
	//  init parameters for c-k-AFN search
	// -------------------------------------------------------------------------
	freq_    = new int[n_pts_];
	lpos_    = new int[m_];
	rpos_    = new int[m_];
	checked_ = new bool[n_pts_];
	flag_    = new bool[m_];
	q_val_   = new float[m_];	
}

// -------------------------------------------------------------------------
float RQALSH::calc_l2_prob(			// calc <p1> and <p2> for L2 distance
	float x)							// x = w / (2.0 * r)
{
	return 1.0f - new_gaussian_prob(x);
}

// -------------------------------------------------------------------------
void RQALSH::gen_hash_func()		// generate hash functions
{
	int size = m_ * dim_;
	a_array_ = new float[size];
	for (int i = 0; i < size; i++) {
		a_array_[i] = gaussian(0.0f, 1.0f);
	}
}

// -------------------------------------------------------------------------
int RQALSH::bulkload()				// build hash tables for selected objs
{
	tables_ = new Result*[m_];
	for (int i = 0; i < m_; ++i) {
		tables_[i] = new Result[n_pts_];
		for (int j = 0; j < n_pts_; ++j) {
			tables_[i][j].id_ = j;
			tables_[i][j].key_ = calc_hash_value(i, data_[j]);
		}
		qsort(tables_[i], n_pts_, sizeof(Result), ResultComp);
	}

	return 0;
}

// -------------------------------------------------------------------------
float RQALSH::calc_hash_value(		// calc hash value
	int table_id,						// hash table id
	const float *data)					// one data object
{
	int base = table_id * dim_;
	float ret = 0.0f;
	for (int i = 0; i < dim_; ++i) {
		ret += (a_array_[base + i] * data[i]);
	}
	return ret;
}

// -------------------------------------------------------------------------
void RQALSH::display()				// display the parameters
{
	printf("Parameters of RQALSH:\n");
	printf("    n          = %d\n", n_pts_);
	printf("    d          = %d\n", dim_);
	printf("    ratio      = %f\n", appr_ratio_);
	printf("    w          = %f\n", w_);
	printf("    p1         = %f\n", p1_);
	printf("    p2         = %f\n", p2_);
	printf("    alpha      = %f\n", alpha_);
	printf("    beta       = %f\n", beta_);
	printf("    delta      = %f\n", delta_);
	printf("    m          = %d\n", m_);
	printf("    l          = %d\n", l_);
	printf("\n");
}

// -----------------------------------------------------------------------------
int RQALSH::kfn(					// k-FN search
	int top_k,							// top-k value
	const float *query,					// query object
	MaxK_List *list)					// k-FN results (return)
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	for (int i = 0; i < n_pts_; ++i) {
		freq_[i] = 0;
		checked_[i] = false;
	}

	for (int i = 0; i < m_; ++i) {
		q_val_[i] = calc_hash_value(i, query);
		lpos_[i]  = 0;  
		rpos_[i]  = n_pts_ - 1;
	}

	// -------------------------------------------------------------------------
	//  c-k-AFN search
	// -------------------------------------------------------------------------
	int   dist_cnt   = 0;			// counter of Eucldiean distance compuation
	int   candidates = CANDIDATES + top_k - 1; // candidate size
	int   num_flag   = -1;			// used for bucket width
	int   id         = -1;			// data object id

	float dist       = -1.0f;		// real distance between data and query
	float ldist      = -1.0f;		// left  projected distance with query
	float rdist      = -1.0f;		// right projected distance with query
	float kfn_dist   = MINREAL;		// k-FN real distance

	float radius = find_radius(q_val_, lpos_, rpos_); // search radius
	float bucket = radius * w_ / 2.0f; // bucket width

	while (true) {
		// ---------------------------------------------------------------------
		//  step 1: initialization
		// ---------------------------------------------------------------------
		num_flag = 0;
		for (int j = 0; j < m_; ++j) flag_[j] = true;

		// ---------------------------------------------------------------------
		//  step 2: (R, c)-FN search
		// ---------------------------------------------------------------------
		int count = 0;
		while (num_flag < m_) {
			for (int j = 0; j < m_; ++j) {
				if (!flag_[j]) continue;

				// -------------------------------------------------------------
				//  step 2.1: scan left part of hash table
				// -------------------------------------------------------------
				count = 0;
				while (count < SCAN_SIZE) {
					ldist = MAXREAL;
					if (lpos_[j] < rpos_[j]) {
						ldist = fabs(q_val_[j] - tables_[j][lpos_[j]].key_);
					}
					else break;
					if (ldist < bucket) break;

					int id = tables_[j][lpos_[j]].id_;
					if (++freq_[id] >= l_ && !checked_[id]) {
						checked_[id] = true;
						dist = calc_l2_dist(dim_, data_[id], query);
						kfn_dist = list->insert(dist, id + 1);

						if (++dist_cnt >= candidates) break;
					}
					lpos_[j]++;
					count++;
				}
				if (dist_cnt >= candidates) break;

				// -------------------------------------------------------------
				//  step 2.2: scan right part of hash table
				// -------------------------------------------------------------
				count = 0;
				while (count < SCAN_SIZE) {
					rdist = MAXREAL;
					if (lpos_[j] < rpos_[j]) {
						rdist = fabs(q_val_[j] - tables_[j][rpos_[j]].key_);
					}
					else break;
					if (rdist < bucket) break;

					int id = tables_[j][rpos_[j]].id_;
					if (++freq_[id] >= l_ && !checked_[id]) {
						checked_[id] = true;
						dist = calc_l2_dist(dim_, data_[id], query);
						kfn_dist = list->insert(dist, id + 1);

						if (++dist_cnt >= candidates) break;
					}
					rpos_[j]--;
					count++;
				}
				if (dist_cnt >= candidates) break;

				// -------------------------------------------------------------
				//  step 2.3: check whether this bucket is finished scanned
				// -------------------------------------------------------------
				if (lpos_[j] >= rpos_[j] || (ldist < bucket && rdist < bucket)) {
					flag_[j] = false;
					num_flag++;
				}
				if (num_flag >= m_) break;
			}
			if (num_flag >= m_ || dist_cnt >= candidates) break;
		}

		// ---------------------------------------------------------------------
		//  step 3: stop conditions 1 & 2
		// ---------------------------------------------------------------------
		if (kfn_dist > radius / appr_ratio_ && dist_cnt >= top_k) break;
		if (dist_cnt >= candidates) break;

		// ---------------------------------------------------------------------
		//  step 4: update radius
		// ---------------------------------------------------------------------
		radius = radius / appr_ratio_;
		bucket = radius * w_ / 2.0f;
	}
	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH::kfn(					// k-FN search
	int top_k,							// top-k value
	const float *query,					// query object
	const int *object_id,				// object id mapping
	MaxK_List *list)					// k-FN results (return)
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	for (int i = 0; i < n_pts_; ++i) {
		freq_[i] = 0;
		checked_[i] = false;
	}

	for (int i = 0; i < m_; ++i) {
		q_val_[i] = calc_hash_value(i, query);
		lpos_[i]  = 0;  
		rpos_[i]  = n_pts_ - 1;
	}

	// -------------------------------------------------------------------------
	//  c-k-AFN search
	// -------------------------------------------------------------------------
	int   dist_cnt   = 0;			// counter of Eucldiean distance compuation
	int   candidates = CANDIDATES + top_k - 1; // candidate size
	int   num_flag   = -1;			// used for bucket width
	int   id         = -1;			// current object id
	int   oid        = -1;			// actual data object id

	float dist       = -1.0f;		// real distance between data and query
	float ldist      = -1.0f;		// left  projected distance with query
	float rdist      = -1.0f;		// right projected distance with query
	float kfn_dist   = MINREAL;		// k-FN real distance

	float radius = find_radius(q_val_, lpos_, rpos_); // search radius
	float bucket = radius * w_ / 2.0f; // bucket width

	while (true) {
		// ---------------------------------------------------------------------
		//  step 1: initialization
		// ---------------------------------------------------------------------
		num_flag = 0;
		for (int j = 0; j < m_; ++j) flag_[j] = true;

		// ---------------------------------------------------------------------
		//  step 2: (R, c)-FN search
		// ---------------------------------------------------------------------
		int count = 0;
		while (num_flag < m_) {
			for (int j = 0; j < m_; ++j) {
				if (!flag_[j]) continue;

				// -------------------------------------------------------------
				//  step 2.1: scan left part of hash table
				// -------------------------------------------------------------
				count = 0;
				while (count < SCAN_SIZE) {
					ldist = MAXREAL;
					if (lpos_[j] < rpos_[j]) {
						ldist = fabs(q_val_[j] - tables_[j][lpos_[j]].key_);
					}
					else break;
					if (ldist < bucket) break;

					int id = tables_[j][lpos_[j]].id_;
					if (++freq_[id] >= l_ && !checked_[id]) {
						checked_[id] = true;
						dist = calc_l2_dist(dim_, data_[id], query);
						oid = object_id[id];
						kfn_dist = list->insert(dist, oid + 1);

						if (++dist_cnt >= candidates) break;
					}
					lpos_[j]++;
					count++;
				}
				if (dist_cnt >= candidates) break;

				// -------------------------------------------------------------
				//  step 2.2: scan right part of hash table
				// -------------------------------------------------------------
				count = 0;
				while (count < SCAN_SIZE) {
					rdist = MAXREAL;
					if (lpos_[j] < rpos_[j]) {
						rdist = fabs(q_val_[j] - tables_[j][rpos_[j]].key_);
					}
					else break;
					if (rdist < bucket) break;

					int id = tables_[j][rpos_[j]].id_;
					if (++freq_[id] >= l_ && !checked_[id]) {
						checked_[id] = true;
						dist = calc_l2_dist(dim_, data_[id], query);
						oid = object_id[id];
						kfn_dist = list->insert(dist, oid + 1);

						if (++dist_cnt >= candidates) break;
					}
					rpos_[j]--;
					count++;
				}
				if (dist_cnt >= candidates) break;

				// -------------------------------------------------------------
				//  step 2.3: check whether this bucket is finished scanned
				// -------------------------------------------------------------
				if (lpos_[j] >= rpos_[j] || (ldist < bucket && rdist < bucket)) {
					flag_[j] = false;
					num_flag++;
				}
				if (num_flag >= m_) break;
			}
			if (num_flag >= m_ || dist_cnt >= candidates) break;
		}

		// ---------------------------------------------------------------------
		//  step 3: stop conditions 1 & 2
		// ---------------------------------------------------------------------
		if (kfn_dist > radius / appr_ratio_ && dist_cnt >= top_k) break;
		if (dist_cnt >= candidates) break;

		// ---------------------------------------------------------------------
		//  step 4: update radius
		// ---------------------------------------------------------------------
		radius = radius / appr_ratio_;
		bucket = radius * w_ / 2.0f;
	}
	return 0;
}

// -----------------------------------------------------------------------------
float RQALSH::find_radius(			// find proper radius
	const float *q_val,					// hash value of query
	const int *lpos,					// left position of hash table
	const int *rpos)					// right position of hash table
{
	// -------------------------------------------------------------------------
	//  find an array of projected distance which is closest to the query in
	//  each of <m> hash tables 
	// -------------------------------------------------------------------------
	vector<float> list;
	for (int i = 0; i < m_; ++i) {
		if (lpos[i] < rpos[i]) {
			list.push_back(fabs(tables_[i][lpos[i]].key_ - q_val[i]));
			list.push_back(fabs(tables_[i][rpos[i]].key_ - q_val[i]));
		}
	}

	// -------------------------------------------------------------------------
	//  sort the array in ascending order 
	// -------------------------------------------------------------------------
	sort(list.begin(), list.end());

	// -------------------------------------------------------------------------
	//  find the median distance and return the new radius
	// -------------------------------------------------------------------------
	int num  = (int) list.size();
	float dist = -1.0f;
	if (num % 2 == 0) dist = (list[num / 2 - 1] + list[num / 2]) / 2.0f;
	else dist = list[num / 2];

	int kappa = (int) ceil(log(2.0f * dist / w_) / log(appr_ratio_));
	return pow(appr_ratio_, kappa);
}
