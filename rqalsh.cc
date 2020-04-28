#include "rqalsh.h"

// -----------------------------------------------------------------------------
RQALSH::RQALSH(						// constructor
	int   n,							// cardinality
	int   d,							// dimensionality
	float ratio)						// approximation ratio
{
	n_pts_ = n;
	dim_   = d;
	ratio_ = ratio;
	data_  = NULL;
	
	init();
}

// -----------------------------------------------------------------------------
RQALSH::RQALSH(						// constructor
	int   n,							// cardinality
	int   d,							// dimensionality
	float ratio,						// approximation ratio
	const float **data)					// data objects
{
	n_pts_ = n;
	dim_   = d;
	ratio_ = ratio;
	data_  = data;
	
	init();

	// -------------------------------------------------------------------------
	//  build hash tables
	// -------------------------------------------------------------------------
	for (int i = 0; i < m_; ++i) {
		for (int j = 0; j < n_pts_; ++j) {
			tables_[i][j].id_ = j;
			tables_[i][j].key_ = calc_hash_value(i, data_[j]);
		}
		qsort(tables_[i], n_pts_, sizeof(Result), ResultComp);
	}
}

// -----------------------------------------------------------------------------
void RQALSH::init()					// init parameters
{
	// -------------------------------------------------------------------------
	//  init <w_> <m_> and <l_> (auto tuning-w)
	// -------------------------------------------------------------------------
	w_  = sqrt((8.0f * log(ratio_)) / (ratio_ * ratio_ - 1.0f));
	
	float p1 = calc_l2_prob(w_ / 2.0f);
	float p2 = calc_l2_prob(w_ * ratio_ / 2.0f);

	float beta  = (float) CANDIDATES / (float) n_pts_;
	float delta = 0.49f;

	float para1 = sqrt(log(2.0f / beta));
	float para2 = sqrt(log(1.0f / delta));
	float para3 = 2.0f * (p1 - p2) * (p1 - p2);

	float eta   = para1 / para2;
	float alpha = (eta * p1 + p2) / (1.0f + eta);

	m_ = (int) ceil((para1 + para2) * (para1 + para2) / para3);
	l_ = (int) ceil(alpha * m_);

	// -------------------------------------------------------------------------
	//  calc parameters and generate hash functions
	// -------------------------------------------------------------------------
	g_memory += SIZEFLOAT * m_ * dim_;
	a_ = new float*[m_];
	for (int i = 0; i < m_; ++i) {
		a_[i] = new float[dim_];
		for (int j = 0; j < dim_; ++j) {
			a_[i][j] = gaussian(0.0f, 1.0f);
		}
	}
	
	// -------------------------------------------------------------------------
	//  bulkloading
	// -------------------------------------------------------------------------
	g_memory += sizeof(Result) * m_ * n_pts_;
	tables_ = new Result*[m_];
	for (int i = 0; i < m_; ++i) tables_[i] = new Result[n_pts_];
}

// -------------------------------------------------------------------------
inline float RQALSH::calc_l2_prob(	// calc <p1> and <p2> for L2 distance
	float x)							// x = w / (2.0 * r)
{
	return 1.0f - new_gaussian_prob(x);
}

// -------------------------------------------------------------------------
float RQALSH::calc_hash_value( 		// calc hash value
	int   tid,							// hash table id
	const float *data)					// one data object
{
	return calc_inner_product(dim_, a_[tid], data);
}

// -----------------------------------------------------------------------------
RQALSH::~RQALSH()					// destructor
{
	for (int i = 0; i < m_; ++i) {
		delete[] a_[i]; a_[i] = NULL;
		delete[] tables_[i]; tables_[i] = NULL;
	}
	delete[] a_; a_ = NULL;
	delete[] tables_; tables_ = NULL;

	g_memory -= SIZEFLOAT * m_ * dim_;
	g_memory -= sizeof(Result) * m_ * n_pts_;
}

// -------------------------------------------------------------------------
void RQALSH::display()				// display the parameters
{
	printf("Parameters of RQALSH:\n");
	printf("    n          = %d\n", n_pts_);
	printf("    d          = %d\n", dim_);
	printf("    ratio      = %f\n", ratio_);
	printf("    w          = %f\n", w_);
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
	int   *freq      = new int[n_pts_];
	int   *left_pos  = new int[m_];
	int   *right_pos = new int[m_];
	bool  *checked   = new bool[n_pts_];
	bool  *flag      = new bool[m_];
	float *q_val     = new float[m_];
	
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	memset(freq, 0, n_pts_ * SIZEFLOAT);
	memset(checked, false, n_pts_ * SIZEBOOL);

	for (int i = 0; i < m_; ++i) {
		q_val[i]     = calc_hash_value(i, query);
		left_pos[i]  = 0;  
		right_pos[i] = n_pts_ - 1;
	}

	// -------------------------------------------------------------------------
	//  c-k-AFN search
	// -------------------------------------------------------------------------
	int candidates = CANDIDATES + top_k - 1; // candidate size
	int cand_cnt   = 0;				// candidate counter

	float kdist  = MINREAL;
	float radius = find_radius(left_pos, right_pos, q_val); // search radius
	float width  = radius * w_ / 2.0f; // bucket width

	while (true) {
		// ---------------------------------------------------------------------
		//  step 1: initialization
		// ---------------------------------------------------------------------
		int num_flag = 0;
		memset(flag, true, m_ * SIZEBOOL);

		// ---------------------------------------------------------------------
		//  step 2: (R,c)-FN search
		// ---------------------------------------------------------------------
		while (num_flag < m_) {
			for (int j = 0; j < m_; ++j) {
				if (!flag[j]) continue;

				int    cnt = -1, lpos = -1, rpos = -1;
				float  q_v = q_val[j], ldist = -1.0f, rdist = -1.0f;
				Result *table = tables_[j];
				
				// -------------------------------------------------------------
				//  step 2.1: scan the left part of bucket
				// -------------------------------------------------------------
				cnt = 0;
				lpos = left_pos[j]; rpos = right_pos[j];
				while (cnt < SCAN_SIZE) {
					ldist = MAXREAL;
					if (lpos < rpos) {
						ldist = fabs(q_v - table[lpos].key_);
					}
					else break;
					if (ldist < width) break;

					int id = table[lpos].id_;
					if (++freq[id] >= l_ && !checked[id]) {
						checked[id] = true;
						float dist = calc_l2_dist(dim_, data_[id], query);
						kdist = list->insert(dist, id + 1);

						if (++cand_cnt >= candidates) break;
					}
					++lpos; ++cnt;
				}
				if (cand_cnt >= candidates) break;
				left_pos[j] = lpos;

				// -------------------------------------------------------------
				//  step 2.2: scan the right part of bucket
				// -------------------------------------------------------------
				cnt = 0;
				while (cnt < SCAN_SIZE) {
					rdist = MAXREAL;
					if (lpos < rpos) {
						rdist = fabs(q_v - table[rpos].key_);
					}
					else break;
					if (rdist < width) break;

					int id = table[rpos].id_;
					if (++freq[id] >= l_ && !checked[id]) {
						checked[id] = true;
						float dist = calc_l2_dist(dim_, data_[id], query);
						kdist = list->insert(dist, id + 1);

						if (++cand_cnt >= candidates) break;
					}
					--rpos; ++cnt;
				}
				if (cand_cnt >= candidates) break;
				right_pos[j] = rpos;

				// -------------------------------------------------------------
				//  step 2.3: check whether this bucket is finished scanning
				// -------------------------------------------------------------
				if (lpos >= rpos || (ldist < width && rdist < width)) {
					flag[j] = false;
					++num_flag;
				}
				if (num_flag >= m_) break;
			}
			if (num_flag >= m_ || cand_cnt >= candidates) break;
		}

		// ---------------------------------------------------------------------
		//  step 3: stop conditions 1 & 2
		// ---------------------------------------------------------------------
		if (kdist > radius / ratio_ && cand_cnt >= top_k) break;
		if (cand_cnt >= candidates) break;

		// ---------------------------------------------------------------------
		//  step 4: update radius
		// ---------------------------------------------------------------------
		radius = radius / ratio_;
		width  = radius * w_ / 2.0f;
	}
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] freq;      freq      = NULL;
	delete[] left_pos;  left_pos  = NULL;
	delete[] right_pos; right_pos = NULL;
	delete[] checked;   checked   = NULL;
	delete[] flag;      flag      = NULL;
	delete[] q_val;     q_val     = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH::kfn(					// k-FN search (for ML_RQALSH)
	int   top_k,						// top-k value
	float R,							// limited search range
	const float *query,					// input query
	std::vector<int> &cand) 			// candidate id (return)
{
	int   *freq        = new int[n_pts_];
	int   *left_pos    = new int[m_];
	int   *right_pos   = new int[m_];
	bool  *checked     = new bool[n_pts_];
	bool  *bucket_flag = new bool[m_];
	bool  *range_flag  = new bool[m_];
	float *q_val       = new float[m_];

	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	memset(freq, 0, n_pts_ * SIZEFLOAT);
	memset(checked, false, n_pts_ * SIZEBOOL);
	memset(range_flag, true, m_ * SIZEBOOL);

	for (int i = 0; i < m_; ++i) {
		q_val[i]     = calc_hash_value(i, query);
		left_pos[i]  = 0;  
		right_pos[i] = n_pts_ - 1;
	}

	// -------------------------------------------------------------------------
	//  c-k-AFN search
	// -------------------------------------------------------------------------
	int candidates = CANDIDATES + top_k - 1; // candidate size
	int cand_cnt  = 0;				// candidate counter
	int num_range = 0; 				// number of search range flag

	float radius = find_radius(left_pos, right_pos, q_val); // search radius
	float width  = radius * w_ / 2.0f; // bucket width
	float range  = R < FLOATZERO ? 0.0f : R * w_ / 2.0f; // search range

	while (true) {
		// ---------------------------------------------------------------------
		//  step 1: initialization
		// ---------------------------------------------------------------------
		int num_bucket = 0;
		memset(bucket_flag, true, m_ * SIZEBOOL);

		// ---------------------------------------------------------------------
		//  step 2: (R,c)-FN search
		// ---------------------------------------------------------------------
		while (num_bucket < m_ && num_range < m_) {
			for (int j = 0; j < m_; ++j) {
				if (!bucket_flag[j]) continue;

				int    cnt = -1, lpos = -1, rpos = -1;
				float  q_v = q_val[j], ldist = -1.0f, rdist = -1.0f;
				Result *table = tables_[j];

				// -------------------------------------------------------------
				//  step 2.1: scan left part of hash table
				// -------------------------------------------------------------
				cnt = 0;
				lpos = left_pos[j]; rpos = right_pos[j];
				while (cnt < SCAN_SIZE) {
					ldist = MAXREAL;
					if (lpos < rpos) {
						ldist = fabs(q_v - table[lpos].key_);
					}
					else break;
					if (ldist < width || ldist < range) break;

					int id = table[lpos].id_;
					if (++freq[id] >= l_ && !checked[id]) {
						checked[id] = true;
						cand.push_back(id);

						if (++cand_cnt >= candidates) break;
					}
					++lpos; ++cnt;
				}
				if (cand_cnt >= candidates) break;
				left_pos[j] = lpos;

				// -------------------------------------------------------------
				//  step 2.2: scan right part of hash table
				// -------------------------------------------------------------
				cnt = 0;
				while (cnt < SCAN_SIZE) {
					rdist = MAXREAL;
					if (lpos < rpos) {
						rdist = fabs(q_v - table[rpos].key_);
					}
					else break;
					if (rdist < width || rdist < range) break;

					int id = table[rpos].id_;
					if (++freq[id] >= l_ && !checked[id]) {
						checked[id] = true;
						cand.push_back(id);

						if (++cand_cnt >= candidates) break;
					}
					--rpos; ++cnt;
				}
				if (cand_cnt >= candidates) break;
				right_pos[j] = rpos;

				// -------------------------------------------------------------
				//  step 2.3: check whether this width is finished scanned
				// -------------------------------------------------------------
				if (lpos >= rpos || (ldist < width && rdist < width)) {
					bucket_flag[j] = false;
					if (++num_bucket > m_) break;
				}
				if (lpos >= rpos || (ldist < range && rdist < range)) {
					if (bucket_flag[j]) {
						bucket_flag[j] = false;
						if (++num_bucket > m_) break;
					}
					if (range_flag[j]) {
						range_flag[j] = false;
						if (++num_range > m_) break;
					}
				}
			}
			if (num_bucket >= m_ || num_range >= m_) break;
			if (cand_cnt >= candidates) break;
		}
		// ---------------------------------------------------------------------
		//  step 3: stop condition
		// ---------------------------------------------------------------------
		if (num_range >= m_ || cand_cnt >= candidates) break;

		// ---------------------------------------------------------------------
		//  step 4: update radius
		// ---------------------------------------------------------------------
		radius = radius / ratio_;
		width  = radius * w_ / 2.0f;
	}
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] freq;        freq        = NULL;
	delete[] left_pos;    left_pos    = NULL;
	delete[] right_pos;   right_pos   = NULL;
	delete[] checked;     checked     = NULL;
	delete[] bucket_flag; bucket_flag = NULL;
	delete[] range_flag;  range_flag  = NULL;
	delete[] q_val;       q_val       = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
float RQALSH::find_radius(			// find proper radius
	const int *left_pos,				// left  position of query in hash table
	const int *right_pos,				// right position of query in hash table
	const float *q_val)					// hash value of query
{
	// -------------------------------------------------------------------------
	//  find an array of projected distance which is closest to the query in
	//  each of <m> hash tables 
	// -------------------------------------------------------------------------
	std::vector<float> list;
	for (int i = 0; i < m_; ++i) {
		if (left_pos[i] < right_pos[i]) {
			list.push_back(fabs(tables_[i][left_pos[i]].key_  - q_val[i]));
			list.push_back(fabs(tables_[i][right_pos[i]].key_ - q_val[i]));
		}
	}

	// -------------------------------------------------------------------------
	//  sort the array in ascending order 
	// -------------------------------------------------------------------------
	std::sort(list.begin(), list.end());

	// -------------------------------------------------------------------------
	//  find the median distance and return the new radius
	// -------------------------------------------------------------------------
	int num  = (int) list.size();
	float dist = -1.0f;
	if (num % 2 == 0) dist = (list[num / 2 - 1] + list[num / 2]) / 2.0f;
	else dist = list[num / 2];

	int kappa = (int) ceil(log(2.0f * dist / w_) / log(ratio_));
	return pow(ratio_, kappa);
}
