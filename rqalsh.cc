#include "rqalsh.h"

// -----------------------------------------------------------------------------
RQALSH::RQALSH(						// constructor
	int   n,							// cardinality
	int   d,							// dimensionality
	float ratio,						// approximation ratio
	const int *index,					// index of data objects
	const float *data)					// data objects
	: n_pts_(n), dim_(d), ratio_(ratio), index_(index), data_(data)
{
	if (n <= N_THRESHOLD) {
		w_      = 0.0f;
		m_      = 0;
		l_      = 0;
		proj_a_ = NULL;
		tables_ = NULL;
	}
	else {
		// auto tuning w and determine m and l
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

		// generate hash functions
		proj_a_ = new float[m_ * d];
		for (int i = 0; i < m_ * d; ++i) proj_a_[i] = gaussian(0.0f, 1.0f);
		
		// build hash tables
		int id = -1;
		tables_ = new Result[m_ * n];
		for (int i = 0; i < m_; ++i) {
			for (int j = 0; j < n; ++j) {
				if (index_) id = index_[j];
				else id = j;

				tables_[i*n+j].id_  = j;
				tables_[i*n+j].key_ = calc_hash_value(i, &data_[id*d]);
			}
			qsort(&tables_[i*n], n, sizeof(Result), ResultComp);
		}
	}
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
	return calc_inner_product(dim_, &proj_a_[tid*dim_], data);
}

// -----------------------------------------------------------------------------
RQALSH::~RQALSH()					// destructor
{
	if (proj_a_ != NULL) { delete[] proj_a_; proj_a_ = NULL; }
	if (tables_ != NULL) { delete[] tables_; tables_ = NULL; }
}

// -------------------------------------------------------------------------
void RQALSH::display()				// display the parameters
{
	printf("Parameters of RQALSH:\n");
	printf("    n     = %d\n",   n_pts_);
	printf("    d     = %d\n",   dim_);
	printf("    ratio = %.1f\n", ratio_);
	printf("    w     = %f\n",   w_);
	printf("    m     = %d\n",   m_);
	printf("    l     = %d\n\n", l_);
}

// -----------------------------------------------------------------------------
int RQALSH::kfn(					// k-FN search (for ML_RQALSH)
	int   top_k,						// top-k value
	float R,							// limited search range
	const float *query,					// input query
	MaxK_List *list)					// c-k-AFN results (return)
{
	if (n_pts_ <= N_THRESHOLD) {
		int   id   = -1;
		float dist = -1.0f;
		for (int i = 0; i < n_pts_; ++i) {
			if (index_) id = index_[i];
			else id = i;

			dist = calc_l2_dist(dim_, query, &data_[id*dim_]);
			list->insert(dist, id + 1);
		}
		return n_pts_;
	}

	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	int   *freq    = new int[n_pts_];	// separation counting frequency
	int   *l_pos   = new int[m_];		// left  position of query
	int   *r_pos   = new int[m_];		// right position of query
	bool  *checked = new bool[n_pts_];	// is check for data id?
	bool  *b_flag  = new bool[m_];		// bucket flag
	bool  *r_flag  = new bool[m_];		// range  flag
	float *q_val   = new float[m_];		// hash value of query

	memset(freq,    0,     n_pts_ * SIZEFLOAT);
	memset(checked, false, n_pts_ * SIZEBOOL);
	memset(b_flag,  true,  m_     * SIZEBOOL);
	memset(r_flag,  true,  m_     * SIZEBOOL);

	for (int i = 0; i < m_; ++i) {
		q_val[i] = calc_hash_value(i, query);
		l_pos[i] = 0;  
		r_pos[i] = n_pts_ - 1;
	}

	// -------------------------------------------------------------------------
	//  c-k-AFN search
	// -------------------------------------------------------------------------
	int   cand      = CANDIDATES + top_k - 1; // candidate size
	int   cand_cnt  = 0;			// candidate counter
	int   num_range = 0; 			// number of search range flag
	float radius    = find_radius(l_pos, r_pos, q_val); 	// search radius
	float width     = radius * w_ / 2.0f; 					// bucket width
	float range     = R < CHECK_ERROR ? 0.0f : R*w_/2.0f; 	// search range

	while (true) {
		// ---------------------------------------------------------------------
		//  step 1: initialization
		// ---------------------------------------------------------------------
		int num_bucket = 0;
		memset(b_flag, true, m_ * SIZEBOOL);

		// ---------------------------------------------------------------------
		//  step 2: (R,c)-FN search
		// ---------------------------------------------------------------------
		while (num_bucket < m_ && num_range < m_) {
			for (int j = 0; j < m_; ++j) {
				// CANNOT add !r_flag[j] as condition, because the
				// r_flag[j] for large radius will affect small radius
				if (!b_flag[j]) continue;

				int    cnt = -1, lpos = -1, rpos = -1;
				float  q_v = q_val[j], ldist = -1.0f, rdist = -1.0f;
				Result *table = &tables_[j * n_pts_];

				// -------------------------------------------------------------
				//  step 2.1: scan left part of hash table
				// -------------------------------------------------------------
				cnt = 0;
				lpos = l_pos[j]; rpos = r_pos[j];
				while (cnt < SCAN_SIZE) {
					ldist = MINREAL;
					if (lpos < rpos) ldist = fabs(q_v - table[lpos].key_);
					else break;
					if (ldist < width || ldist < range) break;

					int id = table[lpos].id_;
					if (++freq[id] >= l_ && !checked[id]) {
						checked[id] = true;

						int did = id;
						if (index_ != NULL) did = index_[id];
						// printf("did = %d\n", did);
						float dist = calc_l2_dist(dim_, query, &data_[did*dim_]);
						list->insert(dist, did + 1);
						if (++cand_cnt >= cand) break;
					}
					++lpos; ++cnt;
				}
				if (cand_cnt >= cand) break;
				l_pos[j] = lpos;

				// -------------------------------------------------------------
				//  step 2.2: scan right part of hash table
				// -------------------------------------------------------------
				cnt = 0;
				while (cnt < SCAN_SIZE) {
					rdist = MINREAL;
					if (lpos < rpos) rdist = fabs(q_v - table[rpos].key_);
					else break;
					if (rdist < width || rdist < range) break;

					int id = table[rpos].id_;
					if (++freq[id] >= l_ && !checked[id]) {
						checked[id] = true;
						
						int did = id;
						if (index_ != NULL) did = index_[id];
						// printf("did = %d\n", did);
						float dist = calc_l2_dist(dim_, query, &data_[did*dim_]);
						list->insert(dist, did + 1);
						if (++cand_cnt >= cand) break;
					}
					--rpos; ++cnt;
				}
				if (cand_cnt >= cand) break;
				r_pos[j] = rpos;

				// -------------------------------------------------------------
				//  step 2.3: check whether this width is finished scanned
				// -------------------------------------------------------------
				if (lpos >= rpos || (ldist < width && rdist < width)) {
					if (b_flag[j]) { b_flag[j] = false; ++num_bucket; }
				}
				if (lpos >= rpos || (ldist < range && rdist < range)) {
					if (b_flag[j]) { b_flag[j] = false; ++num_bucket; }
					if (r_flag[j]) { r_flag[j] = false; ++num_range;  }
				}
				// use break after checking both b_flag and r_flag
				if (num_bucket >= m_ || num_range >= m_) break;
			}
			if (num_bucket >= m_ || num_range >= m_) break;
			if (cand_cnt >= cand) break;
		}
		// ---------------------------------------------------------------------
		//  step 3: stop condition
		// ---------------------------------------------------------------------
		if (num_range >= m_ || cand_cnt >= cand) break;

		// ---------------------------------------------------------------------
		//  step 4: update radius
		// ---------------------------------------------------------------------
		radius = radius / ratio_;
		width  = radius * w_ / 2.0f;
	}
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] freq;
	delete[] l_pos;
	delete[] r_pos;
	delete[] checked;
	delete[] b_flag;
	delete[] r_flag;
	delete[] q_val;

	return cand_cnt;
}

// -----------------------------------------------------------------------------
float RQALSH::find_radius(			// find proper radius
	const int   *l_pos,					// left  position of query in hash table
	const int   *r_pos,					// right position of query in hash table
	const float *q_val)					// hash value of query
{
	// find projected distance closest to the query in each hash tables 
	std::vector<float> list;
	for (int i = 0; i < m_; ++i) {
		if (l_pos[i] < r_pos[i]) {
			list.push_back(fabs(tables_[i*n_pts_+l_pos[i]].key_ - q_val[i]));
			list.push_back(fabs(tables_[i*n_pts_+r_pos[i]].key_ - q_val[i]));
		}
	}
	// sort the array in ascending order 
	std::sort(list.begin(), list.end());

	// find the median distance and return the new radius
	int   num  = (int) list.size();
	float dist = -1.0f;
	if (num % 2 == 0) dist = (list[num / 2 - 1] + list[num / 2]) / 2.0f;
	else dist = list[num / 2];

	int kappa = (int) ceil(log(2.0f * dist / w_) / log(ratio_));
	return pow(ratio_, kappa);
}
