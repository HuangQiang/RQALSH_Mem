#include <algorithm>

#include "def.h"
#include "util.h"
#include "pri_queue.h"
#include "qdafn.h"

// -----------------------------------------------------------------------------
//  Ziggurat Method standard normal pseudorandom number generator code from 
//  George Marsaglia and Wai Wan Tsang (2000).
// 
//  "The Ziggurat Method for Generating Random Variables". Journal of
//  Statistical Software 5 (8).
// -----------------------------------------------------------------------------
long hz;
unsigned long jsr = 123456789;
unsigned long jz;
unsigned long iz, kn[128], ke[256];
float wn[128],fn[128], we[256],fe[256];

// -----------------------------------------------------------------------------
//  nfix() generates variates from the residue when rejection in RNOR occurs
// -----------------------------------------------------------------------------
float nfix()
{
	const float r = 3.442620;		// The start of the right tail
	static float x, y;

	while (1) {
		x = hz * wn[iz];
		if (iz == 0) {				// iz == 0, handles the base strip
			do {
				x = -log(UNI)*0.2904764;
				y = -log(UNI);
			} while (y+y < x*x);

			return (hz>0) ? (r+x) : (-r-x);
		}
									// iz > 0, handle the wedges of other strips
		if (fn[iz] + UNI*(fn[iz-1]-fn[iz]) < exp(-0.5*x*x)) {
			return x;
		}

		hz = SHR3;					// initiate, try to exit loop
		iz = hz & 127;
		if (abs(hz) < kn[iz]) return hz*wn[iz];
	}
}

// -----------------------------------------------------------------------------
//  efix() generates variates from the residue when rejection in REXP occurs
// -----------------------------------------------------------------------------
float efix()
{
	float x = -1.0f;
	while (1) {
		if (iz == 0) {
			return 7.69711-log(UNI);// iz == 0
		}

		x = jz * we[iz];
		if (fe[iz] + UNI*(fe[iz-1]-fe[iz]) < exp(-x)) {
			return x;
		}

		jz = SHR3;					// initiate, try to exit loop
		iz = (jz & 255);
		if (jz < ke[iz]) return jz*we[iz];
	}
}

// -----------------------------------------------------------------------------
//  this procedure sets the seed and creates the tables
// -----------------------------------------------------------------------------
void zigset(						// set the seed and create the tables
	unsigned long jsrseed)				// new seed
{
	const double m1 = 2147483648.0;
	const double m2 = 4294967296.0;

	double dn = 3.442619855899;
	double tn = dn;
	double vn = 9.91256303526217e-3;
	double q;

	double de = 7.697117470131487;
	double te = de;
	double ve = 3.949659822581572e-3;
	int i;

	jsr^=jsrseed;

	// -------------------------------------------------------------------------
	//  set up tables for RNOR
	// -------------------------------------------------------------------------
	q = vn / exp(-0.5 * dn * dn);
	kn[0] = (dn / q) * m1;
	kn[1] = 0;

	wn[0] = q / m1;
	wn[127] = dn / m1;

	fn[0] = 1.0;
	fn[127] = exp(-0.5 * dn * dn);

	for (i = 126; i >= 1; i--) {
		dn = sqrt(-2.0 * log(vn / dn + exp(-0.5 * dn * dn)));
		kn[i+1] = (dn / tn) * m1;
		tn = dn;
		fn[i] = exp(-0.5 * dn * dn);
		wn[i] = dn / m1;
	}

	// -------------------------------------------------------------------------
	//  set up tables for REXP
	// -------------------------------------------------------------------------
	q = ve / exp(-de);
	ke[0] = (de / q) * m2;
	ke[1] = 0;

	we[0] = q / m2;
	we[255] = de / m2;

	fe[0] = 1.0;
	fe[255] = exp(-de);

	for (i = 254; i >= 1; i--) {
		de = -log(ve / de + exp(-de));
		ke[i+1] = (de / te) * m2;
		te = de;
		fe[i] = exp(-de);
		we[i] = de / m2;
	}
}

// -----------------------------------------------------------------------------
static int PDISTComp(               // compare func for quick (ascending)   
	const void *xv,                     // 1st element
	const void *yv)                     // 2nd element
{
	float x,y;

	x = ((PDIST_PAIR *)xv)->u.pdist;
	y = ((PDIST_PAIR *)yv)->u.pdist;

	if (x < y) return -1;
	if (x > y) return 1;
	return 0;
}

// -----------------------------------------------------------------------------
static int PDISTCompDesc(           // compare func for quick (descending)
	const void *xv,                     // 1st element
	const void *yv)                     // 2nd element
{
	float x,y;

	x = ((PDIST_PAIR *)xv)->u.pdist;
	y = ((PDIST_PAIR *)yv)->u.pdist;

	if (x < y) return 1;
	if (x > y) return -1;
	return 0;
}

// -----------------------------------------------------------------------------
//  Ascending order for <rank>, if tie, descending order for <times_achieved>
// -----------------------------------------------------------------------------
static int RTAComp(                 // compare func for quick
	const void *xv,                     // 1st element
	const void *yv)                     // 2nd element
{
	int x, y;

	x = ((PDIST_PAIR *)xv)->u.rta.rank;
	y = ((PDIST_PAIR *)yv)->u.rta.rank;

	if (x < y) return -1;
	if (x > y) return 1;

	x = ((PDIST_PAIR *)xv)->u.rta.times_achieved;
	y = ((PDIST_PAIR *)yv)->u.rta.times_achieved;

	if (x < y) return 1;
	if (x > y) return -1;

	return 0;
}

// -----------------------------------------------------------------------------
QDAFN::QDAFN(						// constructor
	int   n,							// number of data objects
	int   d,							// number of dimensions
	int   L,							// number of projections
	int   M,							// number of candidates
    int   algo,							// which algorithm
	float ratio,						// approximation ratio
	const float **data)			       	// data objects
{
	// -------------------------------------------------------------------------
	//  init the input parameters
	// -------------------------------------------------------------------------
	n_pts_      = n;
	dim_        = d;
	L_          = L;
	M_          = M;
	algo_       = algo;
	appr_ratio_ = ratio;
	data_       = data;

	// -------------------------------------------------------------------------
	//  calc parameters 
	// -------------------------------------------------------------------------
	if (L_ == 0 || M_ == 0) {
		L_ = 2 * (int) ceil(pow((float) n_pts_, 
			1.0F/(appr_ratio_*appr_ratio_)));	
		if (L_ < 1) {
			printf("bad number of projections: L = %d\n", L_);
			exit(1);
		}

		float x = pow(log((float) n_pts_), (appr_ratio_*appr_ratio_/2.0F 
			- 1.0F/3.0F));
		M_ = 1 + (int) ceil(E * E * L_ * x);
		if (M_ < 1) {
			printf("bad number of candidates: M = %d\n", M_);
			exit(1);
		}
	}
	
	// -------------------------------------------------------------------------
	//  generate hash functions
	// -------------------------------------------------------------------------
	// zigset(MAGIC + time(NULL));
	zigset(MAGIC + 17);

	int size = L_ * dim_;
	proj_ = new float[size]; 
	for (int i = 0; i < size; ++i) {
		proj_[i] = RNOR / sqrt((float) dim_);
	}

	// -------------------------------------------------------------------------
	//  bulkloading
	// -------------------------------------------------------------------------
	bulkload();
}

// -----------------------------------------------------------------------------
QDAFN::~QDAFN()						// destrcutor
{
	delete[] pdp_; pdp_ = NULL; 
	delete[] proj_; proj_ = NULL;
}

// -----------------------------------------------------------------------------
int QDAFN::bulkload()				// build index
{
	// -------------------------------------------------------------------------
	//  project all the points
	// -------------------------------------------------------------------------
	pdp_ = new PDIST_PAIR[(L_ + 1) * n_pts_];

	for (int i = 0; i < L_; ++i) {
		for (int j = 0;j < n_pts_; ++j) {
			float x = 0.0f;
			for (int k = 0; k < dim_; ++k) {
				x += proj_[i * dim_ + k] * data_[j][k];
			}
			pdp_[(i+1) * n_pts_ + j].obj     = j + 1;
			pdp_[(i+1) * n_pts_ + j].u.pdist = x;
		}
	}

	// -------------------------------------------------------------------------
	//  compute master ranks
	// -------------------------------------------------------------------------
	if (algo_ == 1) {
		// ---------------------------------------------------------------------
		//  based on ranks and times achieved
		// ---------------------------------------------------------------------

		// ---------------------------------------------------------------------
		//  1. sort within each projection
		// ---------------------------------------------------------------------
		for (int i = 1; i <= L_; ++i) {
			qsort(pdp_+i*n_pts_, n_pts_, sizeof(PDIST_PAIR), PDISTComp);
		}

		// ---------------------------------------------------------------------
		//  2. assign rank numbers
		// ---------------------------------------------------------------------
		for (int i = 1; i <= L_; ++i) {
			for (int j = 0; j < n_pts_; ++j) {
				pdp_[i*n_pts_+j].u.rta.rank = j;
			}
		}

		// ---------------------------------------------------------------------
		//  3. find most extreme rank number for each point
		// ---------------------------------------------------------------------
		for (int j = 0; j < n_pts_; ++j) {
			pdp_[j].obj = j + 1;
			pdp_[j].u.rta.rank = n_pts_;
			pdp_[j].u.rta.times_achieved = 0;
		}

		for (int i = 1; i <= L_; ++i) {
			for (int j = 0; j < n_pts_; ++j) {
				if (pdp_[i*n_pts_+j].u.rta.rank < pdp_[pdp_[i*n_pts_+j].obj-1].u.rta.rank) {
					pdp_[pdp_[i*n_pts_+j].obj-1].u.rta.rank = pdp_[i*n_pts_+j].u.rta.rank;
					pdp_[pdp_[i*n_pts_+j].obj-1].u.rta.times_achieved = 1;
				}
				else if (pdp_[i*n_pts_+j].u.rta.rank == pdp_[pdp_[i*n_pts_+j].obj-1].u.rta.rank) {
					pdp_[pdp_[i*n_pts_+j].obj-1].u.rta.times_achieved++;
				}
			}
		}

		// ---------------------------------------------------------------------
		//  4. sort on those
		// ---------------------------------------------------------------------
		qsort(pdp_, n_pts_, sizeof(PDIST_PAIR), RTAComp);
	}
	else {
		// ---------------------------------------------------------------------
		//  based on projected value
		// ---------------------------------------------------------------------

		// ---------------------------------------------------------------------
		//  1. sort within each projection (extra)
		// ---------------------------------------------------------------------
		for (int i = 1; i <= L_; ++i) {
			qsort(pdp_ + i * n_pts_, n_pts_, sizeof(PDIST_PAIR), PDISTComp);
		}

		// ---------------------------------------------------------------------
		//  2. find most extreme projected value for each point
		// ---------------------------------------------------------------------
		for (int j = 0; j < n_pts_; ++j) {
			pdp_[j].obj = j + 1;
			pdp_[j].u.pdist = 1.0e38;
		}
		for (int i = 1; i <= L_; ++i) {
			for (int j = 0; j < n_pts_; ++j) {
				if (pdp_[i*n_pts_+j].u.pdist < pdp_[pdp_[i*n_pts_+j].obj-1].u.pdist) {
					pdp_[pdp_[i*n_pts_+j].obj-1].u.pdist = pdp_[i*n_pts_+j].u.pdist;
				}
			}
		}

		// ---------------------------------------------------------------------
		//  3. sort on those
		// ---------------------------------------------------------------------
		qsort(pdp_, n_pts_, sizeof(PDIST_PAIR), PDISTComp);
	}
	return 0;
}

// -----------------------------------------------------------------------------
void QDAFN::display() 				// display parameters
{
	printf("Parameters of QDAFN (SISAP2015 paper):\n");
	printf("    n    = %d\n",   n_pts_);
	printf("    d    = %d\n",   dim_);
	printf("    L    = %d\n",   L_);
	printf("    M    = %d\n",   M_);
	printf("    c    = %.1f\n", appr_ratio_);
	printf("    algo = %s\n\n", algoname[algo_]);
}

// -----------------------------------------------------------------------------
int QDAFN::kfn(						// c-k-AFN search
	int   top_k,						// top-k value
	const float *query,					// input query
	MaxK_List *list)					// c-k-AFN results (return)
{
	int candidates = M_ + top_k;
	if (candidates > n_pts_) candidates = n_pts_;

	if (algo_ == 2) {
		// ---------------------------------------------------------------------
		//  query dependent search by projected value
		// ---------------------------------------------------------------------
		int   *next    = new int[L_]; 
		float *proj_q  = new float[L_];
		bool  *checked = new bool[n_pts_];

		for (int i = 0; i < n_pts_; ++i) checked[i] = false;

		for (int i = 0; i < L_; ++i) {
			next[i] = 0;
			float x = 0.0f;
			for (int j = 0; j < dim_; ++j) {
				x += proj_[i*dim_+j] * query[j];
			}
			proj_q[i] = x;
		}

		for (int i = 0; i < candidates; ++i) {
			int   found_next = pdp_[n_pts_ + next[0]].obj;
			float y = fabs(pdp_[n_pts_+next[0]].u.pdist - proj_q[0]);
			int   found_in_proj = 0;

			for (int j = 1; j < L_; ++j) {
				float z = fabs(pdp_[n_pts_*(j+1)+next[j]].u.pdist - proj_q[j]);
				if (z > y) {
					y = z;
					found_next = pdp_[n_pts_*(j+1)+next[j]].obj;
					found_in_proj = j;
				}
			}
			next[found_in_proj]++;

			if (!checked[found_next-1]) {
				float dist = calc_l2_dist(dim_, query, data_[found_next-1]);
				list->insert(dist, found_next);
				checked[found_next-1] = true;
			}
		}
		delete[] proj_q;  proj_q = NULL;
		delete[] next;    next = NULL;
		delete[] checked; checked = NULL;
	}
	else {
		// ---------------------------------------------------------------------
		//  query independent by rank or projected value
		// ---------------------------------------------------------------------
		for (int i = 0; i < candidates; ++i) {
			int found_next = pdp_[i].obj;
			float dist = calc_l2_dist(dim_, query, data_[found_next-1]);
			list->insert(dist, found_next);
		}
	}
	return 0;
}