#include "util.h"

timeval g_start_time;				// global param: start time
timeval g_end_time;					// global param: end time

float g_memory    = -1.0f;			// global param: estimated memory usage (MB)
float g_indextime = -1.0f;			// global param: indexing time (seconds)

float g_runtime   = -1.0f;			// global param: running time (ms)
float g_ratio     = -1.0f;			// global param: overall ratio
float g_recall    = -1.0f;			// global param: recall (%)
float g_fraction  = -1.0f;			// global param: fraction (%)

// -----------------------------------------------------------------------------
void create_dir(					// create directory
	char *path)							// input path
{
	int len = (int) strlen(path);
	for (int i = 0; i < len; ++i) {
		if (path[i] == '/') {
			char ch = path[i + 1];
			path[i + 1] = '\0';

			int ret = access(path, F_OK);
			if (ret != 0) {			// create directory if not exists
				ret = mkdir(path, 0755);
				if (ret != 0) {
					printf("Could not create directory %s\n", path);
					exit(1);
				}
			}
			path[i + 1] = ch;
		}
	}
}

// -----------------------------------------------------------------------------
int read_bin_data(					// read data (binary) from disk
	int   n,							// number of data points
	int   d,							// dimensionality
	bool  is_data,						// is dataset?
	const char *fname,					// address of data
	float *data)						// data (return)
{
	gettimeofday(&g_start_time, NULL);
	FILE *fp = fopen(fname, "rb");
	if (!fp) { printf("Could not open %s\n", fname); return 1; }

	fread(data, SIZEFLOAT, n * d, fp);
	fclose(fp);

	gettimeofday(&g_end_time, NULL);
	float running_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	
	if (is_data) printf("Read Data:  %f Seconds\n", running_time);
	else printf("Read Query: %f Seconds\n", running_time);
	
	return 0;
}

// -----------------------------------------------------------------------------
int read_ground_truth(				// read ground truth results from disk
	int qn,								// number of query objects
	const char *fname,					// address of truth set
	Result *R)							// ground truth results (return)
{
	gettimeofday(&g_start_time, NULL);
	FILE *fp = fopen(fname, "r");
	if (!fp) {
		printf("Could not open %s\n", fname);
		return 1;
	}

	int tmp1 = -1;
	int tmp2 = -1;
	fscanf(fp, "%d %d\n", &tmp1, &tmp2);
	assert(tmp1 == qn && tmp2 == MAXK);

	for (int i = 0; i < qn; ++i) {
		for (int j = 0; j < MAXK; ++j) {
			fscanf(fp, "%d %f ", &R[i*MAXK+j].id_, &R[i*MAXK+j].key_);
		}
		fscanf(fp, "\n");
	}
	fclose(fp);

	gettimeofday(&g_end_time, NULL);
	float running_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Read Truth: %f Seconds\n\n", running_time);

	return 0;
}

// -----------------------------------------------------------------------------
float calc_l2_dist(					// calc L_2 norm (data type is float)
	int   dim,							// dimension
	const float *p1,					// 1st point
	const float *p2)					// 2nd point
{
	float ret  = 0.0F;
	for (int i = 0; i < dim; ++i) {
		ret += SQR(p1[i] - p2[i]);
	}
	return sqrt(ret);
}

// -----------------------------------------------------------------------------
float calc_inner_product(			// calc inner product (data type is float)
	int   dim,							// dimension
	const float *p1,					// 1st point
	const float *p2)					// 2nd point
{
	float ret = 0.0F;
	for (int i = 0; i < dim; ++i) {
		ret += p1[i] * p2[i];
	}
	return ret;
}

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall (percentage)
	int   k,							// top-k value
	const Result *R,					// ground truth results 
	MaxK_List *list)					// results returned by algorithms
{
	int i = k - 1;
	int last = k - 1;
	while (i >= 0 && list->ith_key(i) < R[last].key_) {
		--i;
	}
	return (i + 1) * 100.0f / k;
}

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall (percentage)
	int   k,							// top-k value
	const Result *R,					// ground truth results 
	const Result *result)				// results returned by algorithms
{
	int i = k - 1;
	int last = k - 1;
	while (i >= 0 && result[i].key_ < R[last].key_) {
		--i;
	}
	return (i + 1) * 100.0f / k;
}

// -----------------------------------------------------------------------------
int ground_truth(					// find ground truth
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const float *data,					// data set
	const float *query,					// query set
	const char  *truth_set)				// address of truth set
{
	gettimeofday(&g_start_time, NULL);
	FILE *fp = fopen(truth_set, "w");
	if (!fp) { printf("Could not create %s\n", truth_set); return 1; }

	// -------------------------------------------------------------------------
	//  find ground truth results (using linear scan method)
	// -------------------------------------------------------------------------
	fprintf(fp, "%d %d\n", qn, MAXK);

	MaxK_List *list = new MaxK_List(MAXK);
	for (int i = 0; i < qn; ++i) {
		list->reset();
		k_fn_search(n, d, data, &query[i*d], list);
		for (int j = 0; j < MAXK; ++j) {
			fprintf(fp, "%d %f ", list->ith_id(j), list->ith_key(j));
		}
		fprintf(fp, "\n");
	}
	delete list;
	fclose(fp);

	gettimeofday(&g_end_time, NULL);
	float truth_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Ground Truth: %f Seconds\n\n", truth_time);

	return 0;
}

// -----------------------------------------------------------------------------
int k_fn_search(					// k-FN search
	int   n, 							// cardinality
	int   d, 							// dimensionality
	const float *data,					// data objects
	const float *query,					// query objects
	MaxK_List *list)					// top-k results
{
	// -------------------------------------------------------------------------
	//  k-FN search by linear scan
	// -------------------------------------------------------------------------
	for (int j = 0; j < n; ++j) {
		float dist = calc_l2_dist(d, &data[j*d], query);
		list->insert(dist, j + 1);
	}
	return n;
}