#include "headers.h"

timeval g_start_time;
timeval g_end_time;

float g_runtime = -1.0f;
float g_ratio   = -1.0f;
float g_recall  = -1.0f;

// -----------------------------------------------------------------------------
int ResultComp(						// compare function for qsort (ascending)
	const void *e1,						// 1st element
	const void *e2)						// 2nd element
{
	int ret = 0;
	Result *item1 = (Result*) e1;
	Result *item2 = (Result*) e2;

	if (item1->key_ < item2->key_) {
		ret = -1;
	} 
	else if (item1->key_ > item2->key_) {
		ret = 1;
	} 
	else {
		if (item1->id_ < item2->id_) ret = -1;
		else if (item1->id_ > item2->id_) ret = 1;
	}
	return ret;
}

// -----------------------------------------------------------------------------
int ResultCompDesc(					// compare function for qsort (descending)
	const void *e1,						// 1st element
	const void *e2)						// 2nd element
{
	int ret = 0;
	Result *item1 = (Result*) e1;
	Result *item2 = (Result*) e2;

	if (item1->key_ < item2->key_) {
		ret = 1;
	} 
	else if (item1->key_ > item2->key_) {
		ret = -1;
	} 
	else {
		if (item1->id_ < item2->id_) ret = -1;
		else if (item1->id_ > item2->id_) ret = 1;
	}
	return ret;
}

// -------------------------------------------------------------------------
void create_dir(					// create directory
	char *path)							// input path
{
	int len = (int)strlen(path);
	for (int i = 0; i < len; ++i) {
		if (path[i] == '/') {
			char ch = path[i + 1];
			path[i + 1] = '\0';
									// check whether the directory exists
			int ret = access(path, F_OK);
			if (ret != 0) {			// create the directory
				ret = mkdir(path, 0755);
				if (ret != 0) printf("Could not create directory %s\n", path);
			}
			path[i + 1] = ch;
		}
	}
}

// -----------------------------------------------------------------------------
int read_data(						// read data/query set from disk
	int   n,							// number of data/query objects
	int   d,							// dimensionality
	const char *fname,					// address of data/query set
	float **data)						// data/query objects (return)
{
	gettimeofday(&g_start_time, NULL);
	FILE *fp = fopen(fname, "r");
	if (!fp) {
		printf("Could not open %s.\n", fname);
		return 1;
	}

	int i = 0;
	int j = 0;
	while (!feof(fp) && i < n) {
		fscanf(fp, "%d", &j);
		for (j = 0; j < d; ++j) {
			fscanf(fp, " %f", &data[i][j]);
		}
		fscanf(fp, "\n");
		++i;
	}
	assert(feof(fp) && i == n);
	fclose(fp);

	gettimeofday(&g_end_time, NULL);
	float running_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Read Data: %f Seconds\n\n", running_time);

	return 0;
}

// -----------------------------------------------------------------------------
int read_ground_truth(				// read ground truth results from disk
	int qn,								// number of query objects
	const char *fname,					// address of truth set
	Result **R)							// ground truth results (return)
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
			fscanf(fp, "%d %f ", &R[i][j].id_, &R[i][j].key_);
		}
		fscanf(fp, "\n");
	}
	fclose(fp);

	gettimeofday(&g_end_time, NULL);
	float running_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Read Ground Truth: %f Seconds\n\n", running_time);

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
	const float **data,					// data set
	const float **query,				// query set
	const char *truth_set)				// address of truth set
{
	gettimeofday(&g_start_time, NULL);
	FILE *fp = fopen(truth_set, "w");
	if (!fp) {
		printf("Could not create %s\n", truth_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  find ground truth results (using linear scan method)
	// -------------------------------------------------------------------------
	Result **result = new Result*[qn];
	for (int i = 0; i < qn; ++i) {
		result[i] = new Result[MAXK];
	}
	k_fn_search(n, qn, d, MAXK, data, query, result);

	fprintf(fp, "%d %d\n", qn, MAXK);
	for (int i = 0; i < qn; ++i) {
		for (int j = 0; j < MAXK; ++j) {
			fprintf(fp, "%d %f ", result[i][j].id_, result[i][j].key_);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	for (int i = 0; i < qn; ++i) {
		delete[] result[i]; result[i] = NULL;
	}
	delete[] result; result = NULL;

	gettimeofday(&g_end_time, NULL);
	float truth_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Ground Truth: %f Seconds\n\n", truth_time);

	return 0;
}

// -----------------------------------------------------------------------------
void k_fn_search(					// k-NN search
	int   n, 							// cardinality
	int   qn,							// query number
	int   d, 							// dimensionality
	int   k,							// top-k value
	const float **data,					// data objects
	const float **query,				// query objects
	Result **result)					// k-MIP results (return)
{
	// -------------------------------------------------------------------------
	//  k-FN search by linear scan
	// -------------------------------------------------------------------------
	MaxK_List *list = new MaxK_List(k);
	for (int i = 0; i < qn; ++i) {
		list->reset();
		for (int j = 0; j < n; ++j) {
			float dist = calc_l2_dist(d, data[j], query[i]);
			list->insert(dist, j + 1);
		}

		for (int j = 0; j < k; ++j) {
			result[i][j].id_  = list->ith_id(j);
			result[i][j].key_ = list->ith_key(j);
		}
	}
	delete list; list = NULL;
}