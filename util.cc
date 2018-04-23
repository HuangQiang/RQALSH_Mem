#include "headers.h"

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
int create_dir(						// create directory
	char *path)							// input path
{
	// ---------------------------------------------------------------------
	//  check whether the directory exists. if the directory does not 
	//  exist, we create the directory for each folder)
	// ---------------------------------------------------------------------
	int len = (int)strlen(path);
	for (int i = 0; i < len; ++i) {
		if (path[i] == '/') {
			char ch = path[i + 1];
			path[i + 1] = '\0';

			int ret = access(path, F_OK);
			if (ret != 0) {
				ret = mkdir(path, 0755);
				if (ret != 0) {
					printf("Could not create directory %s\n", path);
					return 1;
				}
			}
			path[i + 1] = ch;
		}
	}
	return 0;
}

// -----------------------------------------------------------------------------
float calc_l2_dist(					// calc L_2 norm (data type is float)
	int dim,							// dimension
	const float* vec1,					// 1st point
	const float* vec2)					// 2nd point
{
	float diff = 0.0F;
	float ret  = 0.0F;

	for (int i = 0; i < dim; i++) {
		diff = vec1[i] - vec2[i];
		ret += diff * diff;
	}
	return sqrt(ret);
}

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall (percentage)
	int  k,								// top-k value
	const Result *R,					// ground truth results 
	MaxK_List *list)					// results returned by algorithms
{
	int i = k - 1;
	int last = k - 1;
	while (i >= 0 && list->ith_key(i) < R[last].key_) {
		i--;
	}
	return (i + 1) * 100.0f / k;
}

// -----------------------------------------------------------------------------
int read_data(						// read data/query set from disk
	int   n,							// number of data/query objects
	int   d,							// dimensionality
	const char *fname,					// address of data/query set
	float **data)						// data/query objects (return)
{
	FILE *fp = fopen(fname, "r");
	if (!fp) {
		printf("Could not open %s\n", fname);
		return 1;
	}

	int i   = 0;
	int tmp = -1;
	while (!feof(fp) && i < n) {
		fscanf(fp, "%d", &tmp);
		for (int j = 0; j < d; ++j) {
			fscanf(fp, " %f", &data[i][j]);
		}
		fscanf(fp, "\n");
		++i;
	}
	assert(feof(fp) && i == n);
	fclose(fp);

	return 0;
}

// -----------------------------------------------------------------------------
int read_ground_truth(				// read ground truth results from disk
	int qn,								// number of query objects
	const char *fname,					// address of truth set
	Result **R)							// ground truth results (return)
{
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

	return 0;
}