#include "headers.h"

// -----------------------------------------------------------------------------
int ground_truth(					// find ground truth
	int   n,							// number of data  objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const char *data_set,				// address of data  set
	const char *query_set,				// address of query set
	const char *truth_set)				// address of truth set
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read data set and query set
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	float **data = new float*[n];
	for (int i = 0; i < n; ++i) data[i] = new float[d];
	if (read_data(n, d, data_set, data) == 1) {
		printf("Reading Dataset Error!\n");
		exit(1);
	}

	float **query = new float*[qn];
	for (int i = 0; i < qn; ++i) query[i] = new float[d];
	if (read_data(qn, d, query_set, query) == 1) {
		printf("Reading Query Set Error!\n");
		exit(1);
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Data and Query: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  find ground truth results (using linear scan method)
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	FILE *fp = fopen(truth_set, "w");
	if (!fp) {
		printf("Could not create %s.\n", truth_set);
		return 1;
	}

	MaxK_List *list = new MaxK_List(MAXK);
	fprintf(fp, "%d %d\n", qn, MAXK);
	for (int i = 0; i < qn; ++i) {
		list->reset();
		for (int j = 0; j < n; ++j) {
			float dist = calc_l2_dist(d, data[j], query[i]);
			list->insert(dist, j + 1);
		}

		for (int j = 0; j < MAXK; ++j) {
			fprintf(fp, "%d %f ", list->ith_id(j), list->ith_key(j));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	gettimeofday(&end_time, NULL);
	float truth_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Ground Truth: %f Seconds\n\n", truth_time);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete list; list = NULL;
	for (int i = 0; i < n; ++i) {
		delete[] data[i]; data[i] = NULL;
	}
	delete[] data; data = NULL;

	for (int i = 0; i < qn; ++i) {
		delete[] query[i]; query[i] = NULL;
	}
	delete[] query; query = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int rqalsh_star(					// c-k-AFN search of RQALSH*
	int   n,							// number of data  objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   L,							// number of projection (drusilla)
	int   M,							// number of candidates (drusilla)
	int   beta,							// false positive percentage
	float delta,						// error probability
	float ratio,						// approximation ratio
	const char *data_set,				// address of data  set
	const char *query_set,				// address of query set
	const char *truth_set,				// address of truth set
	const char *output_folder)			// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read data set, query set and truth set
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	float **data = new float*[n];
	for (int i = 0; i < n; ++i) data[i] = new float[d];
	if (read_data(n, d, data_set, data) == 1) {
		printf("Reading Dataset Error!\n");
		exit(1);
	}

	float **query = new float*[qn];
	for (int i = 0; i < qn; ++i) query[i] = new float[d];
	if (read_data(qn, d, query_set, query) == 1) {
		printf("Reading Query Set Error!\n");
		exit(1);
	}

	Result **R = new Result*[qn];
	for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
	if (read_ground_truth(qn, truth_set, R) == 1) {
		printf("Reading Truth Set Error!\n");
		exit(1);
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Data, Query and Ground Truth: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  Indexing
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	RQALSH_Star* lsh = new RQALSH_Star();
	lsh->build(n, d, L, M, beta, delta, ratio, (const float **) data);
	
	gettimeofday(&end_time, NULL);
	float indexing_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  c-k-Approximate Furthest Neighbor (c-k-AFN) search
	// -------------------------------------------------------------------------
	char output_set[200];
	sprintf(output_set, "%srqalsh_star_mem_L=%d_M=%d.out", output_folder, L, M);

	FILE *fp = fopen(output_set, "w");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		exit(1);
	}

	int kFNs[] = { 1, 2, 5, 10 };
	int maxRound = 4;
	int top_k = -1;

	float runtime = -1.0f;
	float overall_ratio = -1.0f;
	float recall = -1.0f;

	printf("c-k-AFN Search by RQALSH*:\n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int round = 0; round < maxRound; ++round) {
		gettimeofday(&start_time, NULL);
		top_k = kFNs[round];
		MaxK_List *list = new MaxK_List(top_k);

		overall_ratio = 0.0f;
		recall = 0.0f;

		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kfn(top_k, (const float *) query[i], list);
			recall += calc_recall(top_k, (const Result *) R[i], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += R[i][j].key_ / list->ith_key(j);
			}
			overall_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&end_time, NULL);
		runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - 
			start_time.tv_usec) / 1000000.0f;

		overall_ratio = overall_ratio / qn;
		recall        = recall / qn;
		runtime       = (runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.2f\t\t%.2f%%\n", top_k, overall_ratio, 
			runtime, recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, overall_ratio, runtime, recall);
	}
	printf("\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	for (int i = 0; i < n; ++i) {
		delete[] data[i]; data[i] = NULL;
	}
	delete[] data; data = NULL;
	
	for (int i = 0; i < qn; ++i) {
		delete[] query[i]; query[i] = NULL;
		delete[] R[i]; R[i] = NULL;
	}
	delete[] query; query = NULL;
	delete[] R; R = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int rqalsh(							// c-k-AFN search of RQALSH
	int   n,							// number of data  objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   beta,							// false positive percentage
	float delta,						// error probability
	float ratio,						// approximation ratio
	const char *data_set,				// address of data  set
	const char *query_set,				// address of query set
	const char *truth_set,				// address of truth set
	const char *output_folder)			// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read data set, query set and truth set
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	float **data = new float*[n];
	for (int i = 0; i < n; ++i) data[i] = new float[d];
	if (read_data(n, d, data_set, data) == 1) {
		printf("Reading Dataset Error!\n");
		exit(1);
	}

	float **query = new float*[qn];
	for (int i = 0; i < qn; ++i) query[i] = new float[d];
	if (read_data(qn, d, query_set, query) == 1) {
		printf("Reading Query Set Error!\n");
		exit(1);
	}

	Result **R = new Result*[qn];
	for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
	if (read_ground_truth(qn, truth_set, R) == 1) {
		printf("Reading Truth Set Error!\n");
		exit(1);
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Data, Query and Ground Truth: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  Indexing
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	RQALSH* lsh = new RQALSH();
	lsh->build(n, d, beta, delta, ratio, (const float **) data);
	
	gettimeofday(&end_time, NULL);
	float indexing_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  c-k-Approximate Furthest Neighbor (c-k-AFN) search
	// -------------------------------------------------------------------------
	char output_set[200];
	sprintf(output_set, "%srqalsh_mem.out", output_folder);

	FILE *fp = fopen(output_set, "w");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		exit(1);
	}

	int kFNs[] = { 1, 2, 5, 10 };
	int maxRound = 4;
	int top_k = -1;

	float runtime = -1.0f;
	float overall_ratio = -1.0f;
	float recall = -1.0f;

	printf("c-k-AFN Search by RQALSH:\n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int round = 0; round < maxRound; ++round) {
		gettimeofday(&start_time, NULL);
		top_k = kFNs[round];
		MaxK_List *list = new MaxK_List(top_k);

		overall_ratio = 0.0f;
		recall = 0.0f;

		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kfn(top_k, (const float *) query[i], list);
			recall += calc_recall(top_k, (const Result *) R[i], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += R[i][j].key_ / list->ith_key(j);
			}
			overall_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&end_time, NULL);
		runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - 
			start_time.tv_usec) / 1000000.0f;

		overall_ratio = overall_ratio / qn;
		recall        = recall / qn;
		runtime       = (runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.2f\t\t%.2f%%\n", top_k, overall_ratio, 
			runtime, recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, overall_ratio, runtime, recall);
	}
	printf("\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	for (int i = 0; i < n; ++i) {
		delete[] data[i]; data[i] = NULL;
	}
	delete[] data; data = NULL;
	
	for (int i = 0; i < qn; ++i) {
		delete[] query[i]; query[i] = NULL;
		delete[] R[i]; R[i] = NULL;
	}
	delete[] query; query = NULL;
	delete[] R; R = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int linear_scan(					// k-FN search of linear scan
	int   n,							// number of data  objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const char *data_set,				// address of data  set
	const char *query_set,				// address of query set
	const char *truth_set,				// address of truth set
	const char *output_folder)			// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read data set, query set and truth set
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	float **data = new float*[n];
	for (int i = 0; i < n; ++i) data[i] = new float[d];
	if (read_data(n, d, data_set, data) == 1) {
		printf("Reading Dataset Error!\n");
		exit(1);
	}

	float **query = new float*[qn];
	for (int i = 0; i < qn; ++i) query[i] = new float[d];
	if (read_data(qn, d, query_set, query) == 1) {
		printf("Reading Query Set Error!\n");
		exit(1);
	}

	Result **R = new Result*[qn];
	for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
	if (read_ground_truth(qn, truth_set, R) == 1) {
		printf("Reading Truth Set Error!\n");
		exit(1);
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Data, Query and Ground Truth: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  k-FN search
	// -------------------------------------------------------------------------
	char output_set[200];
	sprintf(output_set, "%slinear_mem.out", output_folder);

	FILE *fp = fopen(output_set, "w");
	if (!fp) {
		printf("Could not create %s.\n", output_set);
		return 1;
	}

	int kFNs[] = { 1, 2, 5, 10 };
	int maxRound = 4;
	int top_k = -1;

	float runtime = -1.0f;
	float overall_ratio = -1.0f;
	float recall = -1.0f;

	printf("k-FN Search by Linear Scan:\n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int round = 0; round < maxRound; ++round) {
		gettimeofday(&start_time, NULL);
		top_k = kFNs[round];
		MaxK_List *list = new MaxK_List(top_k);

		overall_ratio = 0.0f;
		recall = 0.0f;
		
		for (int i = 0; i < qn; ++i) {
			list->reset();
			for (int j = 0; j < n; ++j) {
				float dist = calc_l2_dist(d, data[j], query[i]);
				list->insert(dist, j);
			}
			recall += calc_recall(top_k, (const Result *) R[i], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += R[i][j].key_ / list->ith_key(j);
			}
			overall_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&end_time, NULL);
		runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - 
			start_time.tv_usec) / 1000000.0f;

		overall_ratio = overall_ratio / qn;
		recall        = recall / qn;
		runtime       = (runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.2f\t\t%.2f%%\n", top_k, overall_ratio, 
			runtime, recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, overall_ratio, runtime, recall);
	}
	printf("\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	for (int i = 0; i < n; ++i) {
		delete[] data[i]; data[i] = NULL;
	}
	delete[] data; data = NULL;
	
	for (int i = 0; i < qn; ++i) {
		delete[] query[i]; query[i] = NULL;
		delete[] R[i]; R[i] = NULL;
	}
	delete[] query; query = NULL;
	delete[] R; R = NULL;
	
	return 0;
}