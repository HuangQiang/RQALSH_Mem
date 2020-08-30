#include "afn.h"

// -----------------------------------------------------------------------------
int linear_scan(					// k-FN search of linear scan
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const float *data,					// data set
	const float *query,					// query set
	const Result *R, 					// truth set
	const char *out_path)				// output path
{
	char output_set[200]; sprintf(output_set, "%slinear.out", out_path);
	FILE *fp = fopen(output_set, "a+");
	if (!fp) { printf("Could not create %s\n", output_set); return 1; }

	// -------------------------------------------------------------------------
	//  c-k-AFN search
	// -------------------------------------------------------------------------
	fprintf(fp, "Linear Scan:\n");

	printf("Top-k FN Search of Linear Scan:\n");
	printf("Top-k\t\tRatio\t\tTime (ms)\tRecall (%)\tFraction (%)\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List *list = new MaxK_List(top_k);
		
		g_ratio    = 0.0f;
		g_recall   = 0.0f;
		g_fraction = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			int check_k = k_fn_search(n, d, data, &query[i*d], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				if (fabs(list->ith_key(j) - R[i*MAXK+j].key_) < CHECK_ERROR) 
					ratio += 1.0f;
				else 
					ratio += R[i*MAXK+j].key_ / list->ith_key(j);
			}
			g_ratio    += ratio / top_k;
			g_recall   += calc_recall(top_k, &R[i*MAXK], list);
			g_fraction += check_k * 100.0f / n;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

		g_ratio    = g_ratio    / qn;
		g_recall   = g_recall   / qn;
		g_fraction = g_fraction / qn;
		g_runtime  = (g_runtime * 1000.0f) / qn;

		printf("%3d\t\t%.4f\t\t%.4f\t\t%.2f%%\t\t%.2f%%\n", top_k, g_ratio, 
			g_runtime, g_recall, g_fraction);
		fprintf(fp, "%d\t%f\t%f\t%f\t%f\n", top_k, g_ratio, g_runtime, 
			g_recall, g_fraction);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);
	
	return 0;
}

// -----------------------------------------------------------------------------
int qdafn(							// c-k-AFN search of QDAFN
	int   n,							// number of data objects
	int   qn,							// number of query points
	int   d,							// number of dimensions
	int   L,							// number of projections
	int   M,							// number of candidates
	float ratio,						// approximation ratio
	const float *data,					// data set
	const float *query,					// query set
	const Result *R, 					// truth set
	const char *out_path)				// output path
{
	char output_set[200]; sprintf(output_set, "%sqdafn.out", out_path);
	FILE *fp = fopen(output_set, "a+");
	if (!fp) { printf("Could not create %s\n", output_set); return 1; }

	// -------------------------------------------------------------------------
	//  indexing 
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	QDAFN *hash = new QDAFN(n, d, L, M, 2, ratio, data);
	hash->display();

	gettimeofday(&g_end_time, NULL);
	g_indextime = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	g_memory = hash->get_memory_usage() / 1048576.0f;
	
	printf("Indexing Time = %f Seconds\n", g_indextime);
	printf("Memory = %f MB\n\n", g_memory);
	
	fprintf(fp, "QDAFN: L=%d, M=%d\n", L, M);
	fprintf(fp, "Indexing Time: %f Seconds\n", g_indextime);
	fprintf(fp, "Estimated Memory: %f MB\n", g_memory);

	// -------------------------------------------------------------------------
	//  c-k-AFN search of QDAFN
	// -------------------------------------------------------------------------
	printf("Top-k FN Search of QDAFN: \n");
	printf("Top-k\t\tRatio\t\tTime (ms)\tRecall (%)\tFraction (%)\n");	
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List *list = new MaxK_List(top_k);

		g_ratio    = 0.0f;
		g_recall   = 0.0f;
		g_fraction = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			int check_k = hash->kfn(top_k, &query[i*d], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				if (fabs(list->ith_key(j) - R[i*MAXK+j].key_) < CHECK_ERROR) 
					ratio += 1.0f;
				else 
					ratio += R[i*MAXK+j].key_ / list->ith_key(j);
			}
			g_ratio    += ratio / top_k;
			g_recall   += calc_recall(top_k, &R[i*MAXK], list);
			g_fraction += check_k * 100.0f / n;

		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

		g_ratio    = g_ratio    / qn;
		g_recall   = g_recall   / qn;
		g_fraction = g_fraction / qn;
		g_runtime  = (g_runtime * 1000.0f) / qn;

		printf("%3d\t\t%.4f\t\t%.4f\t\t%.2f%%\t\t%.2f%%\n", top_k, g_ratio, 
			g_runtime, g_recall, g_fraction);
		fprintf(fp, "%d\t%f\t%f\t%f\t%f\n", top_k, g_ratio, g_runtime, 
			g_recall, g_fraction);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);
	delete hash;

	return 0;
}

// -----------------------------------------------------------------------------
int drusilla_select(				// c-k-AFN search of Drusilla-Select
	int   n,							// number of data objects
	int   qn,							// number of query points
	int   d,							// number of dimensions
	int   L,							// number of projections
	int   M,							// number of candidates
	const float *data,					// data set
	const float *query,					// query set
	const Result *R, 					// truth set
	const char *out_path)				// output path
{
	char output_set[200]; sprintf(output_set, "%sdrusilla_select.out", out_path);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) { printf("Could not create %s\n", output_set); return 1; }

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	Drusilla_Select *drusilla = new Drusilla_Select(n, d, L, M, data);
	drusilla->display();

	gettimeofday(&g_end_time, NULL);
	g_indextime = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	g_memory = drusilla->get_memory_usage() / 1048576.0f;

	printf("Indexing Time = %f Seconds\n", g_indextime);
	printf("Memory = %f MB\n\n", g_memory);
	
	fprintf(fp, "Drusilla Select: L=%d, M=%d\n", L, M);
	fprintf(fp, "Indexing Time: %f Seconds\n", g_indextime);
	fprintf(fp, "Estimated Memory: %f MB\n", g_memory);

	// -------------------------------------------------------------------------
	//  c-k-AFN search of Drusilla-Select
	// -------------------------------------------------------------------------
	printf("Top-k FN Search of Drusilla-Select: \n");
	printf("Top-k\t\tRatio\t\tTime (ms)\tRecall (%)\tFraction (%)\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List *list = new MaxK_List(top_k);

		g_ratio    = 0.0f;
		g_recall   = 0.0f;
		g_fraction = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			int check_k = drusilla->kfn(&query[i*d], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				if (fabs(list->ith_key(j) - R[i*MAXK+j].key_) < CHECK_ERROR) 
					ratio += 1.0f;
				else 
					ratio += R[i*MAXK+j].key_ / list->ith_key(j);
			}
			g_ratio    += ratio / top_k;
			g_recall   += calc_recall(top_k, &R[i*MAXK], list);
			g_fraction += check_k * 100.0f / n;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

		g_ratio    = g_ratio    / qn;
		g_recall   = g_recall   / qn;
		g_fraction = g_fraction / qn;
		g_runtime  = (g_runtime * 1000.0f) / qn;

		printf("%3d\t\t%.4f\t\t%.4f\t\t%.2f%%\t\t%.2f%%\n", top_k, g_ratio, 
			g_runtime, g_recall, g_fraction);
		fprintf(fp, "%d\t%f\t%f\t%f\t%f\n", top_k, g_ratio, g_runtime, 
			g_recall, g_fraction);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);
	delete drusilla;
	
	return 0;
}

// -----------------------------------------------------------------------------
int rqalsh(							// c-k-AFN search of RQALSH
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	float ratio,						// approximation ratio
	const float *data,					// data set
	const float *query,					// query set
	const Result *R, 					// truth set
	const char *out_path)				// output path
{
	char output_set[200]; sprintf(output_set, "%srqalsh.out", out_path);
	FILE *fp = fopen(output_set, "a+");
	if (!fp) { printf("Could not create %s\n", output_set); return 1; }

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	RQALSH* lsh = new RQALSH(n, d, ratio, NULL, data);
	lsh->display();
	
	gettimeofday(&g_end_time, NULL);
	g_indextime = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	g_memory = lsh->get_memory_usage() / 1048576.0f;

	printf("Indexing Time = %f Seconds\n", g_indextime);
	printf("Memory = %f MB\n\n", g_memory);
	
	fprintf(fp, "RQALSH:\n");
	fprintf(fp, "Indexing Time: %f Seconds\n", g_indextime);
	fprintf(fp, "Estimated Memory: %f MB\n", g_memory);

	// -------------------------------------------------------------------------
	//  c-k-AFN search
	// -------------------------------------------------------------------------
	printf("Top-k FN Search of RQALSH:\n");
	printf("Top-k\t\tRatio\t\tTime (ms)\tRecall (%)\tFraction (%)\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List *list = new MaxK_List(top_k);

		g_ratio    = 0.0f;
		g_recall   = 0.0f;
		g_fraction = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			int check_k = lsh->kfn(top_k, MINREAL, &query[i*d], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				if (fabs(list->ith_key(j) - R[i*MAXK+j].key_) < CHECK_ERROR) 
					ratio += 1.0f;
				else 
					ratio += R[i*MAXK+j].key_ / list->ith_key(j);
			}
			g_ratio    += ratio / top_k;
			g_recall   += calc_recall(top_k, &R[i*MAXK], list);
			g_fraction += check_k * 100.0f / n;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

		g_ratio    = g_ratio    / qn;
		g_recall   = g_recall   / qn;
		g_fraction = g_fraction / qn;
		g_runtime  = (g_runtime * 1000.0f) / qn;

		printf("%3d\t\t%.4f\t\t%.4f\t\t%.2f%%\t\t%.2f%%\n", top_k, g_ratio, 
			g_runtime, g_recall, g_fraction);
		fprintf(fp, "%d\t%f\t%f\t%f\t%f\n", top_k, g_ratio, g_runtime, 
			g_recall, g_fraction);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);
	delete lsh;

	return 0;
}

// -----------------------------------------------------------------------------
int rqalsh_star(					// c-k-AFN search of RQALSH*
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   L,							// number of projection (drusilla)
	int   M,							// number of candidates (drusilla)
	float ratio,						// approximation ratio
	const float *data,					// data set
	const float *query,					// query set
	const Result *R, 					// truth set
	const char *out_path)				// output path
{
	char output_set[200]; sprintf(output_set, "%srqalsh_star.out", out_path);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) { printf("Could not create %s\n", output_set); return 1; }

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	RQALSH_STAR* lsh = new RQALSH_STAR(n, d, L, M, ratio, data);
	lsh->display();
	
	gettimeofday(&g_end_time, NULL);
	g_indextime = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	g_memory = lsh->get_memory_usage() / 1048576.0f;

	printf("Indexing Time = %f Seconds\n", g_indextime);
	printf("Memory = %f MB\n\n", g_memory);
	
	fprintf(fp, "RQALSH*: L=%d, M=%d\n", L, M);
	fprintf(fp, "Indexing Time: %f Seconds\n", g_indextime);
	fprintf(fp, "Estimated Memory: %f MB\n", g_memory);

	// -------------------------------------------------------------------------
	//  c-k-AFN search
	// -------------------------------------------------------------------------
	printf("Top-k FN Search of RQALSH*:\n");
	printf("Top-k\t\tRatio\t\tTime (ms)\tRecall (%)\tFraction (%)\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List *list = new MaxK_List(top_k);

		g_ratio    = 0.0f;
		g_recall   = 0.0f;
		g_fraction = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			int check_k = lsh->kfn(top_k, &query[i*d], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				if (fabs(list->ith_key(j) - R[i*MAXK+j].key_) < CHECK_ERROR) 
					ratio += 1.0f;
				else 
					ratio += R[i*MAXK+j].key_ / list->ith_key(j);
			}
			g_ratio    += ratio / top_k;
			g_recall   += calc_recall(top_k, &R[i*MAXK], list);
			g_fraction += check_k * 100.0f / n;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

		g_ratio    = g_ratio    / qn;
		g_recall   = g_recall   / qn;
		g_fraction = g_fraction / qn;
		g_runtime  = (g_runtime * 1000.0f) / qn;

		printf("%3d\t\t%.4f\t\t%.4f\t\t%.2f%%\t\t%.2f%%\n", top_k, g_ratio, 
			g_runtime, g_recall, g_fraction);
		fprintf(fp, "%d\t%f\t%f\t%f\t%f\n", top_k, g_ratio, g_runtime, 
			g_recall, g_fraction);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);
	delete lsh;

	return 0;
}

// -----------------------------------------------------------------------------
int ml_rqalsh(						// c-k-AFN search of ML-RQALSH
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	float ratio,						// approximation ratio
	const float *data,					// data set
	const float *query,					// query set
	const Result *R, 					// truth set
	const char *out_path)				// output path
{
	char output_set[200]; sprintf(output_set, "%sml_rqalsh.out", out_path);
	FILE *fp = fopen(output_set, "a+");
	if (!fp) { printf("Could not create %s\n", output_set); return 1; }

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	ML_RQALSH* lsh = new ML_RQALSH(n, d, ratio, data);
	lsh->display();

	gettimeofday(&g_end_time, NULL);
	g_indextime = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	g_memory = lsh->get_memory_usage() / 1048576.0f;

	printf("Indexing Time = %f Seconds\n", g_indextime);
	printf("Memory = %f MB\n\n", g_memory);
	
	fprintf(fp, "ML_RQALSH:\n");
	fprintf(fp, "Indexing Time: %f Seconds\n", g_indextime);
	fprintf(fp, "Estimated Memory: %f MB\n", g_memory);

	// -------------------------------------------------------------------------
	//  c-k-AFN search
	// -------------------------------------------------------------------------
	printf("Top-k FN Search of ML_RQALSH:\n");
	printf("Top-k\t\tRatio\t\tTime (ms)\tRecall (%)\tFraction (%)\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List *list = new MaxK_List(top_k);

		g_ratio    = 0.0f;
		g_recall   = 0.0f;
		g_fraction = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			int check_k = lsh->kfn(top_k, &query[i*d], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				if (fabs(list->ith_key(j) - R[i*MAXK+j].key_) < CHECK_ERROR) 
					ratio += 1.0f;
				else 
					ratio += R[i*MAXK+j].key_ / list->ith_key(j);
			}
			g_ratio    += ratio / top_k;
			g_recall   += calc_recall(top_k, &R[i*MAXK], list);
			g_fraction += check_k * 100.0f / n;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

		g_ratio    = g_ratio    / qn;
		g_recall   = g_recall   / qn;
		g_fraction = g_fraction / qn;
		g_runtime  = (g_runtime * 1000.0f) / qn;

		printf("%3d\t\t%.4f\t\t%.4f\t\t%.2f%%\t\t%.2f%%\n", top_k, g_ratio, 
			g_runtime, g_recall, g_fraction);
		fprintf(fp, "%d\t%f\t%f\t%f\t%f\n", top_k, g_ratio, g_runtime, 
			g_recall, g_fraction);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);
	delete lsh; 

	return 0;
}