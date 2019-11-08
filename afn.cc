#include <algorithm>
#include <cstring>
#include <sys/time.h>

#include "def.h"
#include "util.h"
#include "pri_queue.h"
#include "qdafn.h"
#include "drusilla_select.h"
#include "rqalsh.h"
#include "rqalsh_star.h"
#include "ml_rqalsh.h"
#include "afn.h"


// -----------------------------------------------------------------------------
int linear_scan(					// k-FN search of linear scan
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const float **data,					// data set
	const float **query,				// query set
	const Result **R, 					// truth set
	const char *out_path)				// output path
{
	char output_set[200];
	sprintf(output_set, "%slinear.out", out_path);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  c-k-AFN search
	// -------------------------------------------------------------------------
	printf("Top-k FN Search of Linear Scan:\n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		
		int top_k = TOPK[num];
		Result **result = new Result*[qn];
		for (int i = 0; i < qn; ++i) {
			result[i] = new Result[top_k];
		}
		k_fn_search(n, qn, d, top_k, data, query, result);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		for (int i = 0;i < qn; ++i) {
			g_recall += calc_recall(top_k, R[i], (const Result *) result[i]);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += R[i][j].key_ / R[i][j].key_;
			}
			g_ratio += ratio / top_k;
		}
		for (int i = 0; i < qn; ++i) {
			delete[] result[i]; result[i] = NULL;
		}
		delete[] result; result = NULL;

		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime = (g_runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, g_ratio, 
			g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, g_ratio, g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);
	
	return 0;
}

// -----------------------------------------------------------------------------
int ml_rqalsh(						// c-k-AFN search of ML-RQALSH
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	float ratio,						// approximation ratio
	const float **data,					// data set
	const float **query,				// query set
	const Result **R, 					// truth set
	const char *out_path)				// output path
{
	char output_set[200];
	sprintf(output_set, "%sml_rqalsh.out", out_path);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	ML_RQALSH* lsh = new ML_RQALSH(n, d, ratio, data);
	lsh->display();

	gettimeofday(&g_end_time, NULL);
	float indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f seconds\n\n", indexing_time);
	fprintf(fp, "Indexing Time = %f Seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  c-k-AFN search
	// -------------------------------------------------------------------------
	printf("Top-k FN Search of ML_RQALSH:\n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List *list = new MaxK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kfn(top_k, query[i], list);
			g_recall += calc_recall(top_k, R[i], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				if (list->ith_key(j) > FLOATZERO) {
					ratio += R[i][j].key_ / list->ith_key(j);
				}
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime = (g_runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, g_ratio, 
			g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, g_ratio, g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;

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
	const float **data,					// data set
	const float **query,				// query set
	const Result **R, 					// truth set
	const char *out_path)				// output path
{
	char output_set[200];
	sprintf(output_set, "%srqalsh_star.out", out_path);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	RQALSH_STAR* lsh = new RQALSH_STAR(n, d, L, M, ratio, data);
	lsh->display();
	
	gettimeofday(&g_end_time, NULL);
	float indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f seconds\n\n", indexing_time);
	fprintf(fp, "L = %d, M = %d\n", L, M);
	fprintf(fp, "Indexing Time = %f Seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  c-k-AFN search
	// -------------------------------------------------------------------------
	printf("Top-k FN Search of RQALSH*:\n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List *list = new MaxK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kfn(top_k, query[i], list);
			g_recall += calc_recall(top_k, R[i], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				if (list->ith_key(j) > FLOATZERO) {
					ratio += R[i][j].key_ / list->ith_key(j);
				}
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime = (g_runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, g_ratio, 
			g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, g_ratio, g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int rqalsh(							// c-k-AFN search of RQALSH
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	float ratio,						// approximation ratio
	const float **data,					// data set
	const float **query,				// query set
	const Result **R, 					// truth set
	const char *out_path)				// output path
{
	char output_set[200];
	sprintf(output_set, "%srqalsh.out", out_path);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	RQALSH* lsh = new RQALSH(n, d, ratio, data);
	lsh->display();
	
	gettimeofday(&g_end_time, NULL);
	float indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f seconds\n\n", indexing_time);
	fprintf(fp, "Indexing Time = %f Seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  c-k-AFN search
	// -------------------------------------------------------------------------
	printf("Top-k FN Search of RQALSH:\n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List *list = new MaxK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->kfn(top_k, query[i], list);
			g_recall += calc_recall(top_k, R[i], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				if (list->ith_key(j) > FLOATZERO) {
					ratio += R[i][j].key_ / list->ith_key(j);
				}
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime = (g_runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, g_ratio, 
			g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, g_ratio, g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int drusilla_select(				// c-k-AFN search of Drusilla-Select
	int   n,							// number of data objects
	int   qn,							// number of query points
	int   d,							// number of dimensions
	int   L,							// number of projections
	int   M,							// number of candidates
	const float **data,					// data set
	const float **query,				// query set
	const Result **R, 					// truth set
	const char *out_path)				// output path
{
	char output_set[200];
	sprintf(output_set, "%sdrusilla_select.out", out_path);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	Drusilla_Select *drusilla = new Drusilla_Select(n, d, L, M, data);
	drusilla->display();

	gettimeofday(&g_end_time, NULL);
	float indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f seconds\n\n", indexing_time);
	fprintf(fp, "L = %d, M = %d\n", L, M);
	fprintf(fp, "Indexing Time = %f Seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  c-k-AFN search of Drusilla-Select
	// -------------------------------------------------------------------------
	printf("Top-k FN Search of Drusilla-Select: \n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List *list = new MaxK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			drusilla->kfn(query[i], list);
			g_recall += calc_recall(top_k, R[i], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				if (list->ith_key(j) > FLOATZERO) {
					ratio += R[i][j].key_ / list->ith_key(j);
				}
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime = (g_runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, g_ratio, 
			g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, g_ratio, g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete drusilla; drusilla = NULL;

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
	const float **data,					// data set
	const float **query,				// query set
	const Result **R, 					// truth set
	const char *out_path)				// output path
{
	char output_set[200];
	sprintf(output_set, "%sqdafn.out", out_path);
	
	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  indexing 
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	QDAFN *hash = new QDAFN(n, d, L, M, 2, ratio, data);
	hash->display();

	gettimeofday(&g_end_time, NULL);
	float indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f seconds\n\n", indexing_time);
	fprintf(fp, "L = %d, M = %d\n", L, M);
	fprintf(fp, "Indexing Time = %f Seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  c-k-AFN search of QDAFN
	// -------------------------------------------------------------------------
	printf("Top-k FN Search of QDAFN: \n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");	
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MaxK_List *list = new MaxK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			hash->kfn(top_k, query[i], list);
			g_recall += calc_recall(top_k, R[i], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				if (list->ith_key(j) > FLOATZERO) {
					ratio += R[i][j].key_ / list->ith_key(j);
				}
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime = (g_runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, g_ratio, 
			g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, g_ratio, g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete hash; hash = NULL;

	return 0;
}