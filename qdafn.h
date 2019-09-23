#ifndef __QDAFN_H
#define __QDAFN_H

// -----------------------------------------------------------------------------
//  Ziggurat Method standard normal pseudorandom number generator code from 
//  George Marsaglia and Wai Wan Tsang (2000).
// 
//  "The Ziggurat Method for Generating Random Variables". Journal of
//  Statistical Software 5 (8).
// -----------------------------------------------------------------------------
extern long hz;
extern unsigned long jz,jsr;
extern unsigned long iz, kn[128], ke[256];
extern float wn[128],fn[128], we[256],fe[256];

static char *algoname[3] = { "By Value", "By Rank", "Query Dependent" };

#define SHR3 (jz=jsr, jsr^=(jsr<<13), jsr^=(jsr>>17), jsr^=(jsr<<5), jz+jsr)
#define UNI  (0.5 + (signed) SHR3 * 0.2328306e-9)
#define IUNI SHR3
#define RNOR (hz=SHR3, iz=hz&127, (abs(hz)<kn[iz])? hz*wn[iz] : nfix())
#define REXP (jz=SHR3, iz=jz&255, (jz <ke[iz])? jz*we[iz] : efix())

// -----------------------------------------------------------------------------
float nfix();
float efix();
void  zigset(unsigned long jsrseed);

// -----------------------------------------------------------------------------
//  data structure for QDAFN
// -----------------------------------------------------------------------------
struct RTA {
	int rank;
	int times_achieved;
};

// -----------------------------------------------------------------------------
struct PDIST_PAIR {
	int obj;
	union {
		float pdist;
		RTA rta;
	} u;
};

// -----------------------------------------------------------------------------
static int PDISTComp(               // compare func for quick (ascending)   
	const void *xv,                     // 1st element
	const void *yv);                    // 2nd element

// -----------------------------------------------------------------------------
static int PDISTCompDesc(           // compare func for quick (descending)
	const void *xv,                     // 1st element
	const void *yv);                    // 2nd element

// -----------------------------------------------------------------------------
static int RTAComp(                 // compare func for quick
	const void *xv,                     // 1st element
	const void *yv);                    // 2nd element

// -----------------------------------------------------------------------------
//  QDAFN: a hashing scheme for high-dimensional c-k-AFN search
// -----------------------------------------------------------------------------
class QDAFN {
public:
    QDAFN(                          // constructor
        int   n,						// cardinality
        int   d,						// dimensionality
        int   L,						// number of projections
        int   M,						// number of candidates
        int   algo,						// which algorithm
        float ratio,					// approximation ratio
        const float **data);	       	// data objects
    
    // -------------------------------------------------------------------------
    ~QDAFN();                       // destructor

    // -------------------------------------------------------------------------
    void display();                 // display parameters

    // -------------------------------------------------------------------------
    int kfn(                        // c-k-AFN search
        int   top_k,					// top-k value
	    const float *query,				// input query
	    MaxK_List *list);				// c-k-AFN results (return)

protected:
	int   n_pts_;					// cardinality
	int   dim_;				        // dimensionality
    int   L_;			            // number of random projections
	int   M_;				        // number of candidates
	int   algo_;	    	    	// which algorithm
    float appr_ratio_;              // approximation ratio
    const float **data_;			// data objects

    float *proj_;			        // projection vectors
	PDIST_PAIR *pdp_;				// projected info after random projection

	// -------------------------------------------------------------------------
    int bulkload();                 // build index    
};

#endif // __QDAFN_H