#include "headers.h"

// -----------------------------------------------------------------------------
//  use Box-Muller transform to generate a r.v. from Gaussian(mean, sigma)
//  standard Gaussian distr. is Gaussian(0, 1), where mean = 0 and sigma = 1
// -----------------------------------------------------------------------------
float gaussian(						// r.v. from N(mean, sigma)
	float mu,							// mean (location)
	float sigma)						// stanard deviation (scale > 0)
{
	// assert(sigma > 0.0f);
	float u1 = -1.0f;
	float u2 = -1.0f;
	do {
		u1 = uniform(0.0f, 1.0f);
	} while (u1 < FLOATZERO);
	u2 = uniform(0.0f, 1.0f);

	return mu + sigma * sqrt(-2.0f * log(u1)) * cos(2.0f * PI * u2);
	// return mu + sigma * sqrt(-2.0f * log(u1)) * sin(2.0f * PI * u2);
}

// -----------------------------------------------------------------------------
float gaussian_cdf(					// cdf of N(0, 1) in range (-inf, x]
	float x,							// integral border
	float step)							// step increment
{
	float ret = 0.0f;
	for (float i = -10.0f; i < x; i += step) {
		ret += step * gaussian_pdf(i);
	}
	return ret;
}

// -----------------------------------------------------------------------------
float new_gaussian_cdf(				// cdf of N(0, 1) in range [-x, x]
	float x,							// integral border (x > 0)
	float step)							// step increment
{
	// assert(x > 0.0f);
	float ret = 0.0f;
	for (float i = -x; i <= x; i += step) {
		ret += step * gaussian_pdf(i);
	}
	return ret;
}

// -----------------------------------------------------------------------------
//  query-oblivious and query-aware collision probability under gaussian
//  distribution and cauchy distribution
// -----------------------------------------------------------------------------
float orig_gaussian_prob(			// calc original gaussian probability
	float x)							// x = w / r
{
	float norm = gaussian_cdf(-x, 0.001F);
	float tmp  = 2.0F * (1.0F - exp(-x * x / 2.0F)) / (sqrt(2.0F * PI) * x);

	return 1.0F - 2.0F * norm - tmp;
}

// -----------------------------------------------------------------------------
float new_gaussian_prob(			// calc new gaussian probability
	float x)							// x = w / (2 * r)
{
	return new_gaussian_cdf(x, 0.001F);
}

// -----------------------------------------------------------------------------
//  probability vs. w for a fixed ratio c
// -----------------------------------------------------------------------------
void prob_of_gaussian()				// curve of p1, p2 vs. w under gaussian
{
	printf("probability vs. w for c = {2.0, 3.0} under gaussian\n");
	float c[2] = { 2.0f, 3.0f };
	float orig_p1, orig_p2, new_p1, new_p2;

	for (int i = 0; i < 2; ++i) {
		printf("c = %.1f\n", c[i]);

		for (float w = 0.5f; w < 10.1f; w += 0.5f) {
			orig_p1 = orig_gaussian_prob(w);
			orig_p2 = orig_gaussian_prob(w / c[i]);

			new_p1 = new_gaussian_prob(w / 2.0f);
			new_p2 = new_gaussian_prob(w / (2.0f * c[i]));

			printf("%.1f\t%.4f\t%.4f\t%.4f\t%.4f\n",
				w, orig_p1, orig_p2, new_p1, new_p2);
		}
		printf("\n");
	}
}

// -----------------------------------------------------------------------------
//  the difference (p1 - p2) vs. w for a fixed ratio
// -----------------------------------------------------------------------------
void diff_prob_of_gaussian()		// curve of p1 - p2 vs. w under gaussian
{
	printf("prob of diff vs. w for c = {2.0, 3.0} under gaussian\n");
	float c[2] = { 2.0f, 3.0f };
	float orig_p1, orig_p2, new_p1, new_p2, orig_diff, new_diff;

	for (int i = 0; i < 2; ++i) {
		printf("c = %.1f\n", c[i]);

		for (float w = 0.5f; w < 10.1f; w += 0.5f) {
			orig_p1 = orig_gaussian_prob(w);
			orig_p2 = orig_gaussian_prob(w / c[i]);
			orig_diff = orig_p1 - orig_p2;

			new_p1 = new_gaussian_prob(w / 2.0f);
			new_p2 = new_gaussian_prob(w / (2.0f * c[i]));
			new_diff = new_p1 - new_p2;

			printf("%.1f\t%.4f\t%.4f\n", w, orig_diff, new_diff);
		}
		printf("\n");
	}
}

// -----------------------------------------------------------------------------
//  rho = log(1/p1) / log(1/p2) vs. w for a fixed ratio c
// -----------------------------------------------------------------------------
void rho_of_gaussian()				// curve of rho vs. w under gaussian
{
	printf("rho vs. w for c = {2.0, 3.0} under gaussian\n");
	float c[2] = { 2.0f, 3.0f };
	float orig_p1, orig_p2, new_p1, new_p2, orig_rho, new_rho;

	for (int i = 0; i < 2; ++i) {
		printf("c = %.1f\n", c[i]);

		for (float w = 0.5f; w < 10.1f; w += 0.5f) {
			orig_p1 = orig_gaussian_prob(w);
			orig_p2 = orig_gaussian_prob(w / c[i]);
			orig_rho = log(1.0f / orig_p1) / log(1.0f / orig_p2);

			new_p1 = new_gaussian_prob(w / 2.0f);
			new_p2 = new_gaussian_prob(w / (2.0f * c[i]));
			new_rho = log(1.0f / new_p1) / log(1.0f / new_p2);

			printf("%.1f\t%.4f\t%.4f\t%.4f\n", w, orig_rho, new_rho, 1.0f / c[i]);
		}
		printf("\n");
	}
}
