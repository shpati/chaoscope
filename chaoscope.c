/* chaoscope.c -- a pure-C "Chaoscope"-style strange attractor renderer.
 *
 * Reimplementation of the idea behind Nicolas Desprez's Chaoscope (Windows):
 * integrate a chaotic dynamical system, rotate/zoom the 3D view, cycle a
 * colour palette, and render the attractor as a point cloud with density
 * accumulation.  Output is a 24-bit BMP image (and optional animated frame
 * sequences).  No third-party libraries: only standard C99 (<stdio.h>,
 * <stdlib.h>, <math.h>, <string.h>, <time.h>).  Compiles on Windows, Linux,
 * macOS with any C99 compiler:
 *
 *   gcc -O2 -o chaoscope chaoscope.c -lm
 *   ./chaoscope --system lorenz --output attractor.bmp
 *
 * Equation sources (checked against references, then numerically: every
 * default below stays bounded and has a positive largest Lyapunov exponent
 * with the listed initial condition and step):
 *   Chen, Halvorsen, Rucklidge  arXiv:2601.21117 (Attractors section)
 *   Dadras, Aizawa, Burke-Shaw  ChaoticDynamicalSystemLibrary (Julia) docs
 *   Four-wing                   arXiv:2507.07577, Eq. 3.1
 *   Thomas                      cyclic form, chaos onset b ~ 0.208186
 *   Burke-Shaw parameters       s = 10, v = 4.272 (P. Bourke)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* Program version, kept in sync with the Win32 GUI front end. */
#define CHAOSCOPE_VERSION "v1.0.0 by Shpati Koleka"

/* ------------------------------------------------------------------ */
/* Small helpers                                                       */
/* ------------------------------------------------------------------ */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DEG2RAD(d) ((d) * (M_PI / 180.0))

#define MAX_PARAMS  8
#define MAX_DIM     16384    /* largest accepted image width/height     */
#define NORM_PCT    0.998    /* density percentile mapped to full colour */
#define STATE_LIMIT 1e6      /* |state| beyond this counts as divergence */
/* A unit sphere seen through the perspective below (d = 1 + 0.5*z) projects
 * to at most 2/sqrt(3) ~ 1.155 times the focal length.  FIT scales the focal
 * length so that worst case lands just inside the image (zoom = 1). */
#define FIT         0.82

static double clamp01(double v) { return v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v); }

/* NaN-safe sanity check (NaN compares false with everything). */
static int sane(double v) { return v > -STATE_LIMIT && v < STATE_LIMIT; }

static unsigned long long rng_state; /* xorshift64* PRNG (portable) */

/* Seeds the generator and returns the seed actually used (0 => time). */
static unsigned long long rng_seed(unsigned long long s) {
    if (s == 0)
        s = (unsigned long long)(time(NULL) & 0x7fffffff) * 2654435761ULL + 12345ULL;
    if (s == 0) s = 88172645463325252ULL;
    rng_state = s;
    return s;
}

static unsigned long long rng_next(void) {
    unsigned long long x = rng_state;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    rng_state = x;
    return x * 2685821657736338717ULL;
}

static double rng_double(void) { /* in [0,1) */
    return (double)(rng_next() >> 11) / 9007199254740992.0;
}

/* ------------------------------------------------------------------ */
/* Chaotic systems                                                     */
/*                                                                     */
/* Each system has a name, up to MAX_PARAMS parameters, a nominal       */
/* initial condition, a good default step, and a derivative function.   */
/* ------------------------------------------------------------------ */

typedef struct {
    const char *name;
    const char *desc;
    int    nparams;                 /* how many parameters are used      */
    double def[MAX_PARAMS];         /* default parameters                */
    const char *pnames[MAX_PARAMS]; /* parameter display names           */
    double ic[3];                   /* nominal initial condition (x,y,z) */
    double dt;                      /* default integration step          */
    void (*deriv)(const double *p, const double *s, double *d); /* s->d  */
} SystemDef;

/* --- Lorenz: p = sigma, beta, rho --- */
static void lorenz_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    double sigma = p[0], beta = p[1], rho = p[2];
    d[0] = sigma * (y - x);
    d[1] = x * (rho - z) - y;
    d[2] = x * y - beta * z;
}

/* --- Rossler: p = a, b, c --- */
static void rossler_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    double a = p[0], b = p[1], c = p[2];
    d[0] = -y - z;
    d[1] = x + a * y;
    d[2] = b + z * (x - c);
}

/* --- Thomas (cyclically symmetric): p = b --- */
static void thomas_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    double b = p[0];
    d[0] = sin(y) - b * x;
    d[1] = sin(z) - b * y;
    d[2] = sin(x) - b * z;
}

/* --- Aizawa: p = a, b, c, d, e, f --- */
static void aizawa_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    double a = p[0], b = p[1], c = p[2], w = p[3], e = p[4], f = p[5];
    d[0] = (z - b) * x - w * y;
    d[1] = w * x + (z - b) * y;
    d[2] = c + a * z - z * z * z / 3.0 - (x * x + y * y) * (1.0 + e * z)
           + f * z * x * x * x;
}

/* --- Halvorsen: p = a --- */
static void halvorsen_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    double a = p[0];
    d[0] = -a * x - 4.0 * y - 4.0 * z - y * y;
    d[1] = -a * y - 4.0 * z - 4.0 * x - z * z;
    d[2] = -a * z - 4.0 * x - 4.0 * y - x * x;
}

/* --- Chen: p = a, b, c --- */
static void chen_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    double a = p[0], b = p[1], c = p[2];
    d[0] = a * (y - x);
    d[1] = (c - a) * x - x * z + c * y;
    d[2] = x * y - b * z;
}

/* --- Lu (Lu & Chen): p = a, b, c --- */
static void lu_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    double a = p[0], b = p[1], c = p[2];
    d[0] = a * (y - x);
    d[1] = -x * z + c * y;
    d[2] = x * y - b * z;
}

/* --- Rucklidge: p = kappa, lambda --- */
static void rucklidge_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    double kappa = p[0], lambda = p[1];
    d[0] = -kappa * x + lambda * y - y * z;
    d[1] = x;
    d[2] = -z + y * y;
}

/* --- Four-wing: p = a, b, c, d, e, f --- */
static void fourwing_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    double a = p[0], b = p[1], c = p[2], w = p[3], e = p[4], f = p[5];
    d[0] = a * x + c * y * z;
    d[1] = b * x + w * y - x * z;
    d[2] = e * z + f * x * y;
}

/* --- Dadras: p = p, o, r, c, e --- */
static void dadras_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    double pp = p[0], o = p[1], r = p[2], c = p[3], e = p[4];
    d[0] = y - pp * x + o * y * z;
    d[1] = r * y - x * z + z;
    d[2] = c * x * y - e * z;
}

/* --- Burke-Shaw: p = s, v --- */
static void burkeshaw_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    double sg = p[0], v = p[1];
    d[0] = -sg * (x + y);
    d[1] = -y - sg * x * z;
    d[2] = sg * x * y + v;
}

/* --- Sprott's 19 simplest chaotic flows (cases A..S).  A..E have no free
 *     parameters; F..S use Sprott's (a, b) parametrisation. --- */
static void sprott_a_deriv(const double *p, const double *s, double *d) {  /* Nose-Hoover */
    double x = s[0], y = s[1], z = s[2]; (void)p; (void)x;
    d[0] = y;
    d[1] = -x + y * z;
    d[2] = 1.0 - y * y;
}
static void sprott_b_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2]; (void)p;
    d[0] = y * z;
    d[1] = x - y;
    d[2] = 1.0 - x * y;
}
static void sprott_c_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2]; (void)p;
    d[0] = y * z;
    d[1] = x - y;
    d[2] = 1.0 - x * x;
}
static void sprott_d_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2]; (void)p;
    d[0] = -y;
    d[1] = x + z;
    d[2] = x * z + 3.0 * y * y;
}
static void sprott_e_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2]; (void)p;
    d[0] = y * z;
    d[1] = x * x - y;
    d[2] = 1.0 - 4.0 * x;
}
static void sprott_f_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    d[0] = y + z;
    d[1] = -x + p[0] * y;
    d[2] = x * x - p[1] * z;
}
static void sprott_g_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    d[0] = p[0] * x + p[1] * z;
    d[1] = x * z - y;
    d[2] = -x + y;
}
static void sprott_h_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    d[0] = -y + z * z;
    d[1] = x + p[0] * y;
    d[2] = x - p[1] * z;
}
static void sprott_i_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    d[0] = -p[0] * y;
    d[1] = x + z;
    d[2] = x + y * y - p[1] * z;
}
static void sprott_j_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    d[0] = p[0] * z;
    d[1] = -p[1] * y + z;
    d[2] = -x + y + y * y;
}
static void sprott_k_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    d[0] = x * y - z;
    d[1] = x - y;
    d[2] = p[0] * x + p[1] * z;
}
static void sprott_l_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    d[0] = y + p[0] * z;
    d[1] = 0.9 * x * x - y;
    d[2] = p[1] - x;
}
static void sprott_m_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    d[0] = -z;
    d[1] = -x * x - y;
    d[2] = p[0] + p[1] * x + y;
}
static void sprott_n_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    d[0] = -p[0] * y;
    d[1] = x + z * z;
    d[2] = 1.0 + y - p[1] * z;
}
static void sprott_o_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    d[0] = y;
    d[1] = x - z;
    d[2] = p[0] * x + x * z + p[1] * y;
}
static void sprott_p_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    d[0] = p[0] * y + z;
    d[1] = -x + y * y;
    d[2] = x + p[1] * y;
}
static void sprott_q_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    d[0] = -z;
    d[1] = x - y;
    d[2] = p[0] * x + y * y + p[1] * z;
}
static void sprott_r_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    d[0] = p[0] - y;
    d[1] = p[1] + z;
    d[2] = x * y - z;
}
static void sprott_s_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    d[0] = -x - p[0] * y;
    d[1] = x + z * z;
    d[2] = p[1] + x;
}

/* --- Arneodo / Coullet jerk system: p = a, b, c, d --- */
static void arneodo_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    d[0] = y;
    d[1] = z;
    d[2] = -p[0] * x - p[1] * y - p[2] * z + p[3] * x * x * x;
}

/* --- Lorenz-84 / Hadley circulation: p = a, b, F, G --- */
static void lorenz84_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    double a = p[0], b = p[1], f = p[2], g = p[3];
    d[0] = -y * y - z * z - a * x + a * f;
    d[1] = x * y - b * x * z - y + g;
    d[2] = b * x * y + x * z - z;
}

/* --- Lorenz mod 2: p = a, b, c, d --- */
static void lorenzmod2_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    double a = p[0], b = p[1], c = p[2], w = p[3];
    d[0] = -a * x + y * y - z * z + a * c;
    d[1] = x * (y - b * z) + w;
    d[2] = -z + x * (b * y + z);
}

/* --- Chua double scroll: p = alpha, beta, m0, m1 --- */
static void chua_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    double ramp = p[3] * x + 0.5 * (p[2] - p[3]) * (fabs(x + 1.0) - fabs(x - 1.0));
    d[0] = p[0] * (y - x - ramp);
    d[1] = x - y + z;
    d[2] = -p[1] * y;
}

/* --- Chen-Lee: p = a, b, c --- */
static void chenlee_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    d[0] = p[0] * x - y * z;
    d[1] = p[1] * y + x * z;
    d[2] = p[2] * z + x * y / 3.0;
}

/* --- Dequan Li: p = a, c, d, eps, k, f --- */
static void dequanli_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    double a = p[0], c = p[1], w = p[2], eps = p[3], k = p[4], f = p[5];
    d[0] = a * (y - x) + w * x * z;
    d[1] = k * x + f * y - x * z;
    d[2] = c * z + x * y - eps * x * x;
}

/* --- Three-scroll unified system: p = a, b, c, d, e --- */
static void threescroll_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    d[0] = p[0] * (y - x) + p[3] * x * z;
    d[1] = p[2] * y - x * z;
    d[2] = p[1] * z + x * y - p[4] * x * x;
}

/* --- Finance (y shifted by 1/b): p = a, b, c --- */
static void finance_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    d[0] = (1.0 / p[1] - p[0]) * x + z + x * y;
    d[1] = -p[1] * y - x * x;
    d[2] = -x - p[2] * z;
}

/* --- Rabinovich-Fabrikant: p = alpha, gamma --- */
static void rabfab_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    d[0] = y * (z - 1.0 + x * x) + p[1] * x;
    d[1] = x * (3.0 * z + 1.0 - x * x) + p[1] * y;
    d[2] = -2.0 * z * (p[0] + x * y);
}

/* --- Newton-Leipnik: p = alpha, beta --- */
static void leipnik_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    d[0] = -p[0] * x + y + 10.0 * y * z;
    d[1] = -x - 0.4 * y + 5.0 * x * z;
    d[2] = p[1] * z - 5.0 * x * y;
}

/* --- Rikitake dynamo: p = mu, a --- */
static void rikitake_deriv(const double *p, const double *s, double *d) {
    double x = s[0], y = s[1], z = s[2];
    d[0] = -p[0] * x + y * z;
    d[1] = -p[0] * y + (z - p[1]) * x;
    d[2] = 1.0 - x * y;
}

static const SystemDef SYSTEMS[] = {
    { "lorenz",    "Lorenz butterfly", 3,
      {10.0, 8.0 / 3.0, 28.0}, {"sigma", "beta", "rho"},
      {0.1, 0.0, 0.0}, 0.004, lorenz_deriv },
    { "rossler",   "Rossler spiral chaos", 3,
      {0.2, 0.2, 5.7}, {"a", "b", "c"},
      {0.1, 0.0, 0.0}, 0.02, rossler_deriv },
    { "thomas",    "Thomas cyclically symmetric (chaotic from b ~ 0.208186)", 1,
      {0.208186}, {"b"},
      {0.1, 0.0, 0.0}, 0.05, thomas_deriv },
    { "aizawa",    "Aizawa toroidal chaos", 6,
      {0.95, 0.7, 0.6, 3.5, 0.25, 0.1}, {"a", "b", "c", "d", "e", "f"},
      {0.1, 0.0, 0.0}, 0.01, aizawa_deriv },
    { "halvorsen", "Halvorsen cyclic chaos", 1,
      {1.4}, {"a"},
      {-5.0, 0.0, 0.0}, 0.005, halvorsen_deriv },
    { "chen",      "Chen system", 3,
      {35.0, 3.0, 28.0}, {"a", "b", "c"},
      {0.1, 0.0, 0.0}, 0.005, chen_deriv },
    { "lu",        "Lu system (Lu & Chen)", 3,
      {36.0, 3.0, 20.0}, {"a", "b", "c"},
      {0.1, 0.1, 0.1}, 0.005, lu_deriv },
    { "rucklidge", "Rucklidge double convection", 2,
      {2.0, 6.7}, {"kappa", "lambda"},
      {1.0, 0.0, 4.5}, 0.01, rucklidge_deriv },
    { "fourwing",  "Four-wing attractor", 6,
      {0.2, -0.01, 1.0, -0.4, -1.0, -1.0}, {"a", "b", "c", "d", "e", "f"},
      {1.0, 1.0, 1.0}, 0.05, fourwing_deriv },
    { "dadras",    "Dadras multi-scroll", 5,
      {3.0, 2.7, 1.7, 2.0, 9.0}, {"p", "o", "r", "c", "e"},
      {1.1, 2.1, -2.0}, 0.005, dadras_deriv },
    { "burkeshaw", "Burke-Shaw", 2,
      {10.0, 4.272}, {"s", "v"},
      {0.6, 0.5, 0.5}, 0.005, burkeshaw_deriv },
    { "sprott_a",  "Sprott A (Nose-Hoover, conservative)", 0,
      {0}, {0},
      {0.0, 5.0, 0.0}, 0.02, sprott_a_deriv },
    { "sprott_b",  "Sprott B", 0,
      {0}, {0},
      {0.1, 0.1, 0.1}, 0.02, sprott_b_deriv },
    { "sprott_c",  "Sprott C", 0,
      {0}, {0},
      {0.1, 0.1, 0.1}, 0.02, sprott_c_deriv },
    { "sprott_d",  "Sprott D", 0,
      {0}, {0},
      {0.11, 0.08, 0.14}, 0.02, sprott_d_deriv },
    { "sprott_e",  "Sprott E", 0,
      {0}, {0},
      {0.1, 0.1, 0.1}, 0.02, sprott_e_deriv },
    { "sprott_f",  "Sprott F", 2,
      {0.5, 1.0}, {"a", "b"},
      {0.1, 0.1, 0.1}, 0.02, sprott_f_deriv },
    { "sprott_g",  "Sprott G", 2,
      {0.4, 1.0}, {"a", "b"},
      {0.1, 0.1, 0.1}, 0.02, sprott_g_deriv },
    { "sprott_h",  "Sprott H", 2,
      {0.5, 1.0}, {"a", "b"},
      {0.1, 0.1, 0.1}, 0.02, sprott_h_deriv },
    { "sprott_i",  "Sprott I", 2,
      {0.2, 1.0}, {"a", "b"},
      {0.1, 0.1, 0.1}, 0.02, sprott_i_deriv },
    { "sprott_j",  "Sprott J", 2,
      {2.0, 2.0}, {"a", "b"},
      {0.1, 0.1, 0.1}, 0.02, sprott_j_deriv },
    { "sprott_k",  "Sprott K", 2,
      {1.0, 0.3}, {"a", "b"},
      {0.1, 0.1, 0.1}, 0.02, sprott_k_deriv },
    { "sprott_l",  "Sprott L", 2,
      {3.9, 1.0}, {"a", "b"},
      {0.1, 0.1, 0.1}, 0.02, sprott_l_deriv },
    { "sprott_m",  "Sprott M", 2,
      {1.7, 1.7}, {"a", "b"},
      {0.1, 0.1, 0.1}, 0.02, sprott_m_deriv },
    { "sprott_n",  "Sprott N", 2,
      {2.0, 2.0}, {"a", "b"},
      {0.1, 0.1, 0.1}, 0.02, sprott_n_deriv },
    { "sprott_o",  "Sprott O", 2,
      {1.0, 2.7}, {"a", "b"},
      {-0.19, 0.24, -0.12}, 0.02, sprott_o_deriv },
    { "sprott_p",  "Sprott P", 2,
      {2.7, 1.0}, {"a", "b"},
      {0.1, 0.1, 0.1}, 0.02, sprott_p_deriv },
    { "sprott_q",  "Sprott Q", 2,
      {3.1, 0.5}, {"a", "b"},
      {0.1, 0.1, 0.1}, 0.02, sprott_q_deriv },
    { "sprott_r",  "Sprott R", 2,
      {0.9, 0.4}, {"a", "b"},
      {0.1, 0.1, 0.1}, 0.02, sprott_r_deriv },
    { "sprott_s",  "Sprott S", 2,
      {4.0, 1.0}, {"a", "b"},
      {0.1, 0.1, 0.1}, 0.02, sprott_s_deriv },
    { "arneodo",   "Arneodo jerk system", 4,
      {-5.5, 3.5, 1.0, -1.0}, {"a", "b", "c", "d"},
      {0.1, 0.1, 0.1}, 0.02, arneodo_deriv },
    { "coullet",   "Coullet jerk system (Arneodo form)", 4,
      {-0.8, 1.1, 0.45, -1.0}, {"a", "b", "c", "d"},
      {0.1, 0.1, 0.1}, 0.02, arneodo_deriv },
    { "lorenz84",  "Lorenz-84 / Hadley circulation", 4,
      {0.25, 4.0, 8.0, 1.0}, {"a", "b", "F", "G"},
      {0.1, 0.1, 0.1}, 0.01, lorenz84_deriv },
    { "lorenzmod2","Lorenz mod 2", 4,
      {0.9, 5.0, 9.9, 1.0}, {"a", "b", "c", "d"},
      {0.1, 0.1, 0.1}, 0.01, lorenzmod2_deriv },
    { "chua",      "Chua double scroll", 4,
      {15.6, 28.0, -8.0 / 7.0, -5.0 / 7.0}, {"alpha", "beta", "m0", "m1"},
      {0.7, 0.0, 0.0}, 0.005, chua_deriv },
    { "chenlee",   "Chen-Lee", 3,
      {5.0, -10.0, -3.8}, {"a", "b", "c"},
      {0.1, 0.1, 0.1}, 0.01, chenlee_deriv },
    { "dequanli",  "Dequan Li three-scroll", 6,
      {40.0, 11.0 / 6.0, 0.16, 0.65, 55.0, 20.0}, {"a", "c", "d", "eps", "k", "f"},
      {0.349, 0.0, -0.16}, 0.0005, dequanli_deriv },
    { "3scroll",   "Three-scroll unified system", 5,
      {40.0, 0.833, 20.0, 0.5, 0.65}, {"a", "b", "c", "d", "e"},
      {0.1, 0.1, 0.1}, 0.0005, threescroll_deriv },
    { "finance",   "Finance system (y shifted by 1/b)", 3,
      {3.0, 0.1, 1.0}, {"a", "b", "c"},
      {1.0, -8.0, 0.5}, 0.01, finance_deriv },
    { "rabfab",    "Rabinovich-Fabrikant", 2,
      {0.14, 0.10}, {"alpha", "gamma"},
      {-1.0, 0.0, 0.5}, 0.01, rabfab_deriv },
    { "leipnik",   "Newton-Leipnik", 2,
      {0.4, 0.175}, {"alpha", "beta"},
      {0.349, 0.0, -0.16}, 0.01, leipnik_deriv },
    { "rikitake",  "Rikitake dynamo", 2,
      {2.0, 5.0}, {"mu", "a"},
      {0.1, 0.1, 0.1}, 0.02, rikitake_deriv },
};
#define NSYS ((int)(sizeof(SYSTEMS) / sizeof(SYSTEMS[0])))

/* ------------------------------------------------------------------ */
/* RK4 integrator                                                      */
/* ------------------------------------------------------------------ */

typedef struct {
    int    w, h;
    int    iterations;   /* total ODE steps to render                 */
    int    burn;         /* steps to skip before accumulating         */
    int    bounds_n;     /* steps used to estimate the bounding box   */
    double dt;
    double yaw, pitch, roll; /* degrees */
    double zoom;
    double panx, pany;
    int    splat;        /* splat radius (0..4)                       */
    double p[MAX_PARAMS];
    const SystemDef *sys;
} RenderParams;

typedef struct {
    double p[MAX_PARAMS];
    double s[3];    /* current state */
    double dt;      /* integration step */
    const SystemDef *sys;
} Solver;

static void solver_init(Solver *sv, const RenderParams *rp) {
    int k;
    sv->sys = rp->sys;
    memcpy(sv->p, rp->p, sizeof(sv->p));
    sv->dt = rp->dt;
    /* small random nudge so different seeds follow different trajectories */
    for (k = 0; k < 3; k++)
        sv->s[k] = rp->sys->ic[k] + (rng_double() - 0.5) * 0.1;
}

static void solver_step(Solver *sv) {
    double k1[3], k2[3], k3[3], k4[3], tmp[3], st[3];
    const double h = sv->dt;
    const double *p = sv->p;

    sv->sys->deriv(p, sv->s, k1);
    tmp[0] = sv->s[0] + 0.5 * h * k1[0];
    tmp[1] = sv->s[1] + 0.5 * h * k1[1];
    tmp[2] = sv->s[2] + 0.5 * h * k1[2];
    sv->sys->deriv(p, tmp, k2);
    tmp[0] = sv->s[0] + 0.5 * h * k2[0];
    tmp[1] = sv->s[1] + 0.5 * h * k2[1];
    tmp[2] = sv->s[2] + 0.5 * h * k2[2];
    sv->sys->deriv(p, tmp, k3);
    tmp[0] = sv->s[0] + h * k3[0];
    tmp[1] = sv->s[1] + h * k3[1];
    tmp[2] = sv->s[2] + h * k3[2];
    sv->sys->deriv(p, tmp, k4);

    st[0] = sv->s[0] + (h / 6.0) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]);
    st[1] = sv->s[1] + (h / 6.0) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]);
    st[2] = sv->s[2] + (h / 6.0) * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]);
    sv->s[0] = st[0]; sv->s[1] = st[1]; sv->s[2] = st[2];
}

/* One step with divergence protection.  Returns 1 on a normal step, 0 if
 * the state blew up and was restarted from the initial condition, and -1
 * once the system has diverged too many times to be worth continuing. */
static int step_checked(Solver *sv, const RenderParams *rp, int *restarts) {
    solver_step(sv);
    if (sane(sv->s[0]) && sane(sv->s[1]) && sane(sv->s[2])) return 1;
    if (++(*restarts) > 50) return -1;
    solver_init(sv, rp);
    return 0;
}

/* ------------------------------------------------------------------ */
/* Palettes: cyclic colour gradients                                   */
/* ------------------------------------------------------------------ */

typedef struct { double pos; unsigned char r, g, b; } PaletteStop;

static const PaletteStop PAL_RAINBOW[] = {
    {0.00,255, 40, 40},{0.20,255,220, 60},{0.40, 60,255, 60},
    {0.60, 40,200,255},{0.80,120, 60,255},{1.00,255, 40, 40}
};
static const PaletteStop PAL_FIRE[] = {
    {0.00,  0,  0,  0},{0.25, 80,  0,  0},{0.50,255, 60,  0},
    {0.75,255,200, 60},{1.00,255,255,220}
};
static const PaletteStop PAL_ICE[] = {
    {0.00,  0,  0, 40},{0.30,  0, 60,120},{0.60, 40,140,220},
    {0.80,150,220,255},{1.00,255,255,255}
};
static const PaletteStop PAL_PLASMA[] = {
    {0.00, 40,  0, 90},{0.25,150, 20,200},{0.50, 60,120,255},
    {0.75, 40,230,220},{1.00, 40,  0, 90}
};
static const PaletteStop PAL_CHAOSCOPE[] = {
    {0.00,  0,  0, 90},{0.14, 60, 40,200},{0.28, 40,200,255},
    {0.42, 60,255, 90},{0.56,255,220, 60},{0.70,255, 90, 40},
    {0.84,230, 40,120},{1.00,  0,  0, 90}
};

static const struct { const char *name; const PaletteStop *stops; int n; } PALETTES[] = {
    { "rainbow",   PAL_RAINBOW,   (int)(sizeof(PAL_RAINBOW)/sizeof(PAL_RAINBOW[0])) },
    { "fire",      PAL_FIRE,      (int)(sizeof(PAL_FIRE)/sizeof(PAL_FIRE[0])) },
    { "ice",       PAL_ICE,       (int)(sizeof(PAL_ICE)/sizeof(PAL_ICE[0])) },
    { "plasma",    PAL_PLASMA,    (int)(sizeof(PAL_PLASMA)/sizeof(PAL_PLASMA[0])) },
    { "chaoscope", PAL_CHAOSCOPE, (int)(sizeof(PAL_CHAOSCOPE)/sizeof(PAL_CHAOSCOPE[0])) },
};
#define NPAL ((int)(sizeof(PALETTES)/sizeof(PALETTES[0])))

static unsigned char lut[3][256]; /* r,g,b lookup tables */
static double lut_gamma = 1.0;    /* output gamma applied to the table */

/* Fill the lookup table from a stop list, rotated cyclically by `phase`. */
static void build_lut(const PaletteStop *stops, int n, double phase) {
    int i;
    for (i = 0; i < 256; i++) {
        double v = (double)i / 255.0 + phase;
        double a, b, f;
        int seg = 0;
        v -= floor(v);                          /* cyclic wrap in [0,1) */
        while (seg < n - 2 && v > stops[seg + 1].pos) seg++;
        a = stops[seg].pos;
        b = stops[seg + 1].pos;
        f = (b > a) ? clamp01((v - a) / (b - a)) : 0.0;
        lut[0][i] = (unsigned char)(stops[seg].r + (stops[seg + 1].r - stops[seg].r) * f + 0.5);
        lut[1][i] = (unsigned char)(stops[seg].g + (stops[seg + 1].g - stops[seg].g) * f + 0.5);
        lut[2][i] = (unsigned char)(stops[seg].b + (stops[seg + 1].b - stops[seg].b) * f + 0.5);
        if (lut_gamma != 1.0) {            /* >1 lifts dim regions, <1 deepens them */
            int c;
            for (c = 0; c < 3; c++)
                lut[c][i] = (unsigned char)(255.0 * pow(lut[c][i] / 255.0, 1.0 / lut_gamma) + 0.5);
        }
    }
}

/* ------------------------------------------------------------------ */
/* BMP output (24-bit, little-endian, row-padded)                      */
/* ------------------------------------------------------------------ */

static void put_le32(unsigned char *p, unsigned long v) {
    p[0] = (unsigned char)( v        & 0xff);
    p[1] = (unsigned char)((v >>  8) & 0xff);
    p[2] = (unsigned char)((v >> 16) & 0xff);
    p[3] = (unsigned char)((v >> 24) & 0xff);
}

static int write_bmp(const char *path, int w, int h, const unsigned char *rgb) {
    FILE *f;
    unsigned long rowsize = ((unsigned long)w * 3UL + 3UL) & ~3UL; /* 4-byte rows */
    unsigned long datasize = rowsize * (unsigned long)h;
    unsigned char hdr[54];
    unsigned char *buf;
    int i, ok = 1;

    memset(hdr, 0, sizeof(hdr));
    hdr[0] = 'B'; hdr[1] = 'M';
    put_le32(hdr + 2, 54UL + datasize);   /* file size           */
    put_le32(hdr + 10, 54UL);             /* pixel data offset   */
    put_le32(hdr + 14, 40UL);             /* info header size    */
    put_le32(hdr + 18, (unsigned long)w);
    put_le32(hdr + 22, (unsigned long)h); /* positive => bottom-up */
    hdr[26] = 1; hdr[28] = 24;            /* planes=1, bpp=24    */
    put_le32(hdr + 34, datasize);         /* image size          */

    buf = (unsigned char *)malloc((size_t)rowsize);
    if (!buf) return 0;

    f = fopen(path, "wb");
    if (!f) { free(buf); return 0; }

    if (fwrite(hdr, 1, 54, f) != 54) ok = 0;
    for (i = h - 1; i >= 0 && ok; i--) {   /* BMP rows stored bottom-up */
        int x;
        memset(buf, 0, (size_t)rowsize);
        for (x = 0; x < w; x++) {
            const unsigned char *px = rgb + ((size_t)i * (size_t)w + (size_t)x) * 3;
            buf[x * 3 + 0] = px[2]; /* B */
            buf[x * 3 + 1] = px[1]; /* G */
            buf[x * 3 + 2] = px[0]; /* R */
        }
        if (fwrite(buf, 1, (size_t)rowsize, f) != (size_t)rowsize) ok = 0;
    }
    if (fclose(f) != 0) ok = 0;
    free(buf);
    return ok;
}

/* ------------------------------------------------------------------ */
/* Rendering                                                           */
/* ------------------------------------------------------------------ */

static int cmp_double(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

/* Add a (2r+1)^2 soft-disk kernel, scaled by wt, centred on pixel (ix,iy). */
static void splat(double *dens, int w, int h, int ix, int iy, int r,
                  const double *kern, double wt) {
    int dx, dy, kn = 2 * r + 1;
    for (dy = -r; dy <= r; dy++) {
        int yy = iy + dy;
        if (yy < 0 || yy >= h) continue;
        for (dx = -r; dx <= r; dx++) {
            int xx = ix + dx;
            double k = kern[(dy + r) * kn + (dx + r)];
            if (k == 0.0 || xx < 0 || xx >= w) continue;
            dens[(size_t)yy * (size_t)w + (size_t)xx] += k * wt;
        }
    }
}

/* Returns 0 on success, 1 on allocation failure, 2 if the system diverges
 * with the given parameters.  *restarts_out receives the number of times
 * the trajectory blew up and had to be restarted. */
static int render_attractor(const RenderParams *rp, unsigned char *rgb,
                            int *restarts_out) {
    Solver sv;
    int i, restarts = 0;
    double *dens;
    double minv[3], maxv[3];
    double center[3], scale;
    double cx, cy, focal, persp;
    double c0, s0, c1, s1, c2, s2;   /* rotation cos/sin */
    double kern[81];
    double ppx = 0.0, ppy = 0.0;     /* previous projected point */
    int have_prev = 0;
    int r = rp->splat, kn = 2 * rp->splat + 1;
    size_t npix = (size_t)rp->w * (size_t)rp->h, k, nz;

    *restarts_out = 0;
    solver_init(&sv, rp);

    /* burn-in */
    for (i = 0; i < rp->burn; i++)
        if (step_checked(&sv, rp, &restarts) < 0) { *restarts_out = restarts; return 2; }

    /* Estimate the framing from a sample of the orbit (advances state too).
     * centre = bounding-box centre; scale = 1 / (distance to the farthest
     * sampled point), so the attractor fits in a unit sphere and therefore
     * stays inside the image under any yaw/pitch/roll. */
    {
        double *pts = (double *)malloc((size_t)rp->bounds_n * 3 * sizeof(double));
        size_t m = 0, q;
        double r2max = 0.0;
        if (!pts) return 1;
        minv[0] = minv[1] = minv[2] =  1e30;
        maxv[0] = maxv[1] = maxv[2] = -1e30;
        for (i = 0; i < rp->bounds_n; i++) {
            int st = step_checked(&sv, rp, &restarts);
            int a;
            if (st < 0) { free(pts); *restarts_out = restarts; return 2; }
            if (st == 0) continue;
            for (a = 0; a < 3; a++) {
                if (sv.s[a] < minv[a]) minv[a] = sv.s[a];
                if (sv.s[a] > maxv[a]) maxv[a] = sv.s[a];
                pts[m * 3 + (size_t)a] = sv.s[a];
            }
            m++;
        }
        if (m == 0) { free(pts); *restarts_out = restarts; return 2; }
        center[0] = (maxv[0] + minv[0]) * 0.5;
        center[1] = (maxv[1] + minv[1]) * 0.5;
        center[2] = (maxv[2] + minv[2]) * 0.5;
        for (q = 0; q < m; q++) {
            double ex = pts[q * 3 + 0] - center[0];
            double ey = pts[q * 3 + 1] - center[1];
            double ez = pts[q * 3 + 2] - center[2];
            double e2 = ex * ex + ey * ey + ez * ez;
            if (e2 > r2max) r2max = e2;
        }
        free(pts);
        scale = (r2max > 1e-18) ? 1.0 / sqrt(r2max) : 1.0;
    }

    /* density buffer */
    dens = (double *)calloc(npix, sizeof(double));
    if (!dens) return 1;

    /* splat kernel */
    {
        int dx, dy;
        for (dy = -r; dy <= r; dy++) {
            for (dx = -r; dx <= r; dx++) {
                int d2 = dx * dx + dy * dy;
                double wgt = 0.0;
                if (d2 <= r * r + r) {
                    wgt = (r == 0) ? 1.0 : 1.0 - sqrt((double)d2) / (double)(r + 1);
                    if (wgt < 0.0) wgt = 0.0;
                }
                kern[(dy + r) * kn + (dx + r)] = wgt;
            }
        }
    }

    /* projection constants */
    cx = (double)rp->w * 0.5 + rp->panx;   /* was: (double)rp->w * 0.5; */
    cy = (double)rp->h * 0.5 + rp->pany;   /* was: (double)rp->h * 0.5; */
    focal = (double)(rp->w < rp->h ? rp->w : rp->h) * 0.5 * FIT * rp->zoom;
    persp = 0.5;

    c0 = cos(DEG2RAD(rp->yaw));   s0 = sin(DEG2RAD(rp->yaw));
    c1 = cos(DEG2RAD(rp->pitch)); s1 = sin(DEG2RAD(rp->pitch));
    c2 = cos(DEG2RAD(rp->roll));  s2 = sin(DEG2RAD(rp->roll));

    /* render pass */
    for (i = 0; i < rp->iterations; i++) {
        double xn, yn, zn;   /* normalised coordinates  */
        double rx, ry, rz;   /* rotated                 */
        double d, px, py;
        int st, n, j;

        st = step_checked(&sv, rp, &restarts);
        if (st < 0) { free(dens); *restarts_out = restarts; return 2; }
        if (st == 0) { have_prev = 0; continue; }

        xn = (sv.s[0] - center[0]) * scale;
        yn = (sv.s[1] - center[1]) * scale;
        zn = (sv.s[2] - center[2]) * scale;

        /* rotate: yaw (z), then pitch (x), then roll (y) */
        {
            double t;
            t = c0 * xn - s0 * yn;  yn = s0 * xn + c0 * yn;  xn = t;      /* z */
            t = c1 * yn - s1 * zn;  zn = s1 * yn + c1 * zn;  yn = t;      /* x */
            t = c2 * xn + s2 * zn;  zn = -s2 * xn + c2 * zn; xn = t;      /* y */
            rx = xn; ry = yn; rz = zn;
        }

        /* perspective projection */
        d = 1.0 + persp * rz;                 /* avoids divide-by-zero   */
        if (d < 0.05) d = 0.05;
        px = cx + focal * rx / d;
        py = cy - focal * ry / d;
        if (!(px > -1e6 && px < 1e6 && py > -1e6 && py < 1e6)) { have_prev = 0; continue; }

        /* Fill the gap to the previous point with up to 64 interpolated
         * splats (weight 1/n each) so fast stretches of the orbit draw as
         * a continuous curve instead of dotted lines, while the total
         * weight per integration step stays 1 (density ~ time spent). */
        n = 1;
        if (have_prev) {
            double ddx = px - ppx, ddy = py - ppy;
            double dist = sqrt(ddx * ddx + ddy * ddy);
            if (dist > 1.0 && dist < (double)(rp->w + rp->h)) {
                n = (int)ceil(dist);
                if (n > 64) n = 64;
            }
        }
        for (j = 1; j <= n; j++) {
            double t = (double)j / (double)n;
            double qx = (have_prev && n > 1) ? ppx + (px - ppx) * t : px;
            double qy = (have_prev && n > 1) ? ppy + (py - ppy) * t : py;
            splat(dens, rp->w, rp->h, (int)floor(qx), (int)floor(qy), r, kern, 1.0 / (double)n);
        }
        ppx = px; ppy = py; have_prev = 1;
    }

    /* normalisation: pick the NORM_PCT percentile of the non-empty pixels
     * as full brightness so a few hot pixels don't crush the contrast */
    nz = 0;
    for (k = 0; k < npix; k++) if (dens[k] > 0.0) nz++;
    {
        double hi = 1.0;
        if (nz > 0) {
            double *tmp = (double *)malloc(nz * sizeof(double));
            size_t m = 0;
            if (!tmp) { free(dens); return 1; }
            for (k = 0; k < npix; k++) if (dens[k] > 0.0) tmp[m++] = dens[k];
            qsort(tmp, nz, sizeof(double), cmp_double);
            hi = tmp[(size_t)((double)(nz - 1) * NORM_PCT)];
            free(tmp);
            if (hi <= 0.0) hi = 1.0;
        }

        /* density -> colour (log scaling, cyclic palette); empty = black */
        {
            double logmax = log(1.0 + hi);
            for (k = 0; k < npix; k++) {
                if (dens[k] <= 0.0) {
                    rgb[k * 3 + 0] = rgb[k * 3 + 1] = rgb[k * 3 + 2] = 0;
                } else {
                    int idx = (int)(clamp01(log(1.0 + dens[k]) / logmax) * 255.0);
                    rgb[k * 3 + 0] = lut[0][idx];
                    rgb[k * 3 + 1] = lut[1][idx];
                    rgb[k * 3 + 2] = lut[2][idx];
                }
            }
        }
    }

    free(dens);
    *restarts_out = restarts;
    return 0;
}

/* ------------------------------------------------------------------ */
/* CLI                                                                 */
/* ------------------------------------------------------------------ */

static void usage(FILE *out) {
    fprintf(out, "chaoscope %s -- a pure-C Chaoscope-style strange attractor renderer\n",
            CHAOSCOPE_VERSION);
    fprintf(out, "Usage: chaoscope [options]\n");
    fprintf(out, "  --system <name>      attractor system (default lorenz)\n");
    fprintf(out, "  --list               list systems, parameters and palettes\n");
    fprintf(out, "  --list-tsv           same, tab-separated, for front ends\n");
    fprintf(out, "  --width <n>          image width  (default 3840)\n");
    fprintf(out, "  --height <n>         image height (default 2160)\n");
    fprintf(out, "  --iterations <n>     ODE steps to render (default 500000)\n");
    fprintf(out, "  --dt <f>             integration step (default per-system)\n");
    fprintf(out, "  --yaw <deg>          rotation around z (default 180)\n");
    fprintf(out, "  --pitch <deg>        rotation around x (default 180)\n");
    fprintf(out, "  --roll <deg>         rotation around y (default 0)\n");
    fprintf(out, "  --zoom <f>           zoom factor (default 1)\n");
    fprintf(out, "  --panx <f>           horizontal pan in pixels (default 0)\n");
    fprintf(out, "  --pany <f>           vertical pan in pixels (default 0)\n");
    fprintf(out, "  --palette <name>     rainbow|fire|ice|plasma|chaoscope (default chaoscope)\n");
    fprintf(out, "  --cycle <f>          palette cycle phase 0..1 (default 0)\n");
    fprintf(out, "  --splat <n>          splat radius 0..4 (default 2)\n");
    fprintf(out, "  --gamma <f>          output gamma 0.1..10 (default 1; >1 brightens dim regions)\n");
    fprintf(out, "  --params a,b,c,...   override leading system parameters (see --list)\n");
    fprintf(out, "  --seed <n>           RNG seed (default time-based)\n");
    fprintf(out, "  --output <file.bmp>  output file (default attractor.bmp)\n");
    fprintf(out, "  --animation <n>      render n frames (camera orbits, palette cycles)\n");
    fprintf(out, "  --prefix <str>       frame filename prefix (default 'frame')\n");
    fprintf(out, "  --version            show version and exit\n");
    fprintf(out, "  --help               show this help\n");
    fprintf(out, "Example: chaoscope --system thomas --yaw 30 --pitch 20 --palette rainbow\n");
}

static int find_system(const char *name) {
    int i;
    for (i = 0; i < NSYS; i++)
        if (strcmp(SYSTEMS[i].name, name) == 0) return i;
    return -1;
}

static int find_palette(const char *name) {
    int i;
    for (i = 0; i < NPAL; i++)
        if (strcmp(PALETTES[i].name, name) == 0) return i;
    return -1;
}

static int parse_double(const char *s, double *out) {
    char *end; double v = strtod(s, &end);
    if (end == s || *end != '\0' || !(v > -1e300 && v < 1e300)) return 0;
    *out = v; return 1;
}

static int parse_int(const char *s, int *out) {
    char *end; long v = strtol(s, &end, 10);
    if (end == s || *end != '\0' || v < -2000000000L || v > 2000000000L) return 0;
    *out = (int)v; return 1;
}

/* Parse "a,b,c" into p[0..]; at most maxn values.  Returns the count or -1. */
static int parse_params(const char *s, double *p, int maxn) {
    char buf[512];
    char *tok;
    int n = 0;
    if (strlen(s) >= sizeof(buf)) return -1;
    strcpy(buf, s);
    for (tok = strtok(buf, ","); tok; tok = strtok(NULL, ",")) {
        double v;
        if (n >= maxn || !parse_double(tok, &v)) return -1;
        p[n++] = v;
    }
    return n > 0 ? n : -1;
}

static int bad_value(const char *opt) {
    fprintf(stderr, "chaoscope: missing or invalid value for %s\n", opt);
    return 1;
}

int main(int argc, char **argv) {
    RenderParams rp;
    int i, si, pi = -1;
    int width = 3840, height = 2160;
    int iterations = 500000;
    double dt = 0.0;                     /* 0 => per-system default */
    double yaw = 180.0, pitch = 180.0, roll = 0.0;
    double zoom = 1.0, cycle = 0.0;
    double panx = 0.0, pany = 0.0;
    double gamma_v = 1.0;
    int splat_r = 2;
    char outfile[256] = "attractor.bmp";
    int anim = 0;
    char prefix[256] = "frame";
    int seed = 42;
    const char *params_arg = NULL;
    unsigned long long base_seed;
    unsigned char *rgb;
    char fname[1024];
    int status = 0, restarts = 0, total_restarts = 0;

    /* default system = lorenz */
    si = find_system("lorenz");

    /* parse args */
    for (i = 1; i < argc; i++) {
        const char *a = argv[i];
        const char *val = (i + 1 < argc) ? argv[i + 1] : NULL;

        if (strcmp(a, "--help") == 0 || strcmp(a, "-h") == 0) {
            usage(stdout);
            return 0;
        }
        else if (strcmp(a, "--version") == 0 || strcmp(a, "-V") == 0) {
            printf("chaoscope %s\n", CHAOSCOPE_VERSION);
            return 0;
        }
        else if (strcmp(a, "--list") == 0) {
            int k, q;
            printf("Systems:\n");
            for (k = 0; k < NSYS; k++) {
                printf("  %-10s %s\n", SYSTEMS[k].name, SYSTEMS[k].desc);
                printf("  %-10s params:", "");
                for (q = 0; q < SYSTEMS[k].nparams; q++)
                    printf(" %s=%g", SYSTEMS[k].pnames[q], SYSTEMS[k].def[q]);
                printf("   dt=%g\n", SYSTEMS[k].dt);
            }
            printf("Palettes:\n");
            for (k = 0; k < NPAL; k++) printf("  %s\n", PALETTES[k].name);
            return 0;
        }
        else if (strcmp(a, "--list-tsv") == 0) {
            /* machine-readable catalogue for front ends (tab separated):
             *   S <name> <desc> <dt> <nparams> {<pname> <default>}...
             *   P <palette name>                                        */
            int k, q;
            for (k = 0; k < NSYS; k++) {
                printf("S\t%s\t%s\t%.17g\t%d", SYSTEMS[k].name, SYSTEMS[k].desc,
                       SYSTEMS[k].dt, SYSTEMS[k].nparams);
                for (q = 0; q < SYSTEMS[k].nparams; q++)
                    printf("\t%s\t%.17g", SYSTEMS[k].pnames[q], SYSTEMS[k].def[q]);
                printf("\n");
            }
            for (k = 0; k < NPAL; k++) printf("P\t%s\n", PALETTES[k].name);
            return 0;
        }
        else if (strcmp(a, "--system") == 0) {
            int s = (val ? find_system(val) : -1);
            if (s < 0) { fprintf(stderr, "chaoscope: unknown system '%s' (try --list)\n", val ? val : ""); return 1; }
            si = s; i++;
        }
        else if (strcmp(a, "--palette") == 0) {
            int p = (val ? find_palette(val) : -1);
            if (p < 0) { fprintf(stderr, "chaoscope: unknown palette '%s' (try --list)\n", val ? val : ""); return 1; }
            pi = p; i++;
        }
        else if (strcmp(a, "--params") == 0)    { if (!val) return bad_value(a); params_arg = val; i++; }
        else if (strcmp(a, "--width") == 0)     { if (!val || !parse_int(val, &width))      return bad_value(a); i++; }
        else if (strcmp(a, "--height") == 0)    { if (!val || !parse_int(val, &height))     return bad_value(a); i++; }
        else if (strcmp(a, "--iterations") == 0){ if (!val || !parse_int(val, &iterations)) return bad_value(a); i++; }
        else if (strcmp(a, "--dt") == 0)        { if (!val || !parse_double(val, &dt))      return bad_value(a); i++; }
        else if (strcmp(a, "--yaw") == 0)       { if (!val || !parse_double(val, &yaw))     return bad_value(a); i++; }
        else if (strcmp(a, "--pitch") == 0)     { if (!val || !parse_double(val, &pitch))   return bad_value(a); i++; }
        else if (strcmp(a, "--roll") == 0)      { if (!val || !parse_double(val, &roll))    return bad_value(a); i++; }
        else if (strcmp(a, "--zoom") == 0)      { if (!val || !parse_double(val, &zoom))    return bad_value(a); i++; }
        else if (strcmp(a, "--panx") == 0)      { if (!val || !parse_double(val, &panx))    return bad_value(a); i++; }
        else if (strcmp(a, "--pany") == 0)      { if (!val || !parse_double(val, &pany))    return bad_value(a); i++; }
        else if (strcmp(a, "--gamma") == 0)     { if (!val || !parse_double(val, &gamma_v)) return bad_value(a); i++; }
        else if (strcmp(a, "--cycle") == 0)     { if (!val || !parse_double(val, &cycle))   return bad_value(a); i++; }
        else if (strcmp(a, "--splat") == 0)     { if (!val || !parse_int(val, &splat_r))    return bad_value(a); i++; }
        else if (strcmp(a, "--seed") == 0)      { if (!val || !parse_int(val, &seed))       return bad_value(a); i++; }
        else if (strcmp(a, "--animation") == 0) { if (!val || !parse_int(val, &anim))       return bad_value(a); i++; }
        else if (strcmp(a, "--output") == 0)    { if (!val) return bad_value(a); strncpy(outfile, val, sizeof(outfile) - 1); outfile[sizeof(outfile) - 1] = 0; i++; }
        else if (strcmp(a, "--prefix") == 0)    { if (!val) return bad_value(a); strncpy(prefix, val, sizeof(prefix) - 1); prefix[sizeof(prefix) - 1] = 0; i++; }
        else { fprintf(stderr, "chaoscope: unknown option '%s'\n", a); usage(stderr); return 1; }
    }

    /* validate */
    if (width < 1 || width > MAX_DIM || height < 1 || height > MAX_DIM) {
        fprintf(stderr, "chaoscope: width and height must be 1..%d\n", MAX_DIM);
        return 1;
    }
    if (iterations < 1) { fprintf(stderr, "chaoscope: --iterations must be >= 1\n"); return 1; }
    if (!(zoom > 0.0))  { fprintf(stderr, "chaoscope: --zoom must be > 0\n"); return 1; }
    if (!(gamma_v >= 0.1 && gamma_v <= 10.0)) { fprintf(stderr, "chaoscope: --gamma must be 0.1..10\n"); return 1; }
    lut_gamma = gamma_v;
    if (dt < 0.0)       { fprintf(stderr, "chaoscope: --dt must be > 0\n"); return 1; }
    if (anim < 0)       { fprintf(stderr, "chaoscope: --animation must be >= 0\n"); return 1; }
    if (splat_r < 0) splat_r = 0;
    if (splat_r > 4) splat_r = 4;

    /* build the render parameters (system defaults first, then overrides) */
    memset(&rp, 0, sizeof(rp));
    rp.sys = &SYSTEMS[si];
    memcpy(rp.p, rp.sys->def, sizeof(rp.p));
    if (params_arg) {
        if (parse_params(params_arg, rp.p, rp.sys->nparams) < 0) {
            fprintf(stderr, "chaoscope: --params expects 1..%d comma-separated numbers for '%s'\n",
                    rp.sys->nparams, rp.sys->name);
            return 1;
        }
    }
    rp.w = width; rp.h = height;
    rp.iterations = iterations;
    rp.burn = 800; rp.bounds_n = 100000;
    rp.dt = (dt > 0.0) ? dt : rp.sys->dt;
    rp.yaw = yaw; rp.pitch = pitch; rp.roll = roll;
    rp.zoom = zoom; rp.splat = splat_r;
    rp.panx = panx; rp.pany = pany;

    if (pi < 0) pi = find_palette("chaoscope");

    /* One seed for the whole run: every animation frame integrates the same
     * trajectory, so the framing does not jitter from frame to frame. */
    base_seed = rng_seed((unsigned long long)(unsigned int)seed);

    rgb = (unsigned char *)malloc((size_t)width * (size_t)height * 3);
    if (!rgb) { fprintf(stderr, "chaoscope: out of memory\n"); return 1; }

    if (anim > 0) {
        int f;
        for (f = 0; f < anim; f++) {
            double ph = (double)f / (double)anim;
            rp.yaw = yaw + ph * 360.0;                 /* full orbit per sequence */
            build_lut(PALETTES[pi].stops, PALETTES[pi].n, cycle + ph);
            rng_seed(base_seed);
            status = render_attractor(&rp, rgb, &restarts);
            total_restarts += restarts;
            if (status != 0) break;
            snprintf(fname, sizeof(fname), "%s_%04d.bmp", prefix, f);
            if (!write_bmp(fname, width, height, rgb)) {
                fprintf(stderr, "chaoscope: failed to write %s\n", fname);
                free(rgb); return 1;
            }
            printf("wrote %s\n", fname);
            fflush(stdout);   /* lets a front end follow progress through a pipe */
        }
    } else {
        build_lut(PALETTES[pi].stops, PALETTES[pi].n, cycle);
        status = render_attractor(&rp, rgb, &restarts);
        total_restarts += restarts;
        if (status == 0) {
            if (!write_bmp(outfile, width, height, rgb)) {
                fprintf(stderr, "chaoscope: failed to write %s\n", outfile);
                free(rgb); return 1;
            }
            printf("wrote %s (%dx%d, %d iterations, system %s)\n",
                   outfile, width, height, iterations, rp.sys->name);
        }
    }

    free(rgb);
    if (status == 1) { fprintf(stderr, "chaoscope: out of memory\n"); return 1; }
    if (status == 2) {
        fprintf(stderr, "chaoscope: system '%s' diverges with these parameters / step size; "
                        "try a smaller --dt or different --params\n", rp.sys->name);
        return 1;
    }
    if (total_restarts > 0)
        fprintf(stderr, "chaoscope: warning: trajectory diverged and was restarted %d time(s)\n",
                total_restarts);
    return 0;
}