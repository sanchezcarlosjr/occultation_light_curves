//
// Created by cest on 1/31/24.
//

#include "diffraction.h"
#include "numpy.h"
#include <string.h>
#include <gsl/gsl_matrix.h>

#include <libgen.h> // For dirname function
typedef double complex** DiffractionPattern;


void getLibDir(char *libdir, const char *filepath) {
    strcpy(libdir, filepath);
    char *dir = dirname(libdir); // Extract directory name
    strcpy(libdir, dir); // Copy the directory path back into libdir
}

// Function to convert Cartesian coordinates to Polar coordinates
void cart2pol(double x, double y, double *phi, double *rho) {
    *rho = sqrt(x * x + y * y);  // Calculate the radius
    *phi = atan2(y, x);          // Calculate the angle in radians
}

void pol2cart(double rho, double phi, double *x, double *y) {
    *x = rho * cos(phi);  // Calculate the x coordinate
    *y = rho * sin(phi);  // Calculate the y coordinate
}

void pupilCO(int M, double D, double d, double **P) {
    double *m = malloc(M * sizeof(double));
    linspace(m, -D / 2, D / 2, M);

    double **a = malloc(M * sizeof(double *));
    double **b = malloc(M * sizeof(double *));
    for (int i = 0; i < M; i++) {
        a[i] = malloc(M * sizeof(double));
        b[i] = malloc(M * sizeof(double));
    }

    meshgrid(m, m, M, a, b);

    double phi, rho;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            cart2pol(a[i][j], b[i][j], &phi, &rho);
            P[i][j] = (rho >= d / 2) ? 1.0 : 0.0;
        }
    }

    for (int i = 0; i < M; i++) {
        free(a[i]);
        free(b[i]);
    }
    free(a);
    free(b);
    free(m);
}


// The function will allocate memory for the output matrix, perform the translation, and then return the translated matrix.
double complex **translate(DiffractionPattern P, int dx, int dy, int smx, int smy) {
    double complex **MM = allocateComplex2DArray(smy, smx); // Allocate memory for the translated matrix
    if (!MM) return NULL;
    zerosComplex(MM, smy, smx); // Initialize MM with zeros

    // Translate in Y-direction
    for (int i = 0; i < smy; i++) {
        int sourceRow = (i - dy + smy) % smy; // Calculate source row considering circular shift
        for (int j = 0; j < smx; j++) {
            MM[i][j] = P[sourceRow][j];
        }
    }

    // Translate in X-direction and create a new matrix for the result
    double complex **M2 = allocateComplex2DArray(smy, smx); // Allocate memory for the final translated matrix
    if (!M2) {
        freeComplex2DArray(MM, smy); // Clean up MM before returning
        return NULL;
    }
    zerosComplex(M2, smy, smx); // Initialize M2 with zeros

    for (int i = 0; i < smy; i++) {
        for (int j = 0; j < smx; j++) {
            int sourceCol = (j - dx + smx) % smx; // Calculate source column considering circular shift
            M2[i][j] = MM[i][sourceCol];
        }
    }

    freeComplex2DArray(MM, smy); // Free the temporary matrix
    return M2; // Return the final translated matrix
}

// Square Matrix size in pixels
// Matrix size in meters
// Central dimming in meters as if it were circular
ComplexMatrix dynamicPupilDoble(int M, double D, double d, int sepX, int sepY, double r1, double r2) {
    ComplexMatrix P1 = allocateComplex2DArray(M, M);
    ComplexMatrix P2 = allocateComplex2DArray(M, M);
    zerosComplex(P1, M, M);
    zerosComplex(P2, M, M);

    double phi, rho;

    // Generate obstructions
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            double x = -D / 2 + D * i / M;
            double y = -D / 2 + D * j / M;
            cart2pol(x, y, &phi, &rho);
            P1[i][j] = (rho >= r1) ? 1.0 : 0.0;  // Large obstruction
            P2[i][j] = (rho >= r2) ? 1.0 : 0.0;  // Small obstruction
        }
    }

    // Translate and combine obstructions
    ComplexMatrix translatedP1 = translate(P1, M, M, -sepX, sepY);
    ComplexMatrix translatedP2 = translate(P2, M, M, sepX, sepY);
    ComplexMatrix P = allocateComplex2DArray(M, M);
    zerosComplex(P, M, M);

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            P[i][j] = (translatedP1[i][j] + translatedP2[i][j] == 2) ? 1.0 : 0.0;  // Binarize
        }
    }

    // Free allocated memory
    freeComplex2DArray(P1, M);
    freeComplex2DArray(P2, M);
    freeComplex2DArray(translatedP1, M);
    freeComplex2DArray(translatedP2, M);

    return P;
}

ComplexMatrix pupilDobleWithCustomSeparator(int M, double D, double d, int sepX, int sepY) {
    double r1 = (d / 2) * 0.65;
    double r2 = sqrt(pow(d / 2, 2) - pow(r1, 2));
    double d1 = r1 * 2;
    double d2 = r2 * 2;
    double Dx = 0.45 * d1 + 0.45 * d2;
    double Dy = 0;
    return dynamicPupilDoble(M,D, d, (int) ((Dx / 2) / D * M) + sepX, (int) ((Dy / 2) / D * M) + sepY, r1, r2);
}

ComplexMatrix pupilDoble(int M, double D, double d) {
    return pupilDobleWithCustomSeparator(M, D, d, 0, 0);
}


// Function to generate a circular aperture
// TODO: Ellipse
ComplexMatrix pupilCA(int M, double D, double d) {
    ComplexMatrix P = allocateComplex2DArray(M, M);
    zerosComplex(P, M, M);
    double step = D / M;
    double phi, rho;

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            double x = -D / 2 + i * step;  // Calculate x coordinate
            double y = -D / 2 + j * step;  // Calculate y coordinate
            cart2pol(x, y, &phi, &rho);    // Convert to polar coordinates
            P[i][j] = (rho <= d / 2) ? 1.0 : 0.0;  // Mark points within the aperture
        }
    }

    return P;
}

// Function to generate a square obstruction
double **pupilSO(int M, double D, double d) {
    double **P = ones(M, M);
    int t = (int) (M * d / D);  // Size of the obstruction in pixels
    int c = M / 2;  // Center of the matrix

    // Calculate the bounds of the obstruction
    int start = c - t / 2;
    int end = c + t / 2;

    // Create the square obstruction
    for (int i = start; i < end; i++) {
        for (int j = start; j < end; j++) {
            if (i >= 0 && i < M && j >= 0 && j < M) {
                P[i][j] = 0.0;  // Mark the obstruction area with zeros
            }
        }
    }

    return P;
}


void pupilSA(ComplexMatrix P, int M, double D, double d) {
    int t = M * d / D; // Calculate the size of the central obscuration in pixels
    int c = M / 2;     // Center of the array

    // Initialize the entire array to 0
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            P[i][j] = 0.0;
        }
    }

    // Set the central square area to 1
    for (int i = c - t / 2; i < c + t / 2; i++) {
        for (int j = c - t / 2; j < c + t / 2; j++) {
            if (i >= 0 && i < M && j >= 0 && j < M) { // Check boundaries
                P[i][j] = 1.0;
            }
        }
    }
}

void fresnel(double complex *U0, int nx, int ny, int M, double plano, double z, double lmda,
             double complex *output_intensity) {
    double k = 2 * M_PI / lmda;
    double x = (plano / M) * nx; // Normally nx = M, so x = plano in meters
    double y = (plano / M) * ny;
    double fx = 1 / x; // Spatial frequency in m**-1
    double fy = 1 / y;

    // Allocate memory for u, v, O, H, and U arrays
    double complex *u = (double complex *) fftw_malloc(sizeof(double complex) * nx);
    double complex *v = (double complex *) fftw_malloc(sizeof(double complex) * ny);
    double complex *O = (double complex *) fftw_malloc(sizeof(double complex) * nx * ny);
    double complex *H = (double complex *) fftw_malloc(sizeof(double complex) * nx * ny);
    double complex *U = (double complex *) fftw_malloc(sizeof(double complex) * nx * ny);

    // Compute u and v arrays
    for (int i = 0; i < nx; i++) {
        u[i] = (i - nx / 2) * fx;
    }
    for (int j = 0; j < ny; j++) {
        v[j] = (j - ny / 2) * fy;
    }

    // Apply FFT to U0 using FFTW
    fftw_plan plan_forward = fftw_plan_dft_2d(nx, ny, U0, O, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

    // Compute H array
    // np.exp(1j*k*z)*np.exp(-1j*np.pi*(lmda*z)*(u**2+v**2))
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double u_squared = creal(u[i]) * creal(u[i]);
            double v_squared = creal(v[j]) * creal(v[j]);
            H[i * ny + j] = cexp(I * k * z) * cexp(-I * M_PI * lmda * z * (u_squared + v_squared));
        }
    }

    // Multiply O and H element-wise
    for (int i = 0; i < nx * ny; i++) {
        U[i] = O[i] * H[i];
    }

    // Apply inverse FFT to U using FFTW
    fftw_plan plan_backward = fftw_plan_dft_2d(nx, ny, U, output_intensity, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_backward);

    // Compute the intensity pattern
    for (int i = 0; i < nx * ny; i++) {
        output_intensity[i] = cpow(cabs(output_intensity[i]), 2);
    }

    // Free memory and destroy FFTW plans
    fftw_free(u);
    fftw_free(v);
    fftw_free(O);
    fftw_free(H);
    fftw_free(U);
    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);
}

/* listadat.txt--> A0=1;A1=2;A2=3;A3=4;A4=5;A5=6;A7=7;F0=8;F2=9;F3=10;F5=11;F6=12;F7=13;F8=14;
 * G0=15;G1=16;G2=17;G5=18;G8=19;K0=20;K1=21;K2=22;K3=23;K4=24;K5=25;K7=26;
 * M0=27;M1=28;M2=29;M3=30;M4=31;M5=32;M6=33;M7=34;M8=35
*/
void
spectra(double complex *U0, int nx, int ny, int M, double plano, double z, int nEst, int nLmdas, double complex *acc) {
    char libdir[1024];
    getLibDir(libdir, __FILE__);
    // Define file paths - you'll need to define how libdir is set
    char refFilePath[1024];
    sprintf(refFilePath, "%s/listadat.txt", libdir);

    FILE *refFile = fopen(refFilePath, "r");
    if (!refFile) {
        perror("Failed to open reference file");
        return;
    }

    // Read the nEst-th line from the reference file to get the data file name
    char line[256];
    char dataFileName[256];
    int currentLine = 0;
    while (fgets(line, sizeof(line), refFile)) {
        currentLine++;
        if (currentLine == nEst) {
            strncpy(dataFileName, line, strlen(line) - 1);  // Remove newline character
            break;
        }
    }
    fclose(refFile);

    // Construct the path to the data file
    char dataFilePath[1024];
    sprintf(dataFilePath, "%s/spectra/%s", libdir, dataFileName);

    // Open the data file
    FILE *dataFile = fopen(dataFilePath, "r");
    if (!dataFile) {
        perror("Failed to open data file");
        return;
    }

    int U0Size = nx * ny;

    // Read data from the file and perform calculations
    double lamda, peso;
    double complex *output_intensity = (double complex *) malloc(U0Size * sizeof(double complex));
    int count = 0;

    // Assuming the data file has two columns: lamda and peso
    while (fscanf(dataFile, "%lf,%lf", &lamda, &peso) == 2) {
        fresnel(U0, nx, ny, M, plano, z, lamda * 1e-10, output_intensity);

        // Weight the result by 'peso' and accumulate
        for (int i = 0; i < U0Size; i++) {
            acc[i] += output_intensity[i] * peso;
        }

        count++;
        if (count >= nLmdas) break;  // Stop after nLmdas iterations
    }

    fclose(dataFile);
    free(output_intensity);

    // Normalize and finalize the result
    for (int i = 0; i < U0Size; i++) {
        acc[i] /= acc[0];  // Normalize by the first element
    }


}

typedef struct {
    char tipo[10];  // Spectral type
    double T;       // Temperature
    double M;       // Absolute magnitude
    double L;       // Luminosity relative to the Sun
} Star;

/*
 * Function to calculate the apparent radii of stars
% mV --> Apparent magnitude
% nEst --> star number
% ua --> Distance to the object in astronomical units
% Absolute magnitudes in order from A0 type stars to M8 type
* Star[] -> The caller must read the data, for example, from "estrellas.csv", and convert it into a Star struct.
% M0 = [1.5 1.7 1.8 2.0 2.1 2.2 2.4 3.0 3.3 3.5 3.7 4.0 4.3 4.4 4.7 4.9 5.0...
%     5.2 2.6 6.0 6.2 6.4 6.7 7.1 7.4 8.1 8.7 9.4 10.1 10.7 11.2 12.3 13.4...
%     13.9 14.4];
OUT --> type, R_star: spectral type chosen and calculated star radius, respectively
 */
void calcRstar(double mV, int nEst, double ua, Star stars[], int numStars, char *tipo, double *R_star) {
    // Constants
    double ua_meters = 1.496e11 * ua;  // Distance in meters
    double Tsol = 5780;                // Solar temperature in Kelvin
    double Rsol = 6.96e8;              // Solar radius in meters

    if (nEst < 1 || nEst > numStars) {
        printf("Invalid star number.\n");
        return;
    }

    // Select the star
    Star star = stars[nEst - 1];

    // Calculations
    double d1 = pow(10, (mV - star.M + 5) / 5);
    double d = 3.085e16 * d1;  // Convert from parsecs to meters
    double Rst = sqrt(star.L) / pow(star.T / Tsol, 2);  // Star radius in Rsol
    double alfa = (Rsol * Rst) / d;  // Angular size of the star in radians
    *R_star = alfa * ua_meters;  // Size of the star in meters

    // Copy spectral type to output
    snprintf(tipo, 10, "%s", star.tipo);
}

void promedioPD(double complex *diffractionPattern, double R_star, double plano, int M, double d,
                double complex *intensityOut) {
    double star_px = (R_star / plano) * M;
    double obj_px = (d / 4.0 / plano) * M;
    int div = (int) ceil(star_px / obj_px);
    double rr = star_px / div;

    // Assuming a maximum number of steps for allocation
    int maxSteps = (int) ceil((star_px - rr) / rr) + 1;

    int co = 1;
    for (int k1 = 0; k1 < maxSteps; k1++) {
        double currentReso = rr + k1 * rr;
        if (currentReso > star_px) break;

        double perim = 2 * M_PI * currentReso;
        int paso = (int) ceil(perim / obj_px);
        double resot = 2 * M_PI / paso;

        for (double teta = resot; teta < 2 * M_PI + .0001; teta += resot) {
            double dx = currentReso * cos(teta);
            double dy = currentReso * sin(teta);

            double complex **tempIntensity = translate(&diffractionPattern, M, M, dx, dy);

            // Accumulate translated intensity into intensityOut
            for (int i = 0; i < M * M; i++) {
                intensityOut[i] += *tempIntensity[i];
            }

            co++;
            free(tempIntensity);
        }
    }

    // Combine the original and translated intensity images and normalize
    for (int i = 0; i < M * M; i++) {
        intensityOut[i] = (intensityOut[i] + diffractionPattern[i]) / co;
    }


}

void linspace2(double start, double end, int num, double *result) {
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; i++) {
        result[i] = start + i * step;
    }
}

void extraerPerfil(double **I0, int M, double D, double T, double b, double *x, double *y) {
    double m2p = M / D;
    T = T * M_PI / 180; // Convert angle to radians

    linspace2(-D / 2, D / 2, M, x);

    for (int k = 0; k < M; k++) {
        double x1 = x[k] * cos(T) - b * sin(T);
        double x2 = x[k] * sin(T) + b * cos(T);

        // Convert to pixel numbers
        int hp = (int) (m2p * x1 + M / 2);
        int vp = (int) (m2p * x2 + M / 2);

        // Ensure hp and vp are within the bounds of the array
        if (hp >= 0 && hp < M && vp >= 0 && vp < M) {
            y[k] = I0[vp][hp]; // Access intensity value at (vp, hp)
        } else {
            y[k] = 0; // Out of bounds, set intensity to 0
        }
    }
}

double calcPlano(double d, double lmda, double ua) {
    const double AU_TO_METERS = 1.496e11; // Conversion factor from AU to meters
    double z = ua * AU_TO_METERS; // Distance in meters
    double fscale = sqrt(lmda * z / 2); // Fresnel scale
    double Rho = d / (2 * fscale);
    double plano = (50 * d) / Rho; // Size of the plane in meters
    return plano;
}

double SNR_TAOS2(double mV) {
    // Polynomial coefficients
    double p1 = 1.5792;
    double p2 = -57.045;
    double p3 = 515.04;

    // Polynomial calculation for SNR
    double SNR = p1 * mV * mV + p2 * mV + p3;

    return SNR;
}


/*
 * Add Poisson noise to an image
diffractionPattern --> image matrix
mV --> apparent magnitude of the star
OUT --> In: matrix with added noise, assuming NOISE=1/SNR calculated from TAOS-II
 */
void addNoise(DiffractionPattern diffractionPattern, int M, double mV) {
    double noise = 1 / SNR_TAOS2(mV);

    // Initialize GSL random number generator
    const gsl_rng_type *T;
    gsl_rng *GSLRandomNumberGenerator;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    GSLRandomNumberGenerator = gsl_rng_alloc(T);

    double mean = 0.0; // For calculating the mean of the noise mask

    // Add Poisson noise to each pixel and calculate mean of the noise mask
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            double poisson_noise = gsl_ran_poisson(GSLRandomNumberGenerator, diffractionPattern[i][j]);
            diffractionPattern[i][j] += poisson_noise; // Add Poisson noise to the original image
            mean += poisson_noise;
        }
    }
    mean /= (M * M); // Calculate mean of the noise mask

    // Normalize the noise mask and add weighted noise according to TAOS-II
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            double normalized_noise = (diffractionPattern[i][j] / mean) * noise - noise;
            diffractionPattern[i][j] += normalized_noise; // Add the normalized and weighted noise to the image
        }
    }

    gsl_rng_free(GSLRandomNumberGenerator); // Free the GSL random number generator
}

/*
 * Function to sample the diffraction profile by obtaining the average point
lc --> diffraction profile or light curve
D --> size of the plane in meters
vr --> object velocity ~5000 m/s (positive if it goes against Earth's velocity)
fps --> frames per second of the camera, 20 for TAOS-2
toff --> Time offset within the sampling period
vE --> translational velocity of the Earth == 29800 m/s
opangle --> angle from the object's opposition: O, S, E
ua --> Distance of the object in Astronomical Units
OUT --> s_lin, lc_lin, s_pun, lc_pun: time vectors for lines, sample in lines, time at points, and sample at points, RESPECTIVELY
 */
void sampling(double *lc, int tam, double D, double vr, double fps, double toff, double vE, double opangle, double ua,
              double **s_lin, double **lc_lin, double **s_pun, double **lc_pun, int *size_s_lin, int *size_lc_lin, int *size_s_pun, int *size_lc_pun) {
    double T = 1 / fps;  // Exposure time
    double OA = opangle * M_PI / 180;  // Opposition angle in radians

    // Tangential velocity of the object relative to Earth
    double Vt = vE * (cos(OA) - sqrt((1 / ua) * (1 - (1 / (ua * ua)) * sin(OA) * sin(OA)))) + vr;

    double t = D / Vt;  // Visibility of the plane in seconds
    int Nm = (int)(t / T);  // Total number of samples in the observation plane
    int dpix = tam / Nm;
    int pixoffset = (int)(toff * fps);  // Convert time offset to pixel offset

    // Allocate memory for output arrays
    *size_s_lin = tam;  // Adjust the size according to your sampling
    *size_lc_lin = tam;
    *size_s_pun = Nm;  // Adjust the size according to your sampling
    *size_lc_pun = Nm;
    *s_lin = (double *)malloc(*size_s_lin * sizeof(double));
    *lc_lin = (double *)malloc(*size_lc_lin * sizeof(double));
    *s_pun = (double *)malloc(*size_s_pun * sizeof(double));
    *lc_pun = (double *)malloc(*size_lc_pun * sizeof(double));

    // Sample the diffraction profile
    // Assuming the sampling logic is implemented similar to the Python version
    // You need to translate the logic of sampling `lc` into `lc_lin` and `lc_pun` here

    // Fill the time vectors `s_lin` and `s_pun`
    for (int i = 0; i < *size_s_lin; i++) {
        (*s_lin)[i] = -t / 2 + i * (t / (*size_s_lin - 1));
    }
    for (int i = 0; i < *size_s_pun; i++) {
        (*s_pun)[i] = -t / 2 + i * (t / (*size_s_pun - 1));
    }
}


// Function to find peaks in the light curve using the derivative method
// x, y: Arrays containing the occultation data (distance and amplitude)
// n: Number of elements in the x and y arrays
// D: Diameter of the object in meters
// fil: Threshold value for identifying peaks, default is 0.005
// peaks: Array to store the indices of the peaks
// peakValues: Array to store the values of the peaks
// Returns the number of peaks found
int searchPeaks(double *x, double *y, int n, double D, double fil, int **peaks, double **peakValues) {
    double *yp = (double *)malloc((n - 1) * sizeof(double));  // Array for the derivative of y
    int *tempPeaks = (int *)malloc(n * sizeof(int));          // Temporary array for peak indices
    double *tempPeakValues = (double *)malloc(n * sizeof(double));  // Temporary array for peak values
    int count = 0;  // Counter for the number of peaks

    // Calculate the derivative of y
    for (int i = 0; i < n - 1; i++) {
        yp[i] = y[i + 1] - y[i];
    }

    // Find indices where the derivative is close to 0 and within the region of interest
    for (int i = 0; i < n - 1; i++) {
        if (fabs(yp[i]) < fil && fabs(x[i]) < (D / 2)) {
            tempPeaks[count] = i;
            tempPeakValues[count] = y[i];
            count++;
        }
    }

    // Allocate memory for the output arrays with the exact number of peaks
    *peaks = (int *)malloc(count * sizeof(int));
    *peakValues = (double *)malloc(count * sizeof(double));

    // Copy data to the output arrays, filtering out duplicates
    int lastIdx = -1;  // Last index added to the peaks array
    for (int i = 0; i < count; i++) {
        if (i == 0 || fabs(tempPeakValues[i] - tempPeakValues[lastIdx]) > fil) {
            (*peaks)[lastIdx + 1] = tempPeaks[i];
            (*peakValues)[lastIdx + 1] = tempPeakValues[i];
            lastIdx++;
        }
    }

    // Free temporary arrays
    free(yp);
    free(tempPeaks);
    free(tempPeakValues);

    return lastIdx + 1;  // Return the number of peaks found
}