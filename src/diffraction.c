//
// Created by cest on 1/31/24.
//

#include "diffraction.h"
#include "numpy.h"
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <fftw3.h>

#include <libgen.h> // For dirname function

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

double **allocate2DArray(int rows, int cols) {
    double **array = (double **)malloc(rows * sizeof(double *));
    if (!array) return NULL;

    for (int i = 0; i < rows; i++) {
        array[i] = (double *)malloc(cols * sizeof(double));
        if (!array[i]) {
            // Free previously allocated memory in case of failure
            for (int j = 0; j < i; j++) free(array[j]);
            free(array);
            return NULL;
        }
    }
    return array;
}

// Function to free a dynamically allocated 2D array
void free2DArray(double **array, int rows) {
    for (int i = 0; i < rows; i++) free(array[i]);
    free(array);
}

// Function to initialize a 2D array with zeros
void zeros2(double **array, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            array[i][j] = 0.0;
        }
    }
}

double complex **allocateComplex2DArray(int rows, int cols) {
    double complex **array = (double complex **)malloc(rows * sizeof(double complex *));
    if (!array) return NULL;

    for (int i = 0; i < rows; i++) {
        array[i] = (double complex *)malloc(cols * sizeof(double complex));
        if (!array[i]) {
            // Free previously allocated memory in case of failure
            for (int j = 0; j < i; j++) free(array[j]);
            free(array);
            return NULL;
        }
    }
    return array;
}

// Function to initialize a 2D array of complex numbers with zeros
void zerosComplex(double complex **array, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            array[i][j] = 0.0 + 0.0 * I; // Initialize to complex zero
        }
    }
}

// Function to free a dynamically allocated 2D array of complex numbers
void freeComplex2DArray(double complex **array, int rows) {
    for (int i = 0; i < rows; i++) free(array[i]);
    free(array);
}

// The function will allocate memory for the output matrix, perform the translation, and then return the translated matrix.
double complex **translate(double complex **P, int dx, int dy, int smx, int smy) {
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
// tamaño de matriz en metros
// oscurecimiento central en metros como si fuera circular
double **pupil_doble(int M, double D, double d) {
    double r1 = (d / 2) * 0.65;
    double r2 = sqrt(pow(d / 2, 2) - pow(r1, 2));
    double d1 = r1 * 2;
    double d2 = r2 * 2;
    double Dx = 0.45 * d1 + 0.45 * d2;
    double Dy = 0;
    // TODO: An extra argument custom sep
    int sepX = (int) ((Dx / 2) / D * M);
    int sepY = (int) ((Dy / 2) / D * M);

    double complex **P1 = allocateComplex2DArray(M, M);
    double complex **P2 = allocateComplex2DArray(M, M);
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
    double complex ** translatedP1 = translate(P1, M, M, -sepX, sepY);
    double complex ** translatedP2 = translate(P2, M, M, sepX, sepY);
    double **P = zeros(M, M);

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


// Function to generate a circular aperture
// TODO: Ellipse
double **pupilCA(int M, double D, double d) {
    double **P = zeros(M, M);
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


void pupilSA(Matrix P, int M, double D, double d) {
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

void fresnel(double complex *U0, int nx, int ny, int M, double plano, double z, double lmda, double complex *output_intensity) {
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
void spectra(double complex *U0, int nx, int ny, int M, double plano, double z, int nEst, int nLmdas, double complex * acc) {
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

    int U0Size = nx*ny;

    // Read data from the file and perform calculations
    double lamda, peso;
    double complex *output_intensity = (double complex *)malloc(U0Size * sizeof(double complex));
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

void promedioPD(double complex *diffractionPattern, double R_star, double plano, int M, double d, double complex *intensityOut) {
    double star_px = (R_star / plano) * M;
    double obj_px = (d / 4.0 / plano) * M;
    int div = (int)ceil(star_px / obj_px);
    double rr = star_px / div;

    // Assuming a maximum number of steps for allocation
    int maxSteps = (int)ceil((star_px - rr) / rr) + 1;

    int co = 1;
    for (int k1 = 0; k1 < maxSteps; k1++) {
        double currentReso = rr + k1 * rr;
        if (currentReso > star_px) break;

        double perim = 2 * M_PI * currentReso;
        int paso = (int)ceil(perim / obj_px);
        double resot = 2 * M_PI / paso;

        for (double teta = resot; teta < 2 * M_PI + .0001; teta += resot) {
            double dx = currentReso * cos(teta);
            double dy = currentReso * sin(teta);

            double complex** tempIntensity = translate(&diffractionPattern, M,  M, dx, dy);

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


double SNR_TAOS2(double mV) {
    // Polynomial coefficients
    double p1 = 1.5792;
    double p2 = -57.045;
    double p3 = 515.04;

    // Polynomial calculation for SNR
    double SNR = p1 * mV * mV + p2 * mV + p3;

    return SNR;
}

double **createSquareAperture(int M, double D, double d) {
    int t = (int) (M * d / D);  // Calculate the size of the aperture in pixels
    int c = M / 2;  // Center of the matrix
    int start = c - t / 2;
    int end = c + t / 2;

    // Allocate memory for the 2D array
    double **P = zeros(M, M);

    // Set the values within the square aperture to 1
    for (int i = start; i < end; i++) {
        for (int j = start; j < end; j++) {
            P[i][j] = 1.0;
        }
    }

    return P;
}



