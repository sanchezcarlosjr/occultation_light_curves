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
    linspace(m, -D/2, D/2, M);

    double **a = malloc(M * sizeof(double*));
    double **b = malloc(M * sizeof(double*));
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
// This C implementation assumes the input matrix P is a square matrix for simplicity.
// The function will allocate memory for the output matrix, perform the translation, and then return the translated matrix.
double** trasladar(double** P, int M, int smx, int smy) {
    double** MM = zeros(M, M);
    int x = M / 2;
    int y = M / 2;
    int mx = smx;
    int my = smy;

    // Y-direction translation
    if (my > 0) {
        for (int i = y + my; i < M; i++) {
            memcpy(MM[i], P[i - my], M * sizeof(double));
        }
        for (int i = 0; i < my; i++) {
            memcpy(MM[i], P[M - my + i], M * sizeof(double));
        }
        for (int i = my; i < y + my; i++) {
            memcpy(MM[i], P[i - my], M * sizeof(double));
        }
    } else if (my < 0) {
        for (int i = 0; i < y + my; i++) {
            memcpy(MM[i], P[i - my], M * sizeof(double));
        }
        for (int i = M + my; i < M; i++) {
            memcpy(MM[i], P[i - (M + my)], M * sizeof(double));
        }
        for (int i = y + my; i < M + my; i++) {
            memcpy(MM[i], P[i - my], M * sizeof(double));
        }
    } else {
        for (int i = 0; i < M; i++) {
            memcpy(MM[i], P[i], M * sizeof(double));
        }
    }

    // X-direction translation
    double** M2 = zeros(M, M);
    if (mx > 0) {
        for (int i = 0; i < M; i++) {
            memcpy(&M2[i][x + mx], &MM[i][x], (M - mx - x) * sizeof(double));
            memcpy(&M2[i][0], &MM[i][M - mx], mx * sizeof(double));
            memcpy(&M2[i][mx], &MM[i][0], (x - mx) * sizeof(double));
        }
    } else if (mx < 0) {
        for (int i = 0; i < M; i++) {
            memcpy(&M2[i][0], &MM[i][-mx], (x + mx) * sizeof(double));
            memcpy(&M2[i][M + mx], &MM[i][0], -mx * sizeof(double));
            memcpy(&M2[i][x + mx], &MM[i][x], (x - mx) * sizeof(double));
        }
    } else {
        for (int i = 0; i < M; i++) {
            memcpy(M2[i], MM[i], M * sizeof(double));
        }
    }

    freeMatrix(MM, M);
    return M2;
}


// Square Matrix size in pixels
// tamaño de matriz en metros
// oscurecimiento central en metros como si fuera circular
double** pupil_doble(int M, double D, double d) {
    double r1 = (d / 2) * 0.65;
    double r2 = sqrt(pow(d / 2, 2) - pow(r1, 2));
    double d1 = r1 * 2;
    double d2 = r2 * 2;
    double Dx = 0.45 * d1 + 0.45 * d2;
    double Dy = 0;
    // TODO: An extra argument custom sep
    int sepX = (int)((Dx / 2) / D * M);
    int sepY = (int)((Dy / 2) / D * M);

    double** P1 = zeros(M, M);
    double** P2 = zeros(M, M);
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
    double** translatedP1 = trasladar(P1, M, -sepX, sepY);
    double** translatedP2 = trasladar(P2, M, sepX, sepY);
    double** P = zeros(M, M);

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            P[i][j] = (translatedP1[i][j] + translatedP2[i][j] == 2) ? 1.0 : 0.0;  // Binarize
        }
    }

    // Free allocated memory
    freeMatrix(P1, M);
    freeMatrix(P2, M);
    freeMatrix(translatedP1, M);
    freeMatrix(translatedP2, M);

    return P;
}


// Function to generate a circular aperture
// TODO: Ellipse
double** pupilCA(int M, double D, double d) {
    double** P = zeros(M, M);
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
double** pupilSO(int M, double D, double d) {
    double** P = ones(M, M);
    int t = (int)(M * d / D);  // Size of the obstruction in pixels
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

double SNR_TAOS2(double mV) {
    // Polynomial coefficients
    double p1 = 1.5792;
    double p2 = -57.045;
    double p3 = 515.04;

    // Polynomial calculation for SNR
    double SNR = p1*mV*mV + p2*mV + p3;

    return SNR;
}

double** createSquareAperture(int M, double D, double d) {
    int t = (int)(M * d / D);  // Calculate the size of the aperture in pixels
    int c = M / 2;  // Center of the matrix
    int start = c - t / 2;
    int end = c + t / 2;

    // Allocate memory for the 2D array
    double** P = zeros(M, M);

    // Set the values within the square aperture to 1
    for (int i = start; i < end; i++) {
        for (int j = start; j < end; j++) {
            P[i][j] = 1.0;
        }
    }

    return P;
}

void fresnel () {

}

typedef struct {
    char tipo[10];  // Spectral type
    double T;       // Temperature
    double M;       // Absolute magnitude
    double L;       // Luminosity relative to the Sun
} Star;

void calc_rstar(double mV, int nEst, double ua, Star stars[], int numStars, char *tipo, double *R_star) {
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

void calculateResolutionSteps(double R_star, double plano, int M, double d, double **reso, int *resoSize) {
    double star_px = (R_star / plano) * M;
    double obj_px = (d / 4.0 / plano) * M;
    double div = ceil(star_px / obj_px);
    double rr = star_px / div;

    arange(rr, star_px+0.0001, rr);

}

// TODO: data writer