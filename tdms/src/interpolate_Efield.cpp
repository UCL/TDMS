#include "interpolate_Efield.h"
#include "interpolation_methods.h"

/**
 * @brief Interpolate the Ex field component to the centre of a Yee cell
 * 
 * @param[in] Exy,Exz split components of the Yee cell 
 * @param[in] i,j,k Yee cell index 
 * @param[in] I Number of Yee cells in the x-dimension
 * @param[out] Ex Interpolated value of the Ex field at centre of Yee cell i,j,k
 */
void interpolateTimeDomainEx(double ***Exy, double ***Exz, int i, int j, int k, int I, double *Ex) {
    // determine the interpolation scheme to use
    interp_scheme scheme = determineInterpScheme(I, i);
}

/**
 * @brief Interpolate a component of the E-field to the centre of a Yee cell
 *
 * We enforce {a,b,c} = {i,j,k}, with the conventions a=i, b=j, c=k; a=j, b=i, c=k; a=k, b=i, c=j.
 * Indices and variables are named assuming the a=i case, EG (a,b,c) will specify the index of the Yee cell.
 *
 * @param[in] dimension {x, y, z} Identifies which direction "a"-direction represents
 * @param[in] Eab, Eac Split components of the Yee cell
 * @param[in] a,b,c Index of the Yee cell
 * @param[in] A Total number of Yee cells in the a-direction
 * @param[out] Ea The interpolated field component at the centre of the Yee cell
 */
/*
void interpolateTimeDomainEcomponent(dimension aID, double ***Eab, double ***Eac, int a, int b, int c, int A, double *Ea)
{
    // which interpolation scheme can be applied?
    interp_scheme interp_type = determineInterpScheme(A, a);

    switch (interp_type)
    {
    // perform BAND_LIMITED interpolation
    case BAND_LIMITED:
    {
        // Array for performing bandwidth limited interpolation obtained using Matlab's interp function
        const int Nbvec = 8;
        const double bvec[Nbvec] = {-0.006777513830606, 0.039457774230186, -0.142658093428622, 0.609836360661632, 0.609836360661632, -0.142658093428622, 0.039457774230186, -0.006777513830606};
        *Ea = 0.;
        switch (aID)
        {
        case x:
            // a is the i-direction, so b=j and c=k
            for (int ind = 0; ind < Nbvec; ind++)
            {
                *Ea += (Eab[c][b][a - Nbvec / 2 + ind] + Eac[c][b][a - Nbvec / 2 + ind]) * bvec[ind];
            }
            break;
        case y:
            // a is the j-direction, so b=i and c=k
            for (int ind = 0; ind < Nbvec; ind++)
            {
                *Ea += (Eab[c][a - Nbvec / 2 + ind][b] + Eac[c][a - Nbvec / 2 + ind][b]) * bvec[ind];
            }
        case z:
            // a is the k-direction, so b=i and c=j
            for (int ind = 0; ind < Nbvec; ind++)
            {
                *Ea += (Eab[a - Nbvec / 2 + ind][c][b] + Eac[a - Nbvec / 2 + ind][c][b]) * bvec[ind];
            }
        default:
            throw out_of_range("Invalid interpolation dimension (" + to_string(aID) + "), expected 0(i), 1(j), or 2(k).\n");
            break;
        }
        break;
    }
    // perform cubic interpolation (centre centre)
    case INTERP1:
    {
        switch (aID)
        {
        case x:
            // a is the i-direction, so b=j and c=k
            *Ea = interp1(Eab[c][b][a - 2] + Eac[c][b][a - 2], Eab[c][b][a - 1] + Eac[c][b][a - 1], Eab[c][b][a] + Eac[c][b][a], Eab[c][b][a + 1] + Eac[c][b][a + 1]);
            break;
        case y:
            // a is the j-direction, so b = i, and c = k
            *Ea = interp1(Eab[c][a - 2][b] + Eac[c][a - 2][b], Eab[c][a - 1][b] + Eac[c][a - 1][b], Eab[c][a][b] + Eac[c][a][b], Eab[c][a + 1][b] + Eac[c][a + 1][b]);
            break;
        case z:
            // a is the k-direction, so b = i, and c = j
            *Ea = interp1(Eab[a - 2][c][b] + Eac[a - 2][c][b], Eab[a - 1][c][b] + Eac[a - 1][c][b], Eab[a][c][b] + Eac[a][c][b], Eab[a + 1][c][b] + Eac[a + 1][c][b]);
            break;
        default:
            throw out_of_range("Invalid interpolation dimension (" + to_string(aID) + "), expected 0(i), 1(j), or 2(k).\n");
            break;
        }
    }
    // perform cubic interpolation (left centre)
    case INTERP2:
    {
        switch (aID)
        {
        case x:
            // a is the i-direction, so b=j and c=k
            *Ea = interp2(Eab[c][b][a - 1] + Eac[c][b][a - 1], Eab[c][b][a] + Eac[c][b][a], Eab[c][b][a + 1] + Eac[c][b][a + 1], Eab[c][b][a + 2] + Eac[c][b][a + 2]);
            break;
        case y:
            // a is the j-direction, so b = i, and c = k
            *Ea = interp2(Eab[c][a - 1][b] + Eac[c][a - 1][b], Eab[c][a][b] + Eac[c][a][b], Eab[c][a + 1][b] + Eac[c][a + 1][b], Eab[c][a + 2][b] + Eac[c][a + 2][b]);
            break;
        case z:
            // a is the k-direction, so b = i, and c = j
            *Ea = interp2(Eab[a - 1][c][b] + Eac[a - 1][c][b], Eab[a][c][b] + Eac[a][c][b], Eab[a + 1][c][b] + Eac[a + 1][c][b], Eab[a + 2][c][b] + Eac[a + 2][c][b]);
            break;
        default:
            throw out_of_range("Invalid interpolation dimension (" + to_string(aID) + "), expected 0(i), 1(j), or 2(k).\n");
            break;
        }
    }
    // perform cubic interpolation (right centre)
    case INTERP3:
    {
        switch (aID)
        {
        case x:
            // a is the i-direction, so b=j and c=k
            *Ea = interp3(Eab[c][b][a - 3] + Eac[c][b][a - 3], Eab[c][b][a - 2] + Eac[c][b][a - 2], Eab[c][b][a - 1] + Eac[c][b][a - 1], Eab[c][b][a] + Eac[c][b][a]);
            break;
        case y:
            // a is the j-direction, so b = i, and c = k
            *Ea = interp3(Eab[c][a - 3][b] + Eac[c][a - 3][b], Eab[c][a - 2][b] + Eac[c][a - 2][b], Eab[c][a - 1][b] + Eac[c][a - 1][b], Eab[c][a][b] + Eac[c][a][b]);
            break;
        case z:
            // a is the k-direction, so b = i, and c = j
            *Ea = interp3(Eab[a - 3][c][b] + Eac[a - 3][c][b], Eab[a - 2][c][b] + Eac[a - 2][c][b], Eab[a - 1][c][b] + Eac[a - 1][c][b], Eab[a][c][b] + Eac[a][c][b]);
            break;
        default:
            throw out_of_range("Invalid interpolation dimension (" + to_string(aID) + "), expected 0(i), 1(j), or 2(k).\n");
            break;
        }
    }
    }
}
*/

    /**
     * @brief Interpolate the E-field to the origin of the Yee cell in the time domain
     *
     * This function calls the interpolation methods for each component of the E-field separately.
     *
     * @param[in] Exy,Exz,Eyx,Eyz,Ezx,Ezy Split components of the Yee cell
     * @param i,j,k Index of the Yee cell to interpolate to the centre of
     * @param I,J,K Total number of Yee cells in the i,j,k directions respectively
     * @param Ex,Ey,Ez Interpolated E-field x,y,z components (respectively)
     */
    /*
void interpolateTimeDomainE(double ***Exy, double ***Exz, double ***Eyx, double ***Eyz, double ***Ezx, double ***Ezy, int i, int j, int k, int I, int J, int K, double *Ex, double *Ey, double *Ez)
{
    // interpolate each of the field components in sequence
    interpolateTimeDomainEcomponent(x, Exy, Exz, i, j, k, I, Ex);
    interpolateTimeDomainEcomponent(y, Eyx, Eyz, j, i, k, J, Ey);
    interpolateTimeDomainEcomponent(z, Ezx, Ezy, k, i, j, K, Ez);
}
*/