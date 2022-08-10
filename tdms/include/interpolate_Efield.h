/* Interpolation functions for each of the E-field components.
Decided to do with 3 separate functions for readability; as this avoids switches with a lot of cases within functions, and the interpolation scheme might be different in each direction anyway, so there is no point having a centralised function for all of them.

If we want to interpolate all the E-field components, there are functions for this which call the individual component functions separately.
*/

void interpolateTimeDomainEx(double ***Exy, double ***Exz, int i, int j, int k, int I, double *Ex);
void interpolateTimeDomainEy(double ***Eyx, double ***Eyz, int i, int j, int k, int J, double *Ey);
void interpolateTimeDomainEz(double ***Ezx, double ***Ezy, int i, int j, int k, int K, double *Ez);
void interpolateTimeDomainEField(double ***Exy, double ***Exz, double ***Eyx,
                                 double ***Eyz, double ***Ezx, double ***Ezy, 
                                 int i, int j, int k, int I, int J, int K, 
                                 double *Ex, double *Ey, double *Ez);