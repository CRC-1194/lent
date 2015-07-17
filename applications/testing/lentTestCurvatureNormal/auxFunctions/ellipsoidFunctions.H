scalar computeU(vector& v);

scalar computeV(vector& v);

scalar ellipsoidDeviation(scalar u, scalar v, vector& axis, vector& p,
                          vector& center);

vector ellipsoidNormal(scalar u, scalar v, vector& axis);

scalar ellipsoidCurvature(scalar u, scalar v, vector& axis);
