/* For now do not include any headers. Compile time errors will show which
 * are needed...
 */

// Pi...
const scalar pi = 3.141592653589793238; 

// Calculate angle between two vectors
scalar getAngle(vector& a, vector& b);

// Cotangent function to resemble algorihm notion from paper
scalar cot(scalar angle);

// Triangle area functions
// Use heron's formula
scalar areaHeron(vector& A, vector& B, vector& C);

// Use cross product property
scalar areaCross(vector& a, vector& b);

// Aux function to get a reference point for consistent normal vectors
// Does only work for convex, closed surfaces
point getRefPoint(const triSurface& front);
