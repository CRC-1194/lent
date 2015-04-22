/* For now do not include any headers. Compile time errors will show which
 * are needed...
 */

// Pi...
const scalar pi = Foam::constant::mathematical::pi; 

// Calculate angle between two vectors
scalar getAngle(vector& a, vector& b);

// Cotangent function to resemble algorihm notion from paper
scalar cot(vector a, vector b);

// Triangle area functions
// Use heron's formula
scalar areaHeron(vector& A, vector& B, vector& C);

// Use cross product property
scalar areaCross(vector& a, vector& b);

// Sort labels of triangle vertices in such a way, that the label passed
// as an argument is the first element of the returned list
labelList orderVertices(labelledTri& tri, label V);
