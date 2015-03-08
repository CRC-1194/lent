/* Abstract class to define an interface to frontCurvatureModels
 */
#include "fvCFD.H"
#include "triSurface.H"
#include "triSurfaceFields.H"

class frontCurvatureModel
{
    private:
        // TODO: replace the copy of the front by a constant reference 
        // 'const triSurface& front;'
        triSurface front_;

        triSurfacePointVectorField curvatureNormals;

        virtual void calcCurvatureNormals();

    public:
    // Constructors
        frontCurvatureModel(const fileName& frontName);

    // Destructors
        virtual ~frontCurvatureModel();

    // Member functions
        const triSurfacePointVectorField& curvatureNormal() const;
        void curvature(triSurfacePointScalarField&) const;
        void normal(triSurfacePointVectorField&) const;
        void update();
};
