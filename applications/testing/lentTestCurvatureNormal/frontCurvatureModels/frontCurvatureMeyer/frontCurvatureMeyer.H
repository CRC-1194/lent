#ifndef frontCurvatureMeyer_H
#define frontCurvatureMeyer_H

#include "frontCurvatureModel"

class frontCurvatureMeyer
    :
    public frontCurvatureModel
{
    private
        void calcCurvatureNormal();

    public
        ~frontCurvatureMeyer();
};
