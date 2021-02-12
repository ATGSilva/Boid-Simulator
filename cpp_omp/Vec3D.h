/**
-*- adamgillard-cpp -*- Advanced Computational Physics -*-

Vec3D.h

Dependancy file to be used in BoidSimOMP.cpp and BoidSimOMPMPI.cpp.
Contains 3D Vector struct, modification functions and operator 
overrides.

FUNCTION SIGNATURE - RETURN TYPE
    Vec3D::EuclidDist(Vec3D) - double
    Vec3D::LimVec(const float) - Vec3D
    Vec3D::Abs() - Vec3D
    Vec3D::Mag() - double
    DotProd(Vec3D, Vec3D) - double
    RandVal(double, double) - double
    RandVec(double, double) - Vec3D
*/

#pragma once
#include <iostream>

struct Vec3D
{
    Vec3D() : x(0.), y(0.), z(0.) {}
    Vec3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {};

    double x, y, z;

    double EuclidDist(const Vec3D) const;
    Vec3D LimVec(const float limit);
    Vec3D Abs();
    double Mag();

    friend Vec3D operator*(const Vec3D, const float);
    friend Vec3D operator*(const float, const Vec3D);
    friend Vec3D operator*(const Vec3D, const Vec3D);
    friend Vec3D operator+(const Vec3D, const Vec3D);
    friend Vec3D operator-(const Vec3D, const Vec3D);
    friend std::ostream& operator<<(std::ostream&, const Vec3D);
};

double DotProd(Vec3D, Vec3D);
double RandVal(double, double);
Vec3D RandVec(double, double);