// -*- adamgillard-cpp -*-
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