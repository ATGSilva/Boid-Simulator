#include <iostream>
#include <cmath>
#include <random>
#include "Vec3D.h"

// Define vector Euclidian Distance function
double Vec3D::EuclidDist(const Vec3D vec2) const
{
    return sqrt(pow(x-vec2.x, 2) + pow(y-vec2.y, 2) + pow(z-vec2.z, 2));
}

Vec3D Vec3D::LimVec(const short limit)
{
    double x_ = x;
    double y_ = y;
    double z_ = z;
    double norm_val = sqrt(pow(x_, 2) + pow(y_, 2) + pow(z_, 2));

    if (norm_val > limit)
    {
        double norm_fact = limit / norm_val;
        x_ = norm_fact * x_;
        y_ = norm_fact * y_;
        z_ = norm_fact * z_;
    }
    return Vec3D(x_, y_, z_);
}

Vec3D Vec3D::Abs()
{
    double x_ = abs(x);
    double y_ = abs(y);
    double z_ = abs(z);
    
    return Vec3D(x_, y_, z_);
}

// Define vector arithmetic operators
Vec3D operator * (const Vec3D vec, const float a)
{
    return Vec3D(a*vec.x, a*vec.y, a*vec.z);
}

Vec3D operator * (const float a, const Vec3D vec)
{
    return Vec3D(a*vec.x, a*vec.y, a*vec.z);
}

Vec3D operator * (const Vec3D vec1, const Vec3D vec2)
{
    return Vec3D(vec1.x * vec2.x, vec1.y * vec2.y, vec1.z * vec2.z);
}

Vec3D operator + (const Vec3D vec1, const Vec3D vec2)
{
    return Vec3D(vec1.x + vec2.x, vec1.y + vec2.y, vec1.z + vec2.z);
}

Vec3D operator - (const Vec3D vec1, const Vec3D vec2)
{
    return Vec3D(vec1.x - vec2.x, vec1.y - vec2.y, vec1.z - vec2.z);
}

bool operator != (const Vec3D vec1, const Vec3D vec2)
{
    return (vec1.x == vec2.x, vec1.y == vec2.y, vec1.z == vec2.z); 
}

std::ostream& operator<<(std::ostream& os, const Vec3D v)
{
    return os << v.x << "," << v.y << "," << v.z;
}


// Random Vector generator
double RandVal(double min, double max)
{
    std::uniform_real_distribution<double> unif(min, max);
    static std::default_random_engine val;

    return unif(val);
}

Vec3D RandVec(double min, double max)
{
    double x = RandVal(min, max);
    double y = RandVal(min, max);
    double z = RandVal(min, max);
    Vec3D random_vect = Vec3D(x, y, z);
    return random_vect;
}