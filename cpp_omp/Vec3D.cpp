// -*- adamgillard-cpp -*-
#include <iostream>
#include <cmath>
#include <random>

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

// Vec3D Inbuilt Functions --------------------------------------------
double Vec3D::EuclidDist(const Vec3D vec2) const
{
    /**
        Find the Euclidean/Cartesian distance between two 3D vectors.

        EQUATION
            [(x0 - x1)^2 + (y0 - y1)^2 + (z0 - z1)^2]^0.5
        
        RETURN
            double - Euclidean distance value
    */

    return sqrt(pow(x-vec2.x, 2) + pow(y-vec2.y, 2) + pow(z-vec2.z, 2));
}

Vec3D Vec3D::LimVec(const float limit)
{
    /**
        Normalise a vector to a given limit only if the magnitude of 
        the vector is larger than the limit.
        
        RETURN
            3D Vector - Normalised (or not) vector
    */

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
    /**
        Produce a 3D vector of the absolute (positive) values of an
        existing 3D vector.
        ie, negative values in a vector are made positive for all
        elements.
        
        RETURN
            3D Vector - elements are their corresponding absolute 
            values.
    */

    double x_ = abs(x);
    double y_ = abs(y);
    double z_ = abs(z);
    
    return Vec3D(x_, y_, z_);
}

double Vec3D::Mag()
{   
    /**
        Produce the magnitude of the 3D vector.
        
        EQUATION
            (x^2 + y^2 + z^2)^0.5

        RETURN
            double - value of vector magnitude.
    */

    double x_ = x;
    double y_ = y;
    double z_ = z;

    return sqrt(pow(x_, 2) + pow(y_, 2) + pow(z_, 2));
}


// Define vector arithmetic operators ---------------------------------
Vec3D operator * (const Vec3D vec, const float a)                   
{
    /**
        Vector-scalar multiplication
        
        RETURN
            3D Vector where each element is multiplied by a.
    */
    return Vec3D(a*vec.x, a*vec.y, a*vec.z);
}

Vec3D operator * (const float a, const Vec3D vec)
{
    /**
        Scalar-vector multiplication
        
        RETURN
            3D Vector where each element is multiplied by a.
    */
    return Vec3D(a*vec.x, a*vec.y, a*vec.z);
}

Vec3D operator * (const Vec3D vec1, const Vec3D vec2)
{
    /**
        Element-wise vector-vector multiplication
        
        RETURN
            3D Vector where corressponding elements from each input 
            vector are multipled.
    */
    return Vec3D(vec1.x * vec2.x, vec1.y * vec2.y, vec1.z * vec2.z);
}

Vec3D operator + (const Vec3D vec1, const Vec3D vec2)
{
    /**
        Element-wise vector-vector addition
        
        RETURN
            3D Vector where corressponding elements from each input 
            vector are added.
    */
    return Vec3D(vec1.x + vec2.x, vec1.y + vec2.y, vec1.z + vec2.z);
}

Vec3D operator - (const Vec3D vec1, const Vec3D vec2)
{
    /**
        Element-wise vector-vector subtraction
        
        RETURN
            3D Vector where elements from vec2 are subtracted from
            elements from vec1.
            
    */
    return Vec3D(vec1.x - vec2.x, vec1.y - vec2.y, vec1.z - vec2.z);
}

std::ostream& operator<<(std::ostream& os, const Vec3D v)
{
    /**
        Vector element readout for print and file-write
        
        RETURN
            std::ostream memory reference to sting readout of x, y, z
            values of the input vector.
    */
    return os << v.x << "," << v.y << "," << v.z;
}


// Vector-Vector Functions --------------------------------------------
double DotProd(Vec3D vec1, Vec3D vec2)
{
    /**
        Vector-vector dot product
        
        RETURN
            double of sum of element-wise vector multiplication.
    */
    return ((vec1.x * vec2.x) + (vec1.y * vec2.y) + (vec1.z * vec2.z));
}


// Random Vector generator --------------------------------------------
double RandVal(double min, double max)
{
    /**
        Generate a random value between min and max inputs using the
        <random> header. Uniform real distribution.
        
        RETURN
            double - random value between input values
    */

    std::uniform_real_distribution<double> unif(min, max);
    static std::default_random_engine val;

    return unif(val);
}

Vec3D RandVec(double min, double max)
{
    /**
        Generate a 3D vector of random values between the input min and
        max values. Uses __RandVal__ function to generate values.
        
        RETURN
            3D Vector of random values between input values
    */

    double x = RandVal(min, max);
    double y = RandVal(min, max);
    double z = RandVal(min, max);

    return Vec3D(x, y, z);
}