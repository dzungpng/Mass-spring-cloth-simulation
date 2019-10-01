#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <Eigen/Core>
#include <Eigen/Dense>

// template<class T, int dim>
class Geometry {
public:
    using T = float;
    int dim = 3;
    using TV = Eigen::Matrix<T,3,1>;
    T m_collision_stiffness;

    Geometry() : m_collision_stiffness(0) {}

    Geometry(const T collision_stiffness) : m_collision_stiffness(collision_stiffness) {}

    virtual TV pointCollision(const TV &point) = 0;
};

#endif // GEOMETRY_H


#ifndef SPHERE_H
#define SPHERE_H

#include "Geometry.h"

class Sphere : public Geometry {
public: 
    using T = float;
    using TV = Eigen::Matrix<T,3,1>;
    TV m_center; 
    T m_radius;

    Sphere(const T collision_stiffness, const TV center, const T radius) : 
        Geometry(collision_stiffness), m_center(center), m_radius(radius) {}

    Sphere() : m_center(TV::Zero()), m_radius(0) {}

    TV pointCollision(const TV &point) {
        /*
            * Using zero-length spring to calculate collision.
            */
        TV d = point - m_center;
        T distance = pow(d(0), 2) + pow(d(1), 2) + pow(d(2), 2);
        if (distance < m_radius*m_radius) {
            TV target = m_center + ((m_radius + 1e-5) / d.norm()) * (point - m_center);
            TV f_collision = - m_collision_stiffness * (point - target);
            return f_collision;  
        }
        return TV::Zero();
    }
};

#endif // SPHERE_H

class Plane : public Geometry {
public:

}
