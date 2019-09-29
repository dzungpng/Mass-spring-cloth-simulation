#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <fstream>

template<class T, int dim>
class MassSpringSystem{
public:
    using TV = Eigen::Matrix<T,dim,1>;
    
    std::vector<Eigen::Matrix<int,2,1> > segments;
    std::vector<T> m;
    std::vector<TV> x;
    std::vector<TV> v;
    T youngs_modulus;
    T damping_coeff;
    std::vector<bool> node_is_fixed;
    std::vector<T> rest_length;

    MassSpringSystem()
    {}

    void evaluateSpringForces(std::vector<TV >& f)
    {
        f.clear();
        f.resize(x.size(), TV::Zero());
        
        for (int i = 0; i < segments.size(); i++) {
            int x1Idx = segments[i].row(0)(0);
            int x2Idx = segments[i].row(1)(0);
            TV x1 = x[x1Idx];
            TV x2 = x[x2Idx];
            T l = (x1 - x2).norm();
            TV n12 = (x1 - x2) / l;
            TV f1 = -youngs_modulus * (l/rest_length[i] - 1) * n12;
            // std::cout << f1.col(0)(0) << ", " << f1.col(0)(1) << ", " << f1.col(0)(2) << "\n";
            f[x1Idx] += f1;
            f[x2Idx] += -f1; 
        }
    }

    void evaluateDampingForces(std::vector<TV >& f)
    {
        f.clear();
        f.resize(x.size(), TV::Zero());

        for (int i = 0; i < segments.size(); i++) {
            int x1Idx = segments[i].row(0)(0);
            int x2Idx = segments[i].row(1)(0);
            TV x1 = x[x1Idx];
            TV x2 = x[x2Idx];
            T l = (x1 - x2).norm();
            TV n12 = (x1 - x2) / l;
            TV f1 = -damping_coeff * (n12 * n12.transpose()) * (v[x1Idx] - v[x2Idx]);
            f[x1Idx] += f1;
            f[x2Idx] += -f1; 
        }
    }

    void dumpPoly(std::string filename)
    {
        std::ofstream fs;
        fs.open(filename);
        fs << "POINTS\n";
        int count = 0;
        for (auto X : x) {
            fs << ++count << ":";
            for (int i = 0; i < dim; i++)
                fs << " " << X(i);
            if (dim == 2)
                fs << " 0";
            fs << "\n";
        }
        fs << "POLYS\n";
        count = 0;
        for (const Eigen::Matrix<int, 2, 1>& seg : segments)
            fs << ++count << ": " << seg(0) + 1 << " " << seg(1) + 1 << "\n"; // poly segment mesh is 1-indexed
        fs << "END\n";
        fs.close();
    }
    
};
